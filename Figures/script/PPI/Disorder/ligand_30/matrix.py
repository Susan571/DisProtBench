import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

src_dir = '/home/susan2000/Benchmark-final/PPI/Disorder/ligand_30'
summary_file = os.path.join(src_dir, 'ligand_metrics_summary.csv')

# Step 1: Calculate means, variances, and 95% CI for each model
summary = []
metrics = ['RR', 'RP', 'LR', 'LP']
for file in glob.glob(os.path.join(src_dir, '*.csv')):
    if file.endswith('ligand_metrics_summary.csv'):
        continue
    model = os.path.splitext(os.path.basename(file))[0]
    df = pd.read_csv(file)
    if not all(col in df.columns for col in metrics):
        continue
    stats = {'Model': model}
    n = len(df)
    for col in metrics:
        mean = df[col].mean()
        var = df[col].var()
        ci = 1.96 * np.sqrt(var / n) if n > 0 else np.nan
        stats[f'{col}_mean'] = round(mean, 4)
        stats[f'{col}_var'] = round(var, 6)
        stats[f'{col}_CI'] = round(ci, 6)
    summary.append(stats)

summary_df = pd.DataFrame(summary)
summary_df.to_csv(summary_file, index=False)
print("Summary saved to ligand_metrics_summary.csv")

# Step 2: Plotting
means = summary_df[[f'{m}_mean' for m in metrics]].values
vars_ = summary_df[[f'{m}_var' for m in metrics]].values
models = summary_df['Model'].values

# Get sample size for each model
sample_sizes = {}
for model in models:
    file = os.path.join(src_dir, model + '.csv')
    df = pd.read_csv(file)
    sample_sizes[model] = len(df)

# Calculate 95% CI for each metric/model (for plotting)
cis = []
for i, model in enumerate(models):
    n = sample_sizes[model]
    ci = 1.96 * np.sqrt(vars_[i] / n)
    cis.append(ci)
cis = np.array(cis)

x = np.arange(len(models))
width = 0.18

fig, ax = plt.subplots(figsize=(12, 8))
colors = ['#4C72B0', '#55A868', '#C44E52', '#8172B2']

# Plot horizontal bars with error bars
for i, (metric, color) in enumerate(zip(metrics, colors)):
    bars = ax.barh(x + i*width, means[:, i], height=width, xerr=cis[:, i], capsize=8, 
                  label=metric, color=color, ecolor='black', linewidth=2, error_kw=dict(lw=2, capthick=2))
    # Add value labels to the right of the error bars
    for j, bar in enumerate(bars):
        width_val = bar.get_width()
        ci = cis[j, i]
        ax.annotate(f'{width_val:.3f}',
                    xy=(width_val + ci, bar.get_y() + bar.get_height() / 2),
                    xytext=(3, 0),  # 3 points horizontal offset to the right of the error bar
                    textcoords="offset points",
                    ha='left', va='center', fontsize=15)

# Set x-axis limits to zoom in
xmin = max(0, means.min() - 0.05)
xmax = min(1, means.max() + 0.08)
ax.set_xlim([xmin, xmax])

ax.set_yticks(x + 1.5*width)
ax.set_yticklabels(models, fontsize=20)
ax.set_xlabel('Metric Value', fontsize=20)
ax.set_ylabel('Model', fontsize=20)
#ax.set_title('Model Metrics with 95% CI', fontsize=15)
ax.legend(fontsize=18)
# Set x-axis tick label font size
ax.tick_params(axis='x', labelsize=20)
plt.tight_layout()
plot_file = os.path.join(src_dir, 'ligand_metrics_barplot_2.png')
plt.savefig(plot_file)
plt.close()
print("Horizontal bar plot with 95% CI saved as ligand_metrics_barplot.png")