import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Paths to your summary files
dirs = [
    ('Full', '/home/susan2000/Benchmark-final/PPI/Full/ligand/ligand_metrics_summary.csv'),
    ('Disorder 30', '/home/susan2000/Benchmark-final/PPI/Disorder/ligand_30/ligand_metrics_summary.csv'),
    ('Disorder 50', '/home/susan2000/Benchmark-final/PPI/Disorder/ligand_50/ligand_metrics_summary.csv')
]

metrics = ['RR', 'RP', 'LR', 'LP']
colors = ['#4C72B0', '#55A868', '#C44E52', '#8172B2']

fig, axes = plt.subplots(1, 3, figsize=(22, 8), sharey=True)

for idx, (title, summary_file) in enumerate(dirs):
    summary_df = pd.read_csv(summary_file)
    means = summary_df[[f'{m}_mean' for m in metrics]].values
    cis = summary_df[[f'{m}_CI' for m in metrics]].values
    models = summary_df['Model'].values
    x = np.arange(len(models))
    width = 0.22

    ax = axes[idx]
    for i, (metric, color) in enumerate(zip(metrics, colors)):
        bars = ax.barh(x + i*width, means[:, i], height=width, xerr=cis[:, i], capsize=8,
                       label=metric if idx == 0 else "", color=color, ecolor='black', linewidth=2, error_kw=dict(lw=2, capthick=2))
        for j, bar in enumerate(bars):
            width_val = bar.get_width()
            ci = cis[j, i]
            ax.annotate(f'{width_val:.3f}',
                        xy=(width_val + ci, bar.get_y() + bar.get_height() / 2),
                        xytext=(3, 0),
                        textcoords="offset points",
                        ha='left', va='center', fontsize=20)
    ax.set_yticks(x + 1.5*width)
    ax.set_yticklabels(models, fontsize=20)
    if idx == 0:
        ax.set_ylabel('Model', fontsize=20)
    else:
        ax.set_ylabel('')
    ax.set_xlabel('')
    # ax.set_title(title, fontsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.set_xlim([0, 1])

# Place the legend on the right top side, outside the last subplot
handles, labels = axes[-1].get_legend_handles_labels()
axes[-1].legend(handles, labels, fontsize=20, loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0.)

# Set a single x-axis label for the whole figure
fig.text(0.5, 0.04, 'Metric Value', ha='center', va='center', fontsize=20)

plt.tight_layout(rect=[0, 0.05, 0.95, 1])  # leave space for the legend
plt.savefig('/home/susan2000/Benchmark-final/PPI/ligand_metrics_overall_barplot.png', bbox_inches='tight')
plt.close()
print("Overall bar plot saved as ligand_metrics_overall_barplot.png")