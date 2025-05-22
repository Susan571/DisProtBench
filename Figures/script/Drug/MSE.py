import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# File paths and dataset labels
files = [
    ('Full', '/home/susan2000/Benchmark-final/Drug/full/train/MSE_results.csv'),
    ('Disorder_30', '/home/susan2000/Benchmark-final/Drug/disorder_30/train/MSE_results.csv'),
    ('Disorder_50', '/home/susan2000/Benchmark-final/Drug/disorder_50/train/MSE_results.csv')
]

# Read and merge data
dfs = []
for label, path in files:
    df = pd.read_csv(path)
    df['Dataset'] = label
    dfs.append(df)
all_df = pd.concat(dfs)

# Pivot for grouped bar plot
pivot = all_df.pivot(index='Model', columns='Dataset', values='MSE')
pivot = pivot.reindex(sorted(pivot.index))  # Sort models alphabetically

# Plot
fig, ax = plt.subplots(figsize=(27, 10))
bar_width = 0.3  # Wider bars
x = np.arange(len(pivot.index))
colors = ['#4C72B0', '#55A868', '#C44E52']

for i, dataset in enumerate(pivot.columns):
    bars = ax.bar(x + i*bar_width, pivot[dataset], width=bar_width, label=dataset, color=colors[i])
    # Add value labels above each bar
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.4f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 15),  # 15 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=13)

ax.set_xticks(x + bar_width)
ax.set_xticklabels(pivot.index, fontsize=20)
ax.set_ylabel('MSE', fontsize=20)
ax.set_xlabel('Model', fontsize=20)
ax.legend(fontsize=20)
#ax.set_title('MSE by Model and Dataset', fontsize=20)
ax.tick_params(axis='y', labelsize=20)  # Make y-axis tick labels larger
plt.tight_layout()
plt.savefig('/home/susan2000/Benchmark-final/Drug/MSE_barplot.png', bbox_inches='tight')
plt.close()
print("MSE bar plot saved as MSE_barplot.png")