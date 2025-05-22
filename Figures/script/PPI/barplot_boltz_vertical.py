import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Paths to your summary files
summary_files = [
    ('Full', 'Full/ligand/ligand_metrics_summary.csv'),
    ('Disorder 30', 'Disorder/ligand_30/ligand_metrics_summary.csv'),
    ('Disorder 50', 'Disorder/ligand_50/ligand_metrics_summary.csv')
]

metrics = ['RR', 'RP', 'LR', 'LP']
colors = ['#335372', '#4F7BA8', '#82A4C8', '#CEDCE9']

# Collect means and CIs for Boltz across datasets
means = []
cis = []
for title, file in summary_files:
    df = pd.read_csv(file)
    boltz_row = df[df['Model'] == 'Boltz']
    means.append([boltz_row[f'{m}_mean'].values[0] for m in metrics])
    cis.append([boltz_row[f'{m}_CI'].values[0] for m in metrics])
means = np.array(means)  # shape: (3, 4)
cis = np.array(cis)      # shape: (3, 4)

x = np.arange(len(summary_files))  # 0, 1, 2
width = 0.18

fig, ax = plt.subplots(figsize=(10, 8))

for i, (metric, color) in enumerate(zip(metrics, colors)):
    bars = ax.bar(x + i*width - 1.5*width, means[:, i], width=width, yerr=cis[:, i], capsize=8,
                  label=metric, color=color, ecolor='black', linewidth=2, error_kw=dict(lw=2, capthick=2))

ax.set_xticks(x)
ax.set_xticklabels([
    'Full',
    r'pLDDT $\geq$ 50',
    r'pLDDT $\geq$ 50'
], fontsize=20)
ax.set_ylabel('Metric Value', fontsize=20)
#ax.set_xlabel('Dataset', fontsize=18)
ax.set_ylim([0.6, 0.85])
ax.legend(fontsize=20, ncol=2)
ax.tick_params(axis='y', labelsize=20)
#plt.title('Boltz Metrics Across Datasets', fontsize=20)
plt.tight_layout()
plt.savefig('boltz_metrics_vertical_barplot.png', bbox_inches='tight')
plt.close()
print("Vertical bar plot for Boltz saved as boltz_metrics_vertical_barplot.png") 