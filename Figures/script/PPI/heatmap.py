import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Paths to your p-value CSVs (now in p_value folders)
pval_files = [
    ('Full', '/home/susan2000/Benchmark-final/PPI/Full/p_value/model_comparison_results.csv'),
    ('Disorder 30', '/home/susan2000/Benchmark-final/PPI/Disorder/p_value_30/model_comparison_results.csv'),
    ('Disorder 50', '/home/susan2000/Benchmark-final/PPI/Disorder/p_value_50/model_comparison_results.csv')
]

# Get all unique models from the first file (assume all have the same models)
df0 = pd.read_csv(pval_files[0][1])
all_models = sorted(set(df0['Model 1']).union(df0['Model 2']))

# Prepare p-value matrices for each file
matrices = []
vmin, vmax = None, None
for _, file in pval_files:
    df = pd.read_csv(file)
    pval_matrix = pd.DataFrame(np.nan, index=all_models, columns=all_models)
    for _, row in df.iterrows():
        m1, m2, p = row['Model 1'], row['Model 2'], row['p-value']
        # Avoid log(0) by setting a minimum value
        p = max(p, 1e-300)
        logp = -np.log10(p)
        pval_matrix.loc[m1, m2] = logp
        pval_matrix.loc[m2, m1] = logp
    matrices.append(pval_matrix)
    # Track global min/max for colorbar
    min_p = np.nanmin(pval_matrix.values)
    max_p = np.nanmax(pval_matrix.values)
    vmin = min_p if vmin is None else min(vmin, min_p)
    vmax = max_p if vmax is None else max(vmax, max_p)

# Plot
fig, axes = plt.subplots(1, 3, figsize=(24, 8), sharey=True)
titles = [x[0] for x in pval_files]

# Create a single colorbar axis
cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])

for i, (ax, matrix, title) in enumerate(zip(axes, matrices, titles)):
    sns.heatmap(matrix, annot=True, fmt='.2g', cmap='Blues', vmin=vmin, vmax=vmax,
                cbar=i==2, cbar_ax=cbar_ax if i==2 else None, ax=ax,
                linewidths=0.5, linecolor='white', annot_kws={'fontsize': 13})
    #ax.set_title(title, fontsize=20)
    ax.set_xlabel('')
    if i == 0:
        ax.set_ylabel('Model', fontsize=23)
    else:
        ax.set_ylabel('')
    ax.set_xticklabels(matrix.columns, fontsize=20, rotation=45, ha='right')
    ax.set_yticklabels(matrix.index, fontsize=20, rotation=45)

# Only the last heatmap gets the colorbar
cbar_ax.set_ylabel('-log10(p-value)', fontsize=20)
cbar_ax.tick_params(labelsize=20)

plt.tight_layout(rect=[0, 0, 0.9, 1])
plt.savefig('/home/susan2000/Benchmark-final/PPI/overall_pvalue_heatmap.png', bbox_inches='tight')
plt.close()
print("Overall p-value heatmap with shared colorbar saved as overall_pvalue_heatmap.png")
