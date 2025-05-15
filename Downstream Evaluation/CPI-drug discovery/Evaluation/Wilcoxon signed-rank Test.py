import os
import glob
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import rankdata, wilcoxon

# Convert p-values to -log10(p-value), handle zeros and NaNs
def safe_neglog10(x):
    if np.isnan(x):
        return np.nan
    elif x <= 0:
        return 300  # or use a large value like 50
    else:
        return -np.log10(x)

# Helper to compute and save p-value matrix and CSV for a given directory
def compute_and_save_pvalues(base_dir):
    gt_path = os.path.join(base_dir, "GPCR_test.json")
    def load_labels(json_path):
        with open(json_path, "r") as f:
            data = json.load(f)
        return {item["data"]: item["label"] for item in data}
    gt_dict = load_labels(gt_path)
    all_keys = sorted(gt_dict.keys())
    gt_labels = np.array([gt_dict[k] for k in all_keys])
    gt_ranks = rankdata(gt_labels, method='average')
    model_files = [f for f in glob.glob(os.path.join(base_dir, "GPCR_test_*.json")) if not f.endswith("GPCR_test.json")]
    model_names = [os.path.splitext(os.path.basename(f))[0].replace("GPCR_test_", "") for f in model_files]
    rank_errors = {}
    for model_file, model_name in zip(model_files, model_names):
        pred_dict = load_labels(model_file)
        pred_labels = np.array([pred_dict[k] for k in all_keys])
        pred_ranks = rankdata(pred_labels, method='average')
        rank_errors[model_name] = np.abs(pred_ranks - gt_ranks)
    n_models = len(model_names)
    pval_matrix = np.zeros((n_models, n_models))
    for i in range(n_models):
        for j in range(n_models):
            if i == j:
                pval_matrix[i, j] = np.nan
            else:
                stat, pval = wilcoxon(rank_errors[model_names[i]], rank_errors[model_names[j]])
                pval_matrix[i, j] = pval
    # Save p-values to CSV for heatmap plotting
    rows = []
    for i, m1 in enumerate(model_names):
        for j, m2 in enumerate(model_names):
            if i < j:
                rows.append({'Model 1': m1, 'Model 2': m2, 'p-value': pval_matrix[i, j]})
    pval_df = pd.DataFrame(rows)
    pval_df.to_csv(os.path.join(base_dir, 'model_comparison_results.csv'), index=False)
    # Plot single heatmap for this dataset (log scale)
    neglog_pval_matrix = np.vectorize(safe_neglog10)(pval_matrix)
    plt.figure(figsize=(10, 8))
    sns.heatmap(neglog_pval_matrix, annot=True, fmt=".2f", xticklabels=model_names, yticklabels=model_names, cmap="viridis", cbar_kws={'label': '-log10(p-value)'})
    plt.title(f"Pairwise Wilcoxon Test -log10(p-value) (Rank Errors)\n{os.path.basename(base_dir)}")
    plt.tight_layout()
    plt.savefig(os.path.join(base_dir, "rank_significance_heatmap.png"))
    plt.close()
    print(f"Heatmap and CSV saved in {base_dir}")
    return model_names

# Main logic
if __name__ == "__main__":
    cwd = os.path.abspath(os.getcwd())
    # If run from the Drug directory, plot the overall heatmap for all three datasets
    if os.path.basename(cwd) == "Drug":
        # Compute p-values for all three datasets
        subdirs = [
            ("Full", os.path.join(cwd, "full/test")),
            ("Disorder 30", os.path.join(cwd, "disorder_30/test")),
            ("Disorder 50", os.path.join(cwd, "disorder_50/test")),
        ]
        all_model_names = None
        for _, d in subdirs:
            model_names = compute_and_save_pvalues(d)
            if all_model_names is None:
                all_model_names = model_names
        # Now plot the overall heatmap as in PPI/heatmap.py (log scale)
        pval_files = [
            (name, os.path.join(path, 'model_comparison_results.csv'))
            for name, path in subdirs
        ]
        matrices = []
        vmin, vmax = None, None
        for _, file in pval_files:
            df = pd.read_csv(file)
            pval_matrix = pd.DataFrame(np.nan, index=all_model_names, columns=all_model_names)
            for _, row in df.iterrows():
                m1, m2, p = row['Model 1'], row['Model 2'], row['p-value']
                p = max(p, 1e-300)
                logp = safe_neglog10(p)
                pval_matrix.loc[m1, m2] = logp
                pval_matrix.loc[m2, m1] = logp
            matrices.append(pval_matrix)
            min_p = np.nanmin(pval_matrix.values)
            max_p = np.nanmax(pval_matrix.values)
            vmin = min_p if vmin is None else min(vmin, min_p)
            vmax = max_p if vmax is None else max(vmax, max_p)
        fig, axes = plt.subplots(1, 3, figsize=(24, 8), sharey=True)
        titles = [x[0] for x in pval_files]
        cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
        for i, (ax, matrix, title) in enumerate(zip(axes, matrices, titles)):
            sns.heatmap(matrix, annot=True, fmt='.1f', cmap='viridis', vmin=vmin, vmax=vmax,
                        cbar=i==2, cbar_ax=cbar_ax if i==2 else None, ax=ax,
                        linewidths=0.5, linecolor='white', annot_kws={'fontsize': 12})
            ax.set_xlabel('')
            if i == 0:
                ax.set_ylabel('Model', fontsize=23)
            else:
                ax.set_ylabel('')
            ax.set_xticklabels(matrix.columns, fontsize=20, rotation=45, ha='right')
            ax.set_yticklabels(matrix.index, fontsize=20, rotation=45)
        cbar_ax.set_ylabel('-log10(p-value)', fontsize=20)
        cbar_ax.tick_params(labelsize=20)
        plt.tight_layout(rect=[0, 0, 0.9, 1])
        plt.savefig(os.path.join(cwd, 'overall_pvalue_heatmap.png'), bbox_inches='tight')
        plt.close()
        print("Overall p-value heatmap with shared colorbar saved as overall_pvalue_heatmap.png")
    else:
        # Just compute and save for this dataset
        compute_and_save_pvalues(cwd)
