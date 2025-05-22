import os
import glob
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import rankdata, wilcoxon

# Paths
base_dir = os.path.dirname(os.path.abspath(__file__))
gt_path = os.path.join(base_dir, "GPCR_test.json")

# Load ground truth
def load_labels(json_path):
    with open(json_path, "r") as f:
        data = json.load(f)
    return {item["data"]: item["label"] for item in data}

gt_dict = load_labels(gt_path)
all_keys = sorted(gt_dict.keys())
gt_labels = np.array([gt_dict[k] for k in all_keys])
gt_ranks = rankdata(gt_labels, method='average')

# Find all model files (exclude ground truth)
model_files = [f for f in glob.glob(os.path.join(base_dir, "GPCR_test_*.json")) if not f.endswith("GPCR_test.json")]
model_names = [os.path.splitext(os.path.basename(f))[0].replace("GPCR_test_", "") for f in model_files]

# Load model predictions and compute rank errors
rank_errors = {}
for model_file, model_name in zip(model_files, model_names):
    pred_dict = load_labels(model_file)
    # Ensure same order as ground truth
    pred_labels = np.array([pred_dict[k] for k in all_keys])
    pred_ranks = rankdata(pred_labels, method='average')
    rank_errors[model_name] = np.abs(pred_ranks - gt_ranks)

# Pairwise Wilcoxon signed-rank test
n_models = len(model_names)
pval_matrix = np.zeros((n_models, n_models))
for i in range(n_models):
    for j in range(n_models):
        if i == j:
            pval_matrix[i, j] = np.nan
        else:
            stat, pval = wilcoxon(rank_errors[model_names[i]], rank_errors[model_names[j]])
            pval_matrix[i, j] = pval

# Plot heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(pval_matrix, annot=True, fmt=".2e", xticklabels=model_names, yticklabels=model_names, cmap="viridis", cbar_kws={'label': 'p-value'})
plt.title("Pairwise Wilcoxon Test p-values (Rank Errors)")
plt.tight_layout()
plt.savefig(os.path.join(base_dir, "rank_significance_heatmap_try.png"))
plt.show()
print("Heatmap saved as rank_significance_heatmap_try.png")
