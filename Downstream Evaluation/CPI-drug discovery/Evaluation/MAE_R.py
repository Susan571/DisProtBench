import json
import glob
import os
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error
from scipy.stats import pearsonr

# Paths
base_dir = "/Benchmark-final/Drug/full/test"
gt_path = os.path.join(base_dir, "GPCR_test.json")

# Load ground truth
with open(gt_path, "r") as f:
    ground_truth = json.load(f)
gt_dict = {item["data"]: item["label"] for item in ground_truth}

# Find all model files (exclude ground truth)
model_files = [f for f in glob.glob(os.path.join(base_dir, "GPCR_test_*.json")) if not f.endswith("GPCR_test.json")]

results = []
for model_file in model_files:
    model_name = os.path.splitext(os.path.basename(model_file))[0].replace("GPCR_test_", "")
    with open(model_file, "r") as f:
        predictions = json.load(f)
    pred_dict = {item["data"]: item["label"] for item in predictions}
    common_keys = set(gt_dict.keys()) & set(pred_dict.keys())
    y_true = np.array([gt_dict[k] for k in common_keys])
    y_pred = np.array([pred_dict[k] for k in common_keys])
    n = len(y_true)

    # MAE and its variance/CI
    abs_err = np.abs(y_true - y_pred)
    mae = abs_err.mean()
    mae_var = abs_err.var(ddof=1)
    mae_ci = 1.96 * np.sqrt(mae_var / n) if n > 0 else np.nan

    # Pearson R and its variance/CI (Fisher transformation)
    r, _ = pearsonr(y_true, y_pred)
    if n > 3:
        fisher_z = np.arctanh(r)
        se = 1 / np.sqrt(n - 3)
        z_ci = 1.96 * se
        r_ci_low = np.tanh(fisher_z - z_ci)
        r_ci_high = np.tanh(fisher_z + z_ci)
        r_var = se**2
        r_ci = (r_ci_high - r_ci_low) / 2
    else:
        r_var = np.nan
        r_ci = np.nan

    results.append({
        "Model": model_name,
        "MAE": round(mae, 6),
        "MAE_Var": round(mae_var, 6),
        "MAE_CI": round(mae_ci, 6),
        "R": round(r, 6),
        "R_Var": round(r_var, 6),
        "R_CI": round(r_ci, 6)
    })

results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(base_dir, "Test_results.csv"), index=False)
print("Results saved to Test_results.csv")
