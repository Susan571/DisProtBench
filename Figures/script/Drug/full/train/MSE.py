import json
import glob
import os
import numpy as np
from sklearn.metrics import mean_squared_error
import pandas as pd

# Paths
base_dir = "/home/susan2000/Benchmark-final/Drug/full/train"
gt_path = os.path.join(base_dir, "GPCR_train.json")

# Load ground truth
with open(gt_path, "r") as f:
    ground_truth = json.load(f)
gt_dict = {item["data"]: item["label"] for item in ground_truth}

# Find all model files (exclude ground truth)
model_files = [f for f in glob.glob(os.path.join(base_dir, "GPCR_train_*.json")) if not f.endswith("GPCR_train.json")]

results = []
for model_file in model_files:
    model_name = os.path.splitext(os.path.basename(model_file))[0].replace("GPCR_train_", "")
    with open(model_file, "r") as f:
        predictions = json.load(f)
    pred_dict = {item["data"]: item["label"] for item in predictions}
    common_keys = set(gt_dict.keys()) & set(pred_dict.keys())
    y_true = np.array([gt_dict[k] for k in common_keys])
    y_pred = np.array([pred_dict[k] for k in common_keys])
    mse = mean_squared_error(y_true, y_pred)
    var = np.var((y_true - y_pred) ** 2, ddof=1)
    n = len(y_true)
    ci = 1.96 * np.sqrt(var / n) if n > 0 else np.nan
    results.append({
        "Model": model_name,
        "MSE": round(mse, 4),
        "Variance": round(var, 4),
        "95%CI": round(ci, 4)
    })

results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(base_dir, "MSE_results.csv"), index=False)
print("Results saved to MSE_results.csv")
