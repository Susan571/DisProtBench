import pandas as pd
import glob
import os
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

# Load true labels
true_labels = pd.read_csv('true_label.csv')
true_labels['key'] = true_labels['Seq_A'] + '_' + true_labels['Seq_B']
true_labels = true_labels.set_index('key')

# Find all model files (excluding true_label.csv)
model_files = [f for f in glob.glob("*.csv") if f != 'true_label.csv']

results = []

for file in model_files:
    model_name = os.path.splitext(os.path.basename(file))[0]
    df = pd.read_csv(file)
    if not {'Seq_A', 'Seq_B', 'Label'}.issubset(df.columns):
        print(f"Skipping {file}: missing required columns")
        continue
    df['key'] = df['Seq_A'] + '_' + df['Seq_B']
    merged = df.set_index('key').join(true_labels[['Label']], rsuffix='_true', how='inner')
    merged = merged.rename(columns={'Label': 'Pred', 'Label_true': 'True'})
    # Map to binary
    y_true = merged['True'].map({'positive': 1, 'negative': 0})
    y_pred = merged['Pred'].map({'positive': 1, 'negative': 0})
    # Calculate metrics
    acc = accuracy_score(y_true, y_pred)
    prec = precision_score(y_true, y_pred, zero_division=0)
    rec = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    results.append({
        'Model': model_name,
        'N': len(merged),
        'Accuracy': acc,
        'Precision': prec,
        'Recall': rec,
        'F1 Score': f1
    })

# Save results
results_df = pd.DataFrame(results)
results_df.to_csv('model_metrics.csv', index=False)
print("Results saved to model_metrics.csv")
print(results_df)