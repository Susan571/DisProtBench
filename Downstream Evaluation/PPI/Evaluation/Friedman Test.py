import pandas as pd
import numpy as np
import glob
from statsmodels.stats.contingency_tables import mcnemar
from itertools import combinations
from sklearn.metrics import accuracy_score
import os
import seaborn as sns
import matplotlib.pyplot as plt

# 1. Load true labels
true_labels = pd.read_csv('true_label.csv')
print(true_labels.columns)
true_labels['key'] = true_labels['Seq_A'] + '_' + true_labels['Seq_B']
true_labels = true_labels.set_index('key')

# 2. Load all model files (excluding true_label.csv)
model_files = [f for f in glob.glob("*.csv") if f != 'true_label.csv']
model_data = {}
for file in model_files:
    df = pd.read_csv(file)
    # Check for required columns
    required_cols = {'Seq_A', 'Seq_B', 'Label'}
    if not required_cols.issubset(df.columns):
        print(f"Skipping {file}: missing required columns {required_cols - set(df.columns)}")
        continue
    model_name = os.path.splitext(os.path.basename(file))[0]
    df['key'] = df['Seq_A'] + '_' + df['Seq_B']
    # Merge with true labels
    merged = df.set_index('key').join(true_labels[['Label']], rsuffix='_true', how='inner')
    merged = merged.rename(columns={'Label': 'Pred', 'Label_true': 'True'})
    model_data[model_name] = merged

# 3. Prepare results storage
results = []

# 4. Pairwise comparison
for (model1, data1), (model2, data2) in combinations(model_data.items(), 2):
    # Find intersection of keys
    common_keys = data1.index.intersection(data2.index)
    if len(common_keys) == 0:
        continue  # No common samples to compare

    # Get predictions and true labels
    y_true = data1.loc[common_keys, 'True']
    y_pred1 = data1.loc[common_keys, 'Pred']
    y_pred2 = data2.loc[common_keys, 'Pred']

    # Map to binary
    y_pred1_bin = y_pred1.map({'positive': 1, 'negative': 0})
    y_pred2_bin = y_pred2.map({'positive': 1, 'negative': 0})
    y_true_bin = y_true.map({'positive': 1, 'negative': 0})

    # Calculate accuracy for each model
    acc1 = accuracy_score(y_true_bin, y_pred1_bin)
    acc2 = accuracy_score(y_true_bin, y_pred2_bin)

    # McNemar's test
    both_correct = np.sum((y_pred1_bin == y_true_bin) & (y_pred2_bin == y_true_bin))
    model1_only = np.sum((y_pred1_bin == y_true_bin) & (y_pred2_bin != y_true_bin))
    model2_only = np.sum((y_pred1_bin != y_true_bin) & (y_pred2_bin == y_true_bin))
    both_wrong = np.sum((y_pred1_bin != y_true_bin) & (y_pred2_bin != y_true_bin))
    table = [[both_correct, model1_only], [model2_only, both_wrong]]

    result = mcnemar(table, exact=True)
    p_value = result.pvalue

    results.append({
        'Model 1': model1,
        'Model 2': model2,
        'N Common': len(common_keys),
        'Acc 1': acc1,
        'Acc 2': acc2,
        'p-value': p_value,
        'Significant (p<0.05)': p_value < 0.05
    })

# 5. Export results
results_df = pd.DataFrame(results)
results_df.to_csv('model_comparison_results.csv', index=False)

print("Results exported to model_comparison_results.csv and model_comparison_results.tex")

# Plot heatmap of p-values
# Pivot the results to a matrix form
pivot = results_df.pivot(index='Model 1', columns='Model 2', values='p-value')
# Make the matrix symmetric
all_models = sorted(set(results_df['Model 1']).union(results_df['Model 2']))
pval_matrix = pd.DataFrame(np.nan, index=all_models, columns=all_models)
for _, row in results_df.iterrows():
    m1, m2, p = row['Model 1'], row['Model 2'], row['p-value']
    pval_matrix.loc[m1, m2] = p
    pval_matrix.loc[m2, m1] = p

# Set font sizes
plt.rcParams.update({'font.size': 15})  # Increase base font size
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18

plt.figure(figsize=(10, 8))
sns.heatmap(pval_matrix, annot=True, fmt='.2g', cmap='viridis', cbar_kws={'label': 'p-value'})
#plt.title('Heatmap of p-values (McNemar test)')
plt.tight_layout()
plt.savefig('pvalue_heatmap.png')
plt.close()
print('Heatmap saved as pvalue_heatmap.png')