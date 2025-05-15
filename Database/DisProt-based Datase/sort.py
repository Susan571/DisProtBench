import pandas as pd

# Read TSV file with tab separator
df = pd.read_csv('/home/xyzeng/Data/Uniprot/filtered_interactions_with_variant_counts.tsv', sep='\t')  # Changed separator to '\t'

# Print column names to verify
print("Column names in the dataframe:")
print(df.columns.tolist())

# Sort by 'appearance' column in descending order
df_sorted = df.sort_values(by='total_variants', ascending=False)

# Save to Excel file
df_sorted.to_excel('sorted_interactions_variants.xlsx', index=False)

# Display the first few rows to confirm sorting
print("\nSorted data (top rows):")
print(df_sorted[['Uniprot_A', 'Uniprot_B', 'Gene_A', 'Gene_B', 'total_variants']].head())