import pandas as pd
import random

def limit_partners(df, max_partners=10):
    """
    Limit the number of partners for each gene to max_partners,
    prioritizing partners with pathogenic variants (appearance > 0)
    """
    # Create a copy of the dataframe to avoid modifying the original
    filtered_df = pd.DataFrame()
    
    # Get unique genes from both Gene_A and Gene_B
    all_genes = set(df['Gene_A'].unique()) | set(df['Gene_B'].unique())
    
    for gene in all_genes:
        # Get all interactions where this gene appears in either Gene_A or Gene_B
        # Sort by appearance to prioritize interactions with higher appearance values
        gene_interactions = df[
            ((df['Gene_A'] == gene) | (df['Gene_B'] == gene))
        ].sort_values('total_variants', ascending=False)
        
        if len(gene_interactions) <= max_partners:
            filtered_df = pd.concat([filtered_df, gene_interactions])
        else:
            # Take the top max_partners interactions (already sorted by appearance)
            selected = gene_interactions.head(max_partners)
            filtered_df = pd.concat([filtered_df, selected])
    
    # Remove duplicates (since each interaction might have been selected from both genes' perspective)
    filtered_df = filtered_df.drop_duplicates()
    
    return filtered_df

try:
    # Read your sorted interactions file
    df = pd.read_excel('sorted_interactions_variants.xlsx')
    
    # Apply the partner limitation
    filtered_df = limit_partners(df, max_partners=10)
    
    # Save the results to Excel
    filtered_df.to_excel('limited_partners_interactions_variants.xlsx', index=False)
    
    # Print summary statistics
    print("\nSummary after limiting partners:")
    print(f"Original number of interactions: {len(df)}")
    print(f"Number of interactions after limiting partners: {len(filtered_df)}")
    print(f"Number of unique genes: {len(set(filtered_df['Gene_A']) | set(filtered_df['Gene_B']))}")
    
    # Statistics about appearance values
    print("\nAppearance statistics in filtered dataset:")
    print(f"Mean appearance: {filtered_df['total_variants'].mean():.2f}")
    print(f"Median appearance: {filtered_df['total_variants'].median():.2f}")
    print(f"Max appearance: {filtered_df['total_variants'].max()}")
    print(f"Min appearance: {filtered_df['total_variants'].min()}")

except ModuleNotFoundError:
    print("The openpyxl module is required to work with Excel files.")
    print("Please install it using:")
    print("pip install openpyxl")
    print("or")
    print("conda install openpyxl")
except Exception as e:
    print(f"An error occurred: {str(e)}")