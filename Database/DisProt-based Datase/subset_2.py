import pandas as pd

def load_filtered_summary_genes(summary_file):
    """Load genes from filtered summary data"""
    summary_data = pd.read_csv(summary_file, sep='\t')
    return set(summary_data['SubmittedGeneSymbol'].unique())

def filter_interactions(interaction_file, gene_set):
    """Filter interactions where either Gene_A or Gene_B is in the gene set"""
    # Load interaction data
    interactions = pd.read_csv(interaction_file, sep='\t')
    
    # Filter rows where either Gene_A or Gene_B is in the gene set
    filtered_interactions = interactions[
        (interactions['Gene_A'].isin(gene_set)) |
        (interactions['Gene_B'].isin(gene_set))
    ]
    
    return filtered_interactions

def main():
    summary_file = '/home/xyzeng/Data/Uniprot/filtered_summary_data_2.txt'
    interaction_file = '/home/xyzeng/Data/Uniprot/selected_interactions_20_unique.txt'
    output_file = '/home/xyzeng/Data/Uniprot/final_seq_unique.csv'
    
    print("Loading summary genes...")
    summary_genes = load_filtered_summary_genes(summary_file)
    print(f"Found {len(summary_genes)} unique genes in summary data")
    
    print("\nFiltering interactions...")
    interactions = pd.read_csv(interaction_file, sep='\t')
    filtered_interactions = filter_interactions(interaction_file, summary_genes)
    
    print("\nFiltering Statistics:")
    print(f"Original interaction rows: {len(interactions)}")
    print(f"Filtered interaction rows: {len(filtered_interactions)}")
    print(f"Rows removed: {len(interactions) - len(filtered_interactions)}")
    
    # Calculate how many unique genes from each source made it to the final dataset
    final_genes_a = set(filtered_interactions['Gene_A'])
    final_genes_b = set(filtered_interactions['Gene_B'])
    final_genes_all = final_genes_a | final_genes_b
    
    print(f"\nUnique genes in filtered interactions: {len(final_genes_all)}")
    print(f"Unique genes from Gene_A: {len(final_genes_a)}")
    print(f"Unique genes from Gene_B: {len(final_genes_b)}")
    
    # Save filtered interactions
    print(f"\nSaving filtered interactions to {output_file}")
    filtered_interactions.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()