import pandas as pd

def count_clinical_significance_rows(gene, filtered_df):
    """
    Count number of rows where the gene appears in the already filtered dataset
    """
    return len(filtered_df[filtered_df['geneId'].apply(lambda x: gene in str(x).split('|') if pd.notna(x) else False)])

def process_variants_and_interactions(variant_file, interaction_file):
    """
    Process variant and interaction files based on specified criteria
    
    Args:
        variant_file (str): Path to the variant TSV file
        interaction_file (str): Path to the interaction TXT file
        
    Returns:
        tuple: (filtered interactions DataFrame, original interaction count, filtered interaction count)
    """
    # Read the variant TSV file
    variants_df = pd.read_csv(variant_file, sep='\t')
    
    # Step 1: Filter by clinical significance (P or LP)
    filtered_df = variants_df[variants_df['_clinSignCode'].isin(['P'])]
    
    # Step 2: Filter for missense variants
    filtered_df = filtered_df[filtered_df['molConseq'] == 'missense variant']
    
    # Step 3: Keep only unique variants based on vcfDesc
    filtered_df = filtered_df.drop_duplicates(subset=['vcfDesc'])
    
    # Extract all gene IDs from the filtered variants
    gene_set = set()
    for genes in filtered_df['geneId']:
        if pd.notna(genes):  # Check if the value is not NaN
            gene_set.update(genes.split('|'))
    
    # Read interaction file
    interactions_df = pd.read_csv(interaction_file, sep='\t')
    total_rows = len(interactions_df)
    
    # Remove rows where Uniprot_A equals Uniprot_B or Gene_A equals Gene_B
    filtered_interactions = interactions_df[
        (interactions_df['Uniprot_A'] != interactions_df['Uniprot_B']) & 
        (interactions_df['Gene_A'] != interactions_df['Gene_B'])
    ]
    rows_after_id_filter = len(filtered_interactions)
    
    # Filter interactions to keep only those where both genes are in our gene set
    filtered_interactions = filtered_interactions[
        filtered_interactions['Gene_A'].isin(gene_set) & 
        filtered_interactions['Gene_B'].isin(gene_set)
    ]
    
    # Count rows with clinical significance for each gene using the filtered dataset
    print("Counting rows with clinical significance for each gene...")
    
    # Create new columns for individual counts
    filtered_interactions['Gene_A_count'] = filtered_interactions['Gene_A'].apply(
        lambda x: count_clinical_significance_rows(x, filtered_df)
    )
    filtered_interactions['Gene_B_count'] = filtered_interactions['Gene_B'].apply(
        lambda x: count_clinical_significance_rows(x, filtered_df)
    )
    
    # Calculate appearance counts for each gene
    gene_a_appearances = filtered_interactions['Gene_A'].value_counts().to_dict()
    gene_b_appearances = filtered_interactions['Gene_B'].value_counts().to_dict()

    # Add new columns for gene appearance counts
    filtered_interactions['Gene_A_appearance'] = filtered_interactions['Gene_A'].map(gene_a_appearances)
    filtered_interactions['Gene_B_appearance'] = filtered_interactions['Gene_B'].map(gene_b_appearances)
    
    # Sum the counts to create the appearance column
    filtered_interactions['appearance'] = filtered_interactions['Gene_A_count'] + filtered_interactions['Gene_B_count']
    
    return filtered_interactions, total_rows, rows_after_id_filter

def main():
    # File paths
    variant_file = '/home/xyzeng/Data/Uniprot/clinvarMain_ucsc.110424.tsv'
    interaction_file = '/home/xyzeng/Data/Uniprot/selected_interactions_20_unique.txt'
    
    try:
        # Process the files
        result, total_rows, rows_after_id_filter = process_variants_and_interactions(variant_file, interaction_file)
        
        # Save the results
        result.to_csv('filtered_interactions_with_counts.tsv', sep='\t', index=False)
        
        # Print summary statistics
        print("\nSummary:")
        print(f"Number of filtered interactions: {len(result)}")
        print(f"Number of unique genes in filtered interactions: {len(set(result['Gene_A']).union(set(result['Gene_B'])))}")
        print(f"Number of unique Uniprot IDs: {len(set(result['Uniprot_A']).union(set(result['Uniprot_B'])))}")
        print(f"\nRows removed due to same IDs: {total_rows - rows_after_id_filter}")
        print(f"Original number of rows: {total_rows}")
        print(f"Rows after removing same IDs: {rows_after_id_filter}")
        print(f"Final number of rows: {len(result)}")
        
        # Print statistics about the appearance counts
        print("\nAppearance Statistics:")
        print(f"Mean appearance count: {result['appearance'].mean():.2f}")
        print(f"Median appearance count: {result['appearance'].median():.2f}")
        print(f"Max appearance count: {result['appearance'].max()}")
        print(f"Min appearance count: {result['appearance'].min()}")
        
        # Print statistics about gene appearances
        print("\nGene Appearance Statistics:")
        print(f"Max Gene_A appearances: {result['Gene_A_appearance'].max()}")
        print(f"Max Gene_B appearances: {result['Gene_B_appearance'].max()}")
        print(f"Mean Gene_A appearances: {result['Gene_A_appearance'].mean():.2f}")
        print(f"Mean Gene_B appearances: {result['Gene_B_appearance'].mean():.2f}")
        
    except KeyError as e:
        print(f"Error: Column {e} not found in the data.")
        print("Please check the column names in your input files and update the script accordingly.")

if __name__ == "__main__":
    main()