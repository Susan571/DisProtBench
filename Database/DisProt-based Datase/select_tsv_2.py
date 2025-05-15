import pandas as pd

def process_variants_and_interactions(variant_file, interaction_file):
    """
    Process variant and interaction files based on specified criteria
    """
    # Read the variant TSV file with low_memory=False to avoid warnings
    print("Reading variant file...")
    variants_df = pd.read_csv(variant_file, sep='\t', low_memory=False)
    print(f"Initial variant rows: {len(variants_df)}")
    
    # Step 1: Filter by clinical significance (P or LP)
    filtered_df = variants_df[variants_df['_clinSignCode'].isin(['PG'])]
    print(f"Rows after clinical significance filter: {len(filtered_df)}")
    print("Unique clinical significance codes:", variants_df['_clinSignCode'].unique())
    
    # Step 2: Filter for missense variants
    filtered_df = filtered_df[filtered_df['molConseq'] == 'missense variant']
    print(f"Rows after missense variant filter: {len(filtered_df)}")
    print("Unique molecular consequences:", variants_df['molConseq'].unique())
    
    # Step 3: Keep only unique variants based on vcfDesc
    filtered_df = filtered_df.drop_duplicates(subset=['vcfDesc'])
    print(f"Rows after removing duplicate variants: {len(filtered_df)}")
    
    # Create a dictionary to store variant counts for each gene
    gene_variant_counts = {}
    
    # Count variants for each gene
    print("\nCounting variants for each gene...")
    for _, row in filtered_df.iterrows():
        if pd.notna(row['geneId']):
            genes = row['geneId'].split('|')
            for gene in genes:
                if gene in gene_variant_counts:
                    gene_variant_counts[gene] += 1
                else:
                    gene_variant_counts[gene] = 1
    
    print(f"Number of genes with variants: {len(gene_variant_counts)}")
    if gene_variant_counts:
        print("Sample of gene counts:", dict(list(gene_variant_counts.items())[:5]))
    
    # Read interaction file
    print("\nReading interaction file...")
    interactions_df = pd.read_csv(interaction_file, sep='\t')
    total_rows = len(interactions_df)
    print(f"Initial interaction rows: {total_rows}")
    
    # Print sample of interaction data
    print("\nSample of interaction data columns:", interactions_df.columns.tolist())
    print("Sample of first few Gene_A values:", interactions_df['Gene_A'].head().tolist())
    print("Sample of first few Gene_B values:", interactions_df['Gene_B'].head().tolist())
    
    # Remove rows where Uniprot_A equals Uniprot_B or Gene_A equals Gene_B
    filtered_interactions = interactions_df[
        (interactions_df['Uniprot_A'] != interactions_df['Uniprot_B']) & 
        (interactions_df['Gene_A'] != interactions_df['Gene_B'])
    ]
    rows_after_id_filter = len(filtered_interactions)
    print(f"\nRows after removing self-interactions: {rows_after_id_filter}")
    
    # Check overlap between interaction genes and variant genes
    interaction_genes_a = set(filtered_interactions['Gene_A'])
    interaction_genes_b = set(filtered_interactions['Gene_B'])
    variant_genes = set(gene_variant_counts.keys())
    
    print(f"\nNumber of unique genes in interactions: {len(interaction_genes_a.union(interaction_genes_b))}")
    print(f"Number of genes with variants: {len(variant_genes)}")
    print(f"Number of overlapping genes: {len(variant_genes.intersection(interaction_genes_a.union(interaction_genes_b)))}")
    
    # Filter interactions to keep only those where both genes have variants
    filtered_interactions = filtered_interactions[
        filtered_interactions['Gene_A'].isin(gene_variant_counts.keys()) & 
        filtered_interactions['Gene_B'].isin(gene_variant_counts.keys())
    ]
    print(f"\nFinal number of interactions after gene filter: {len(filtered_interactions)}")
    
    if len(filtered_interactions) > 0:
        # Add variant counts for each gene
        filtered_interactions['Gene_A_variants'] = filtered_interactions['Gene_A'].map(gene_variant_counts)
        filtered_interactions['Gene_B_variants'] = filtered_interactions['Gene_B'].map(gene_variant_counts)
        
        # Calculate total variants count
        filtered_interactions['total_variants'] = (
            filtered_interactions['Gene_A_variants'] + 
            filtered_interactions['Gene_B_variants']
        )
        
        # Sort by total variants count in descending order
        filtered_interactions = filtered_interactions.sort_values(
            by='total_variants', 
            ascending=False
        )
    
    return filtered_interactions, total_rows, rows_after_id_filter

def main():
    # File paths
    variant_file = '/home/xyzeng/Data/Uniprot/clinvarMain_ucsc.110424.tsv'
    interaction_file = '/home/xyzeng/Data/Uniprot/selected_interactions_20_unique.txt'
    
    try:
        # Process the files
        result, total_rows, rows_after_id_filter = process_variants_and_interactions(variant_file, interaction_file)
        
        if len(result) > 0:
            # Save the results
            result.to_csv('filtered_interactions_with_variant_counts.tsv', sep='\t', index=False)
            
            # Print summary statistics
            print("\nFinal Summary:")
            print(f"Number of filtered interactions: {len(result)}")
            print(f"Number of unique genes in filtered interactions: {len(set(result['Gene_A']).union(set(result['Gene_B'])))}")
            print(f"Number of unique Uniprot IDs: {len(set(result['Uniprot_A']).union(set(result['Uniprot_B'])))}")
            
            # Print statistics about the variant counts
            print("\nVariant Count Statistics:")
            print(f"Mean total variants: {result['total_variants'].mean():.2f}")
            print(f"Median total variants: {result['total_variants'].median():.2f}")
            print(f"Max total variants: {result['total_variants'].max()}")
            print(f"Min total variants: {result['total_variants'].min()}")
            
            # Print top 10 interactions by variant count
            print("\nTop 10 Interactions by Total Variant Count:")
            top_10 = result.head(10)[['Gene_A', 'Gene_B', 'Gene_A_variants', 'Gene_B_variants', 'total_variants']]
            print(top_10.to_string(index=False))
        else:
            print("\nNo matching interactions found after filtering!")
            print("Please check the filtering criteria and the overlap between variant and interaction data.")
        
    except KeyError as e:
        print(f"Error: Column {e} not found in the data.")
        print("Please check the column names in your input files and update the script accordingly.")
    except Exception as e:
        print(f"Error occurred: {str(e)}")

if __name__ == "__main__":
    main()