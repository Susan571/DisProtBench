import pandas as pd
from collections import Counter

def count_gene_appearances(summary_file):
    """Count appearances of each gene in summary data"""
    # Load summary data genes
    summary_data = pd.read_csv(summary_file, sep='\t')
    gene_counts = Counter(summary_data['SubmittedGeneSymbol'])
    return dict(gene_counts)

def rank_interactions(interaction_file, gene_counts):
    """Rank interactions based on total appearances of their genes"""
    # Load interaction data
    interactions = pd.read_csv(interaction_file, sep='\t')
    
    # Calculate total appearances for each row
    def get_total_appearances(row):
        return gene_counts.get(row['Gene_A'], 0) + gene_counts.get(row['Gene_B'], 0)
    
    # Add total appearances column
    interactions['total_appearances'] = interactions.apply(get_total_appearances, axis=1)
    
    # Filter for total appearances > 2 and sort
    filtered_interactions = interactions[interactions['total_appearances'] > 2]
    ranked_interactions = filtered_interactions.sort_values('total_appearances', ascending=False)
    
    return ranked_interactions

def save_to_excel(ranked_interactions, gene_counts, output_file):
    """Save data to Excel file using available engine"""
    try:
        # Try to save with basic formatting
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # Save ranked interactions to first sheet
            ranked_interactions.to_excel(writer, sheet_name='Ranked_Interactions', index=False)
            
            # Save gene counts to second sheet
            counts_df = pd.DataFrame.from_dict(gene_counts, orient='index', columns=['appearances'])
            counts_df.sort_values('appearances', ascending=False, inplace=True)
            counts_df.to_excel(writer, sheet_name='Gene_Appearances')
            
        return True
    except Exception as e:
        print(f"Warning: Could not save as Excel, will save as CSV instead. Error: {e}")
        # Save as CSV if Excel saving fails
        ranked_interactions.to_csv('ranked_interactions.csv', sep='\t', index=False)
        
        counts_df = pd.DataFrame.from_dict(gene_counts, orient='index', columns=['appearances'])
        counts_df.sort_values('appearances', ascending=False, inplace=True)
        counts_df.to_csv('gene_appearances.csv', sep='\t')
        return False

def main():
    summary_file = 'filtered_summary_data_2.txt'
    interaction_file = 'selected_interactions_20_unique.txt'
    output_file = 'ranked_interactions.xlsx'
    
    print("Counting gene appearances in summary data...")
    gene_counts = count_gene_appearances(summary_file)
    
    print("Ranking interactions...")
    interactions = pd.read_csv(interaction_file, sep='\t')
    ranked_interactions = rank_interactions(interaction_file, gene_counts)
    
    print("\nRanking Statistics:")
    print(f"Original interaction rows: {len(interactions)}")
    print(f"Filtered interaction rows (appearances > 2): {len(ranked_interactions)}")
    print(f"Rows removed: {len(interactions) - len(ranked_interactions)}")
    
    # Show distribution of total appearances
    print("\nAppearance Distribution in filtered interactions:")
    print(ranked_interactions['total_appearances'].describe())
    
    # Show top 5 interactions by appearances
    print("\nTop 5 interactions by total appearances:")
    print(ranked_interactions[['Gene_A', 'Gene_B', 'total_appearances']].head())
    
    # Save the results
    success = save_to_excel(ranked_interactions, gene_counts, output_file)
    
    if success:
        print(f"\nResults saved to {output_file} with two sheets:")
        print("1. Ranked_Interactions - Contains the filtered and ranked interactions")
        print("2. Gene_Appearances - Contains the count of appearances for each gene")
    else:
        print("\nResults saved as separate CSV files:")
        print("1. ranked_interactions.csv")
        print("2. gene_appearances.csv")

if __name__ == "__main__":
    main()