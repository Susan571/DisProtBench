import csv

def read_disprot_uniprot_pairs(csv_file):
    """
    Read DisProt-UniProt ID pairs from previously saved CSV file.
    
    Args:
        csv_file (str): Path to CSV file containing DisProt and UniProt IDs
    Returns:
        list: List of UniProt IDs maintaining all entries (including duplicates)
    """
    uniprot_ids = []
    pairs_count = 0
    with open(csv_file, 'r') as f:
        for line in f:
            pairs_count += 1
            disprot_id, uniprot_id = line.strip().split(',')
            uniprot_ids.append(uniprot_id)
    
    print(f"Total pairs in CSV: {pairs_count}")
    print(f"Unique UniProt IDs: {len(set(uniprot_ids))}")
    return uniprot_ids

def match_interactions(interaction_file, uniprot_ids, output_file):
    """
    Match and save interactions where either Uniprot_A or Uniprot_B exists in our dataset.
    
    Args:
        interaction_file (str): Path to interaction file
        uniprot_ids (list): List of UniProt IDs from our sequence dataset
        output_file (str): Path to output file for matched interactions
    Returns:
        tuple: (matched_count, total_count) of interactions
    """
    matched_count = 0
    total_count = 0
    uniprot_set = set(uniprot_ids)  # Convert to set for faster lookup
    
    with open(interaction_file, 'r') as f_in, open(output_file, 'w') as f_out:
        # Write header
        header = f_in.readline()
        f_out.write(header)
        
        # Process interactions
        for line in f_in:
            total_count += 1
            parts = line.strip().split('\t')
            
            if len(parts) >= 2:
                uniprot_a = parts[0]
                uniprot_b = parts[1]
                
                # If either Uniprot_A or Uniprot_B is in our dataset, save the interaction
                if uniprot_a in uniprot_set or uniprot_b in uniprot_set:
                    f_out.write(line)
                    matched_count += 1
    
    return matched_count, total_count

def process_files(csv_file, interaction_file, output_file):
    """
    Process files to match interactions with our sequence dataset.
    
    Args:
        csv_file (str): Path to CSV file with DisProt-UniProt pairs from previous step
        interaction_file (str): Path to interaction dataset
        output_file (str): Path to output file for matched interactions
    """
    # Read UniProt IDs from our sequence dataset
    print("Reading UniProt IDs from sequence dataset...")
    uniprot_ids = read_disprot_uniprot_pairs(csv_file)
    
    # Match interactions
    print("\nMatching interactions...")
    matched_count, total_count = match_interactions(interaction_file, uniprot_ids, output_file)
    print(f"Matched {matched_count} interactions out of {total_count} total interactions")
    
    # Calculate percentage
    if total_count > 0:
        percentage = (matched_count / total_count) * 100
        print(f"Match rate: {percentage:.2f}%")
    
    # Print matching statistics
    matches_a = 0
    matches_b = 0
    both_matches = 0
    uniprot_set = set(uniprot_ids)
    
    with open(output_file, 'r') as f:
        # Skip header
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                a_match = parts[0] in uniprot_set
                b_match = parts[1] in uniprot_set
                if a_match and b_match:
                    both_matches += 1
                elif a_match:
                    matches_a += 1
                elif b_match:
                    matches_b += 1
    
    print("\nMatching Statistics:")
    print(f"Matches in Uniprot_A only: {matches_a}")
    print(f"Matches in Uniprot_B only: {matches_b}")
    print(f"Matches in both columns: {both_matches}")

# Example usage
if __name__ == "__main__":
    # File paths
    csv_file = "/home/xyzeng/Data/Uniprot/disorder_regions_20_ID.csv"  # CSV file from previous step
    interaction_file = "/home/xyzeng/Data/Uniprot/selected_interactions.txt"  # Your new interaction dataset
    output_file = "/home/xyzeng/Data/Uniprot/selected_interactions_20.txt"   # Where to save matched interactions
    
    try:
        process_files(csv_file, interaction_file, output_file)
        print("\nProcessing complete!")
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
    except Exception as e:
        print(f"Error occurred: {e}")