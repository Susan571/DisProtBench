import csv

def filter_sequences(input_file, output_file, min_length=10):
    """
    Filter sequences from a FASTA file based on minimum length requirement.
    
    Args:
        input_file (str): Path to input FASTA file
        output_file (str): Path to output FASTA file
        min_length (int): Minimum sequence length (default: 10)
    Returns:
        int: Number of sequences that met the length requirement
    """
    count = 0
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        header = None
        sequence = None
        
        for line in f_in:
            line = line.strip()
            
            if line.startswith('>'):
                if header and sequence and len(sequence) >= min_length:
                    f_out.write(f"{header}\n{sequence}\n")
                    count += 1
                    
                header = line
                sequence = None
            else:
                sequence = line
                
        if header and sequence and len(sequence) >= min_length:
            f_out.write(f"{header}\n{sequence}\n")
            count += 1
            
    return count

def extract_ids_from_fasta(input_file, output_csv):
    """
    Extract DisProt and UniProt IDs from FASTA file and save to CSV without headers.
    
    Args:
        input_file (str): Path to input FASTA file
        output_csv (str): Path to output CSV file
    Returns:
        list: List of tuples containing the extracted IDs
    """
    extracted_ids = []
    
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip()[1:].split('|')
                if len(parts) >= 2:
                    disprot_id = parts[0]
                    uniprot_id = parts[1]
                    extracted_ids.append([disprot_id, uniprot_id])
    
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(extracted_ids)
    
    return extracted_ids

def read_id_pairs(csv_file):
    """
    Read ID pairs from CSV file.
    
    Args:
        csv_file (str): Path to CSV file containing DisProt and UniProt IDs
    Returns:
        set: Set of tuples containing (DisProt_ID, UniProt_ID)
    """
    id_pairs = set()
    with open(csv_file, 'r') as f:
        for line in f:
            disprot_id, uniprot_id = line.strip().split(',')
            id_pairs.add((disprot_id, uniprot_id))
    return id_pairs

def match_sequences(fasta_file, id_pairs, output_file):
    """
    Extract sequences from FASTA file that match the ID pairs.
    
    Args:
        fasta_file (str): Path to input FASTA file
        id_pairs (set): Set of (DisProt_ID, UniProt_ID) tuples to match
        output_file (str): Path to output FASTA file
    Returns:
        int: Number of matched sequences
    """
    count = 0
    with open(fasta_file, 'r') as f_in, open(output_file, 'w') as f_out:
        write_sequence = False
        
        for line in f_in:
            if line.startswith('>'):
                parts = line.strip()[1:].split('|')
                write_sequence = False
                
                if len(parts) >= 2:
                    disprot_id = parts[0]
                    uniprot_id = parts[1]
                    
                    if (disprot_id, uniprot_id) in id_pairs:
                        write_sequence = True
                        f_out.write(line)
                        count += 1
            elif write_sequence:
                f_out.write(line)
                
    return count

def process_sequences(disorder_regions_file, full_sequences_file, min_length=20):
    """
    Process sequences through all steps: filtering, ID extraction, and matching.
    
    Args:
        disorder_regions_file (str): Path to disorder regions FASTA file
        full_sequences_file (str): Path to full sequences FASTA file
        min_length (int): Minimum sequence length for filtering
    """
    # Step 1: Filter sequences by length
    filtered_file = disorder_regions_file.replace('.fasta', f'_{min_length}.fasta')
    filtered_count = filter_sequences(disorder_regions_file, filtered_file, min_length)
    print(f"Filtered {filtered_count} sequences with length >= {min_length}")
    
    # Step 2: Extract IDs from filtered sequences
    ids_file = filtered_file.replace('.fasta', '_ID.csv')
    extracted_ids = extract_ids_from_fasta(filtered_file, ids_file)
    print(f"Extracted {len(extracted_ids)} ID pairs to {ids_file}")
    
    # Step 3: Match and extract full sequences
    id_pairs = read_id_pairs(ids_file)
    output_full_sequences = full_sequences_file.replace('.fasta', f'_{min_length}.fasta')
    matched_count = match_sequences(full_sequences_file, id_pairs, output_full_sequences)
    print(f"Matched and extracted {matched_count} sequences to {output_full_sequences}")

# Example usage
if __name__ == "__main__":
    # File paths
    disorder_regions_file = "/home/xyzeng/Data/Uniprot/disorder_regions.fasta"
    full_sequences_file = "/home/xyzeng/Data/Uniprot/full_sequences.fasta"
    min_length = 20
    
    # Process all steps
    process_sequences(disorder_regions_file, full_sequences_file, min_length)
    print("\nProcessing complete!")