def format_sequences(input_file, output_file):
    """
    Reads a sequence file and reformats it so each sequence takes exactly two lines:
    1. Header line (starting with >sp)
    2. Sequence line (all on one line without breaks)
    
    Args:
        input_file (str): Path to input file
        output_file (str): Path to output file
    """
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        current_header = None
        current_sequence = []
        
        for line in f_in:
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
                
            # If we find a header line (starts with >sp)
            if line.startswith('>sp'):
                # If we have a previous sequence to write
                if current_header is not None:
                    # Write the previous sequence
                    f_out.write(current_header + '\n')
                    f_out.write(''.join(current_sequence) + '\n')
                    current_sequence = []
                
                # Store the new header
                current_header = line
            else:
                # Add to current sequence
                current_sequence.append(line)
        
        # Write the last sequence
        if current_header is not None:
            f_out.write(current_header + '\n')
            f_out.write(''.join(current_sequence) + '\n')

# Example usage
if __name__ == "__main__":
    input_file = "/home/xyzeng/Data/Uniprot/uniprot_sequences_20241102.fasta"
    output_file = "final_data.txt"
    format_sequences(input_file, output_file)