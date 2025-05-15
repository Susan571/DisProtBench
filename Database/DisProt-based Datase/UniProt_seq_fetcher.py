import pandas as pd
import requests
from time import sleep
from collections import defaultdict
from datetime import datetime

def get_unique_uniprot_ids(df):
    """Extract unique UniProt IDs from both columns."""
    uniprot_ids = set(df['Uniprot_A'].unique()) | set(df['Uniprot_B'].unique())
    # Remove any isoform information (e.g., P38398-6 -> P38398)
    return {id.split('-')[0] for id in uniprot_ids}

def fetch_sequence(uniprot_id):
    """Fetch protein sequence from UniProt API."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            # Get the full FASTA entry including header
            fasta_lines = response.text.strip().split('\n')
            header = fasta_lines[0]
            sequence = ''.join(fasta_lines[1:])
            return header, sequence
        else:
            print(f"Failed to fetch sequence for {uniprot_id}: Status code {response.status_code}")
            return None, None
    except Exception as e:
        print(f"Error fetching sequence for {uniprot_id}: {str(e)}")
        return None, None

def format_sequence(sequence, width=60):
    """Format sequence with line breaks every 60 characters."""
    return '\n'.join(sequence[i:i+width] for i in range(0, len(sequence), width))

def main():
    # Read the CSV file
    df = pd.read_csv('/home/xyzeng/Data/Uniprot/gene_appearances_missense_selected.csv')
    
    # Get unique UniProt IDs
    unique_ids = get_unique_uniprot_ids(df)
    
    # Generate output filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d")
    output_file = f"uniprot_sequences_{timestamp}.fasta"
    
    # Counter for successful retrievals
    success_count = 0
    
    print(f"Fetching sequences for {len(unique_ids)} unique proteins...")
    
    # Open file in write mode
    with open(output_file, 'w') as f:
        for uniprot_id in unique_ids:
            print(f"Fetching sequence for {uniprot_id}")
            header, sequence = fetch_sequence(uniprot_id)
            
            if header and sequence:
                # Write the header and formatted sequence
                f.write(f"{header}\n")
                f.write(f"{format_sequence(sequence)}\n")
                success_count += 1
            
            # Add a small delay to be nice to the UniProt API
            sleep(1)
    
    print(f"Successfully retrieved {success_count} sequences")
    print(f"Results saved to: {output_file}")

if __name__ == "__main__":
    main()