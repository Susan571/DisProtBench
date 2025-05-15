import pandas as pd

def load_summary_data(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Find the start of the actual data
    data_start = 0
    for i, line in enumerate(lines):
        if line.startswith('#VariationID\t'):
            data_start = i
            break

    # Extract header and data
    header = lines[data_start].strip().split('\t')
    data = [line.strip().split('\t') for line in lines[data_start+1:] if line.strip() and not line.startswith('#')]

    # Create DataFrame
    df = pd.DataFrame(data, columns=header)
    return df

def process_data(summary_data):
    # Filter summary data with conditions:
    filtered_data = summary_data[
        # Remove empty/invalid gene symbols
        (summary_data['SubmittedGeneSymbol'] != '-') &
        (summary_data['SubmittedGeneSymbol'].notna()) &
        # Filter for pathogenic variants
        (summary_data['ClinicalSignificance'] == 'Pathogenic') &
        # Filter for variant or missense in description
        (
            summary_data['Description'].str.contains('variant', case=False, na=False) |
            summary_data['Description'].str.contains('missense', case=False, na=False)
        )
    ]
    
    return filtered_data

def main():
    summary_file = '/home/xyzeng/Data/Uniprot/submission_summary.txt'
    output_file = 'filtered_summary_data_2.txt'

    print("Loading summary data...")
    summary_data = load_summary_data(summary_file)

    print("Processing data...")
    filtered_data = process_data(summary_data)

    print("Saving results...")
    filtered_data.to_csv(output_file, sep='\t', index=False)

    print("\nFiltering Statistics:")
    print(f"Original rows: {len(summary_data)}")
    print(f"Pathogenic variants: {len(summary_data[summary_data['ClinicalSignificance'] == 'Pathogenic'])}")
    print(f"Variants with 'variant' or 'missense' in description: {len(summary_data[summary_data['Description'].str.contains('variant|missense', case=False, na=False)])}")
    print(f"Final selected variants: {len(filtered_data)}")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()