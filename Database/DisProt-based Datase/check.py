import csv

def print_vcf_desc(filename, num_entries=10):
    with open(filename, 'r') as file:
        # Read TSV with tab delimiter
        reader = csv.DictReader(file, delimiter='\t')
        
        # Print the specified number of vcfDesc entries
        for i, row in enumerate(reader):
            if i >= num_entries:
                break
            print(f"Entry {i + 1}: {row['vcfDesc']}")

# Assuming the file is saved as 'variants.tsv'
print_vcf_desc('/home/xyzeng/Data/Uniprot/clinvarMain_ucsc.110424.tsv')