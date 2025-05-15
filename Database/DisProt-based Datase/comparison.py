import pandas as pd
import csv

def load_csv_data(csv_file):
    return pd.read_csv(csv_file)

def load_text_file_data(text_file):
    uniprot_ids = set()
    text_data = []
    with open(text_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        header = next(reader)
        for row in reader:
            uniprot_ids.add(row[0])  # Uniprot_A
            uniprot_ids.add(row[1])  # Uniprot_B
            text_data.append(row)
    return uniprot_ids, text_data, header

def select_matching_sequences(csv_data, text_file_uniprot_ids):
    return csv_data[csv_data['UniProt_ID'].isin(text_file_uniprot_ids)]

def select_matching_text_rows(text_data, selected_uniprot_ids):
    return [row for row in text_data if row[0] in selected_uniprot_ids or row[1] in selected_uniprot_ids]

def save_results(selected_csv_data, selected_text_data, text_header, csv_output_file, text_output_file):
    selected_csv_data.to_csv(csv_output_file, index=False)
    
    with open(text_output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(text_header)
        writer.writerows(selected_text_data)

def main():
    csv_file = "/home/xyzeng/Data/Uniprot/protein_ids.csv"  # Update this to your CSV file name
    text_file = "/home/xyzeng/Data/Uniprot/HomoSapiens_binary_hq.txt"  # Update this to your text file name
    csv_output_file = "selected_sequences.csv"
    text_output_file = "selected_interactions.txt"

    print("Loading data...")
    csv_data = load_csv_data(csv_file)
    text_file_uniprot_ids, text_data, text_header = load_text_file_data(text_file)

    print("Selecting matching sequences...")
    selected_csv_data = select_matching_sequences(csv_data, text_file_uniprot_ids)
    
    print("Selecting matching text rows...")
    selected_uniprot_ids = set(selected_csv_data['UniProt_ID'])
    selected_text_data = select_matching_text_rows(text_data, selected_uniprot_ids)

    print("Saving results...")
    save_results(selected_csv_data, selected_text_data, text_header, csv_output_file, text_output_file)

    print(f"Process complete.")
    print(f"Selected {len(selected_csv_data)} sequences, saved to {csv_output_file}")
    print(f"Selected {len(selected_text_data)} interaction rows, saved to {text_output_file}")

if __name__ == "__main__":
    main()