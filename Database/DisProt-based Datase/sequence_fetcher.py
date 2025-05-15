import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from bs4 import BeautifulSoup
import re
import time

# Set up ChromeDriver with Selenium
chrome_driver_path = "/opt/homebrew/bin/chromedriver"  # Update this path to your ChromeDriver location
output_disorder_fasta = "disorder_regions.fasta"
output_full_fasta = "full_sequences.fasta"  # New file for full sequences

# Load the CSV file
csv_file = "protein_ids.csv"  # Update this to your CSV file name
data = pd.read_csv(csv_file)

# Initialize the FASTA outputs
with open(output_disorder_fasta, "w") as disorder_file, open(output_full_fasta, "w") as full_file:

    # Set up Selenium options
    service = Service(chrome_driver_path)
    chrome_options = Options()
    chrome_options.add_argument("--headless")  # Run in headless mode
    browser = webdriver.Chrome(service=service, options=chrome_options)

    # Iterate through each row in the CSV file
    for index, row in data.iterrows():
        disprot_id = row['DisProt_ID']  # Assuming column name in CSV is 'DisProt_ID'
        uniprot_id = row['UniProt_ID']  # Assuming column name in CSV is 'UniProt_ID'

        # Fetch DisProt page
        url = f"https://www.disprot.org/{disprot_id}"
        browser.get(url)
        time.sleep(2)  # Wait for the page to load fully

        # Parse page source
        html_content = browser.page_source
        soup = BeautifulSoup(html_content, 'html.parser')

        # Extract the full protein sequence
        sequence_match = re.search(r'"hasBioPolymerSequence":"([A-Z]+)"', html_content)
        if sequence_match:
            full_sequence = sequence_match.group(1)
            # Write full sequence to FASTA file
            full_file.write(f">{disprot_id}|{uniprot_id}|full\n")
            full_file.write(full_sequence + "\n")
        else:
            print(f"Full sequence not found for {disprot_id}")
            continue

        # Extract disorder regions
        disorder_sequence = ""
        disorder_tags = soup.find_all(id=re.compile(r"^f_disprot_consensus_(\d+)-(\d+)$"))

        for tag in disorder_tags:
            match = re.search(r"f_disprot_consensus_(\d+)-(\d+)", tag['id'])
            if match:
                start_position = int(match.group(1)) - 1  # Convert to zero-based index
                end_position = int(match.group(2))
                disorder_sequence += full_sequence[start_position:end_position]

        # Write disorder regions to FASTA if found
        if disorder_sequence:
            disorder_file.write(f">{disprot_id}|{uniprot_id}|disorder\n")
            disorder_file.write(disorder_sequence + "\n")
        else:
            print(f"No disorder regions found for {disprot_id}")

    # Close browser
    browser.quit()

print(f"Disorder regions saved in {output_disorder_fasta}")
print(f"Full sequences saved in {output_full_fasta}")