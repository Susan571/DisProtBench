from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
import pandas as pd
import time

def extract_uniprot_ids_from_page(driver):
    """
    Extract UniProt IDs from the current page
    """
    # Wait for the table to load
    WebDriverWait(driver, 20).until(
        EC.presence_of_element_located((By.CLASS_NAME, "cdk-table"))
    )
    time.sleep(2)  # Additional small wait to ensure content is loaded
    
    # Find all cells in the UniProt ID column
    uniprot_cells = driver.find_elements(By.CSS_SELECTOR, "td.cdk-cell.cdk-column-acc span.btn-block")
    
    page_ids = []
    for cell in uniprot_cells:
        uniprot_id = cell.text.strip()
        if uniprot_id:
            page_ids.append(uniprot_id)
            
    return page_ids

def extract_all_uniprot_ids(base_url, total_pages=63):
    """
    Extract UniProt IDs from all pages
    """
    chrome_options = Options()
    chrome_options.add_argument('--headless')
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument('--disable-dev-shm-usage')
    
    all_uniprot_ids = []
    
    try:
        print("Starting Chrome...")
        driver = webdriver.Chrome(options=chrome_options)
        
        for page in range(total_pages):
            # Construct URL for current page
            url = f"{base_url}&page={page}"
            print(f"\nProcessing page {page + 1} of {total_pages}")
            
            # Access the page
            driver.get(url)
            
            # Extract IDs from current page
            page_ids = extract_uniprot_ids_from_page(driver)
            all_uniprot_ids.extend(page_ids)
            
            print(f"Found {len(page_ids)} UniProt IDs on page {page + 1}")
            
            # Optional: Add a short delay between pages to be nice to the server
            time.sleep(1)
            
        return all_uniprot_ids
        
    except Exception as e:
        print(f"Error during extraction: {e}")
        return all_uniprot_ids  # Return any IDs collected before the error
        
    finally:
        print("\nClosing browser...")
        driver.quit()

def save_to_csv(uniprot_ids, filename='uniprot_ids.csv'):
    """
    Save UniProt IDs to CSV file with additional information
    """
    df = pd.DataFrame({
        'UniProt_ID': uniprot_ids,
        'Source': ['DisProt'] * len(uniprot_ids)
    })
    df.to_csv(filename, index=False)
    print(f"\nData saved to {filename}")
    return df

if __name__ == "__main__":
    # Base URL without page parameter
    base_url = "https://www.disprot.org/browse?sort_field=disprot_id&sort_value=asc&page_size=20&release=current&show_ambiguous=true&show_obsolete=false&ncbi_taxon_id=9606"
    
    print("Starting UniProt ID extraction from all pages...")
    
    # Extract IDs from all pages
    all_uniprot_ids = extract_all_uniprot_ids(base_url)
    
    # Print results
    if all_uniprot_ids:
        print(f"\nTotal UniProt IDs found: {len(all_uniprot_ids)}")
        for idx, uniprot_id in enumerate(all_uniprot_ids[:5], 1):
            print(f"{idx}. {uniprot_id}")
        
        # Save to CSV with additional information
        df = save_to_csv(all_uniprot_ids)
        
        # Print summary statistics
        print("\nSummary:")
        print(f"Total unique UniProt IDs: {len(df['UniProt_ID'].unique())}")
        print(f"Data retrieval date: {df['Date_Retrieved'].iloc[0]}")
    else:
        print("\nNo UniProt IDs found.")
        
    print("\nExtraction complete!")