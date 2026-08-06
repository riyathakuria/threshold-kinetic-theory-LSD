import pandas as pd
import requests
import time

# 1. Define File Paths (Updated to use the 'Raw data' folder)
ROOT = Path(__file__).resolve().parents[2]

def get_uniprot_sequence(uniprot_id):
    """Fetches the canonical FASTA sequence from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        return "".join(response.text.split('\n')[1:])
    return None

def predict_palmitoylation_sites(sequence):
    """
    Identifies candidate Cysteine (C) sites.
    (Note: Plug in your local MusiteDeep CLI command here if running the deep learning model locally).
    """
    candidates = []
    for i, amino_acid in enumerate(sequence):
        if amino_acid == 'C':
            candidates.append(f"C{i + 1}")
    return candidates

# --- MAIN PIPELINE EXECUTION ---

print("1. Loading deduplicated dataset from 'Raw data' folder...")
df = pd.read_csv(input_csv)
results = []

print("\n2. Fetching sequences and predicting sites...")
print("(This will process one gene per row and compile all sites into a single cell)")

# Iterate through the dataset
for index, row in df.iterrows():
    # Safely extract the UniProt ID
    uniprot_str = str(row['UniProtIds'])
    if uniprot_str == 'nan': 
        continue
    uniprot_id = uniprot_str.split(';')[0]
    
    # Standardize gene name to UPPERCASE for later Allen Brain Atlas mapping
    gene = str(row['Genes']).upper()
    
    print(f"Processing {gene} ({uniprot_id})...")
    
    # Fetch sequence
    seq = get_uniprot_sequence(uniprot_id)
    
    if seq:
        # Predict sites
        sites = predict_palmitoylation_sites(seq)
        
        if sites:
            # Append as a single row per gene
            results.append({
                'Gene': gene,
                'UniProtID': uniprot_id,
                'Total_Predicted_Sites': len(sites),
                'Cysteine_Positions': ", ".join(sites) # Groups all sites into one string separated by commas
            })
            
    # Politeness delay to prevent UniProt from blocking your IP
    time.sleep(0.5) 

print("\n3. Saving final results...")
final_df = pd.DataFrame(results)

if not final_df.empty:
    final_df.to_csv(output_csv, index=False)
    print(f"Success! {len(final_df)} genes processed and saved to:\n{output_csv}")
else:
    print("No sites found.")
