import pandas as pd

# Define your local file paths
musite_path = r'/Users/riyathakuria/Documents/LSD_Insilico_Analysis_24:3:26/Raw data/MusiteDeep_Palmitoylation_Predictions.csv'
swisspalm_path = r'/Users/riyathakuria/Documents/LSD_Insilico_Analysis_24:3:26/Raw data/SwissPalm Release 5 for Homo sapiens.txt'
output_path = r'/Users/riyathakuria/Documents/LSD_Insilico_Analysis_24:3:26/Raw data/Master_LSD_Palmitoylome_Integrated.csv'

def generate_master_file(musite_csv, swisspalm_txt, save_path):
    # 1. Load your MusiteDeep predictions
    print("Loading MusiteDeep predictions...")
    df = pd.read_csv(musite_csv)

    # 2. Parse SwissPalm database to extract verified IDs and Gene Symbols
    print("Parsing SwissPalm database...")
    verified_ids = set()
    verified_genes = set()
    
    with open(swisspalm_txt, 'r') as f:
        for line in f:
            # Skip comment lines and header
            if line.startswith('#') or line.startswith('UniProt AC') or not line.strip():
                continue
            
            parts = line.split('\t')
            if len(parts) >= 2:
                uniprot_ac = parts[0].strip().upper()
                gene_symbol = parts[1].strip().upper()
                verified_ids.add(uniprot_ac)
                verified_genes.add(gene_symbol)

    # 3. Define the logic for verification status
    def identify_status(row):
        # Clean and uppercase identifiers for robust matching
        u_id = str(row['UniProtID']).strip().upper()
        gene = str(row['Gene']).strip().upper()
        
        if u_id in verified_ids or gene in verified_genes:
            return 'SwissProt Verified'
        else:
            return 'Novel Prediction'

    # 4. Apply mapping and create the new column
    print("Mapping proteins to SwissPalm...")
    df['Validation_Status'] = df.apply(identify_status, axis=1)

    # 5. Save the integrated master file
    df.to_csv(save_path, index=False)
    
    # Summary report
    counts = df['Validation_Status'].value_counts()
    print("-" * 30)
    print(f"Master file saved to: {save_path}")
    print(f"Total Proteins: {len(df)}")
    print(f"Verified: {counts.get('SwissProt Verified', 0)}")
    print(f"Novel: {counts.get('Novel Prediction', 0)}")
    print("-" * 30)

# Run the function
if __name__ == "__main__":
    generate_master_file(musite_path, swisspalm_path, output_path)