import pandas as pd
import numpy as np
import re

# --- 1. CONFIGURATION ---
ROOT = Path(__file__).resolve().parents[3]

lineages = {
    'Striatum': ['medial ganglionic eminence', 'lateral ganglionic eminence', 'striatum'],
    'Cerebellum': ['upper (rostral) rhombic lip', 'cerebellum', 'cerebellar cortex'],
    'Hippocampus': ['hippocampus (hippocampal formation)'],
    'Visual Cortex': ['occipital neocortex', 'primary visual cortex (striate cortex, area V1/17)']
}

target_ages = ['8 pcw', '9 pcw', '12 pcw', '13 pcw', '16 pcw', '17 pcw', '19 pcw', 
               '21 pcw', '24 pcw', '25 pcw', '26 pcw', '35 pcw', '37 pcw',
               '4 mos', '10 mos', '1 yrs', '2 yrs', '3 yrs', '4 yrs', '8 yrs', '11 yrs']

def clean_symbol(s):
    """Removes all non-alphanumeric characters and converts to uppercase for strict comparison."""
    return re.sub(r'[^A-Z0-9]', '', str(s).upper())

def perform_smart_mapping():
    # --- STEP 1: LOAD PPI LIST ---
    print("Step 1: Reading PPI cluster genes...")
    ppi_df = pd.read_excel(ppi_path)
    gene_col = [col for col in ppi_df.columns if 'gene' in col.lower()][0]
    
    # Create a mapping of CleanedSymbol -> OriginalSymbol
    target_genes_raw = ppi_df[gene_col].astype(str).str.strip().unique()
    target_map = {clean_symbol(g): g for g in target_genes_raw}
    print(f"-> Target List: {len(target_genes_raw)} unique genes.")

    # --- STEP 2: IDENTIFY COLUMNS ---
    print("\nStep 2: Identifying relevant columns from Atlas...")
    meta = pd.read_excel(matrix_path, nrows=4, header=None)
    regions_row = meta.iloc[0].tolist() 
    ages_row = meta.iloc[3].tolist() 

    valid_cols = [0]
    col_metadata = {} 
    for i in range(1, len(regions_row)):
        reg_val = str(regions_row[i]).strip().lower()
        age_val = str(ages_row[i]).strip().lower()
        if any(age.lower() == age_val for age in target_ages):
            for lineage_name, loc_list in lineages.items():
                if any(loc.lower() in reg_val for loc in loc_list):
                    valid_cols.append(i)
                    col_metadata[i] = {'Lineage': lineage_name, 'Age': age_val}
                    break
    print(f"-> Identified {len(valid_cols)-1} valid columns.")

    # --- STEP 3: ROBUST GENE MATCHING ---
    print("\nStep 3: Finding gene row indices using Robust Matching...")
    df_atlas_symbols = pd.read_excel(matrix_path, usecols=[0], header=None)
    
    match_indices = []
    found_clean_symbols = set()

    for idx, val in df_atlas_symbols[0].items():
        # Clean the Atlas symbol the same way
        atlas_clean = clean_symbol(val)
        if atlas_clean in target_map:
            match_indices.append(idx)
            found_clean_symbols.add(atlas_clean)

    # Report missing for diagnostics
    missing = set(target_map.keys()) - found_clean_symbols
    if missing:
        print(f" Warning: {len(missing)} genes still not found.")
        print("Missing Cleaned Symbols:", [target_map[m] for m in missing])
    else:
        print(" Success: All genes found!")

    # --- STEP 4: EXTRACTION ---
    print(f"\nStep 4: Loading data for {len(match_indices)} rows...")
    rows_to_keep = set([0, 1, 2, 3, 4]) | set(match_indices)
    df_slice = pd.read_excel(matrix_path, usecols=valid_cols, 
                             skiprows=lambda x: x not in rows_to_keep, header=None)

    # --- STEP 5: TRANSFORM ---
    print("Step 5: Formatting results...")
    results = []
    for row_idx in range(5, len(df_slice)):
        # We use the official symbol from your PPI list for the final output
        atlas_val = str(df_slice.iloc[row_idx, 0]).strip()
        final_gene_name = target_map.get(clean_symbol(atlas_val), atlas_val.upper())
        
        for slice_col_idx, original_col_idx in enumerate(valid_cols):
            if original_col_idx == 0: continue
            meta_info = col_metadata[original_col_idx]
            val = pd.to_numeric(df_slice.iloc[row_idx, slice_col_idx], errors='coerce')
            if not np.isnan(val):
                results.append({
                    'Gene': final_gene_name,
                    'Age': meta_info['Age'],
                    'Lineage': meta_info['Lineage'],
                    'Expression': val
                })

    final_df = pd.DataFrame(results)
    final_df.to_csv(output_path, index=False)
    print(f"\n Done! File saved to: {output_path}")

if __name__ == "__main__":
    perform_smart_mapping()
