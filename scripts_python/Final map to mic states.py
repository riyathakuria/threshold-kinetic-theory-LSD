import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import filedialog
import os

def select_file(title):
    root = tk.Tk()
    root.withdraw()
    print(f"Please select: {title}")
    path = filedialog.askopenfilename(title=title, filetypes=[("CSV files", "*.csv")])
    return path

def map_microglia_states():
    # 1. Interactive File Selection
    gene_analysis_path = select_file("Select your 'Final_Gene_Analysis.csv'")
    microglia_states_path = select_file("Select your 'Microglia_States_Named.csv'")

    if not gene_analysis_path or not microglia_states_path:
        print("File selection cancelled.")
        return

    # 2. Load Data (First column as Index/Gene Names)
    # Final_Gene_Analysis has Gene names in the first column (col 0)
    df_genes = pd.read_csv(gene_analysis_path, index_col=0)
    # Microglia_States_Named has Gene names in the first column
    df_states = pd.read_csv(microglia_states_path, index_col=0)

    # 3. Clean Indices (Removing quotes if present)
    df_states.index = df_states.index.str.replace('"', '').str.strip()
    df_genes.index = df_genes.index.str.replace('"', '').str.strip()

    # 4. Z-Score Normalization (to find the 'Dominant' State)
    # We normalize across the states (rows) to see where the gene is most active
    df_states_z = df_states.apply(lambda x: (x - x.mean()) / x.std() if x.std() != 0 else 0, axis=1)

    # 5. Mapping Logic
    results = []
    for gene in df_genes.index:
        if gene in df_states_z.index:
            gene_data = df_states_z.loc[gene]
            # Find the state with the highest Z-score
            primary_state = gene_data.idxmax()
            z_score_val = gene_data.max()
            
            # Combine raw values and the identified primary state
            row = df_states.loc[gene].to_dict()
            row['Primary_Microglia_State'] = primary_state
            row['State_Z_Score'] = round(z_score_val, 3)
        else:
            # If gene not found, mark as N/A
            row = {col: "N/A" for col in df_states.columns}
            row['Primary_Microglia_State'] = "N/A"
            row['State_Z_Score'] = "N/A"
        
        row['Gene'] = gene
        results.append(row)

    # 6. Final Dataframe Assembly
    output_df = pd.DataFrame(results)
    
    # Reorder columns to put Gene and Primary State first
    cols = ['Gene', 'Primary_Microglia_State', 'State_Z_Score'] + [c for c in df_states.columns]
    output_df = output_df[cols]

    # 7. Save to CSV
    output_name = "Lastly_Mapped_Gene_Microglia_States.csv"
    output_df.to_csv(output_name, index=False)
    
    print("\n" + "="*50)
    print(f"SUCCESS: Mapping Complete.")
    print(f"Saved as: {output_name}")
    print("="*50)
    print(output_df.head(10))

if __name__ == "__main__":
    map_microglia_states()