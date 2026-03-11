import pandas as pd
import tkinter as tk
from tkinter import filedialog
import os
import sys

def select_files(title):
    """Opens a pop-up window to select multiple CSV files."""
    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    file_paths = filedialog.askopenfilenames(
        title=title,
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
    )
    root.destroy()
    return list(file_paths)

# --- 1. Interactive File Selection ---
print("Please select the 5 Regional Kinetic Analysis files...")
files = select_files("Select All 5 Analysis CSVs")

if not files:
    print("Error: No files selected. Exiting.")
    sys.exit()

results = []

# --- 2. Process Each File ---
for file_path in files:
    # Get region name from filename
    region_name = os.path.basename(file_path).replace('_Kinetic_Analysis.csv', '').replace('_', ' ')
    
    df = pd.read_csv(file_path)
    
    # Clean column names just in case
    df.columns = df.columns.str.strip()
    
    # Check for required columns
    required = ['Gene_Symbol', 'mu_Reg', 'Enrichment_Score_E', 'Kinetic_Classification']
    if not all(col in df.columns for col in required):
        print(f"Skipping {region_name}: Required columns missing.")
        continue

    # Identify Highest Mean among Stable Genes
    stable_df = df[df['Kinetic_Classification'] == 'Stable']
    top_stable = None
    if not stable_df.empty:
        top_stable = stable_df.loc[stable_df['mu_Reg'].idxmax()]
    
    # Identify Highest Enrichment among Dynamic Genes
    dynamic_df = df[df['Kinetic_Classification'] == 'Dynamic']
    top_dynamic = None
    if not dynamic_df.empty:
        top_dynamic = dynamic_df.loc[dynamic_df['Enrichment_Score_E'].idxmax()]
        
    results.append({
        'Region': region_name,
        'Top_Stable_Gene': top_stable['Gene_Symbol'] if top_stable is not None else "N/A",
        'Stable_mu_Reg': top_stable['mu_Reg'] if top_stable is not None else 0,
        'Top_Dynamic_Gene': top_dynamic['Gene_Symbol'] if top_dynamic is not None else "N/A",
        'Dynamic_Enrichment_E': top_dynamic['Enrichment_Score_E'] if top_dynamic is not None else 0
    })

# --- 3. Display Results ---
summary_df = pd.DataFrame(results)

print("\n" + "="*80)
print("REGIONAL VULNERABILITY SUMMARY")
print("="*80)
print(summary_df.to_string(index=False))
print("="*80)

# Save the summary to a file
summary_df.to_csv("Region_Top_Metrics_Summary.csv", index=False)
print("\nSummary saved to: Region_Top_Metrics_Summary.csv")