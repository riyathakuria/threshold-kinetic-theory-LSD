import pandas as pd
import tkinter as tk
from tkinter import filedialog
import os

# --- CONFIGURATION VARIABLES ---
ENRICH_CUTOFF = 1.0  # Change to 1.1 as needed
TOP_N = 3            # Number of drivers to return per region
# -------------------------------

def process_regional_drivers(file_path, enrich_cutoff, top_n):
    """
    Identifies dynamic drivers for a single region based on kinetic criteria.
    Returns a tuple of (drivers_dataframe, counts_dictionary).
    """
    # 1. Load data and identify region
    region_name = os.path.basename(file_path).split('_')[0]
    df = pd.read_csv(file_path)
    
    # 2. Data Cleaning: Numeric conversion and NaN removal
    cols_to_clean = ['mu_Reg', 'CV_Reg', 'Enrichment_Score_E']
    for col in cols_to_clean:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    total_genes = len(df)
    df = df.dropna(subset=cols_to_clean)
    
    # 3. Stepwise Filtering
    # Step 1: Expression eligibility (mu_Reg >= 0.5)
    pass_mu = df[df['mu_Reg'] >= 0.5].copy()
    
    # Step 2: Dynamic eligibility (CV_Reg >= 0.5)
    pass_cv = pass_mu[pass_mu['CV_Reg'] >= 0.5].copy()
    
    # Step 3: Region-enrichment eligibility (Enrichment > cutoff)
    pass_enrich = pass_cv[pass_cv['Enrichment_Score_E'] > enrich_cutoff].copy()
    
    # 4. Ranking and Selection
    # Step 4: Sort by Enrichment descending
    ranked_df = pass_enrich.sort_values(by='Enrichment_Score_E', ascending=False)
    
    # Step 5: Take top N
    drivers_df = ranked_df.head(top_n).copy()
    
    # Add Region column for the combined output
    drivers_df.insert(0, 'Region', region_name)
    
    # Store counts for the summary report
    counts = {
        'Region': region_name,
        'Total': total_genes,
        'Pass_mu': len(pass_mu),
        'Pass_CV': len(pass_cv),
        'Pass_Enrich': len(pass_enrich),
        'Final_Drivers': len(drivers_df)
    }
    
    return drivers_df, counts

def main():
    # Initialize Tkinter and hide the root window
    root = tk.Tk()
    root.withdraw()

    print(f"Select your Kinetic Analysis CSV files (Cutoff: {ENRICH_CUTOFF}, Top N: {TOP_N})")
    file_paths = filedialog.askopenfilenames(
        title="Select Regional Kinetic Analysis Files",
        filetypes=[("CSV files", "*.csv")]
    )

    if not file_paths:
        print("No files selected. Exiting.")
        return

    all_drivers = []
    summary_stats = []

    for path in file_paths:
        drivers_df, counts = process_regional_drivers(path, ENRICH_CUTOFF, TOP_N)
        
        # Save individual region CSV (even if empty)
        region = counts['Region']
        filename = f"Dynamic_Drivers_Top3_{region}_Enriched.csv"
        drivers_df.to_csv(filename, index=False)
        
        all_drivers.append(drivers_df)
        summary_stats.append(counts)

    # Create and print the summary table
    summary_df = pd.DataFrame(summary_stats)
    print("\n" + "="*70)
    print("STEPWISE FILTERING SUMMARY PER REGION")
    print("="*70)
    print(summary_df.to_string(index=False))
    print("="*70)

    # Combine all regions and save
    combined_drivers_df = pd.concat(all_drivers, ignore_index=True)
    combined_filename = "Dynamic_Drivers_Top3_AllRegions_Enriched.csv"
    combined_drivers_df.to_csv(combined_filename, index=False)

    print(f"\nCombined drivers saved to: {combined_filename}")
    print("\nFINAL COMBINED DRIVERS (HEAD):")
    print("-" * 30)
    if not combined_drivers_df.empty:
        print(combined_drivers_df)
    else:
        print("No genes met the combined criteria across selected regions.")

if __name__ == "__main__":
    main()