import pandas as pd
import tkinter as tk
from tkinter import filedialog
import os

# --- CONFIGURATION ---
TOP_K = 10              # Number of stable anchors per region
MIN_EXPRESSION = 0.5    # Minimum mu_Reg to be considered expressed
MAX_CV = 0.5            # Maximum CV_Reg to be considered stable
# ---------------------

def identify_regional_anchors(file_path, top_k, min_mu, max_cv):
    """
    Processes a single file to identify top stable anchors based on kinetic metrics.
    """
    # 1. Load data and derive region name
    region_name = os.path.basename(file_path).split('_')[0]
    df = pd.read_csv(file_path)
    total_genes_in_file = len(df)

    # 2. Data Cleaning: Numeric conversion and NaN removal
    cols_to_clean = ['mu_Reg', 'CV_Reg']
    for col in cols_to_clean:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    df = df.dropna(subset=cols_to_clean)

    # 3. Apply Filters
    # Step 1: Stability filter (CV_Reg < 0.5)
    stable_df = df[df['CV_Reg'] < max_cv].copy()
    count_stable = len(stable_df)

    # Step 2: Expression eligibility (mu_Reg >= 0.5)
    stable_expressed_df = stable_df[stable_df['mu_Reg'] >= min_mu].copy()
    count_expressed = len(stable_expressed_df)

    # 4. Ranking: Sort by regional mean descending
    anchors_df = stable_expressed_df.sort_values(by='mu_Reg', ascending=False).head(top_k)
    
    # Add Region column for the combined output
    anchors_df.insert(0, 'Region', region_name)

    # 5. Prepare Summary Counts
    counts = {
        'Region': region_name,
        'Total_Genes': total_genes_in_file,
        'Stable_Genes': count_stable,
        'Stable_Expressed': count_expressed,
        'Anchors_Returned': len(anchors_df)
    }

    return anchors_df, counts

def main():
    # Setup Tkinter for file selection
    root = tk.Tk()
    root.withdraw()

    print(f"Opening file selector... (Targeting Top {TOP_K} Stable Anchors)")
    file_paths = filedialog.askopenfilenames(
        title="Select Regional Kinetic Analysis Files",
        filetypes=[("CSV files", "*.csv")]
    )

    if not file_paths:
        print("No files selected. Exiting.")
        return

    all_region_anchors = []
    summary_list = []

    # Process each region
    for path in file_paths:
        region_df, region_counts = identify_regional_anchors(path, TOP_K, MIN_EXPRESSION, MAX_CV)
        
        # Save individual region CSV
        region_name = region_counts['Region']
        individual_filename = f"Stable_Anchors_Top{TOP_K}_{region_name}.csv"
        region_df.to_csv(individual_filename, index=False)
        
        all_region_anchors.append(region_df)
        summary_list.append(region_counts)

    # Create Combined Dataframe
    combined_df = pd.concat(all_region_anchors, ignore_index=True)
    combined_filename = f"Stable_Anchors_Top{TOP_K}_AllRegions.csv"
    combined_df.to_csv(combined_filename, index=False)

    # Compute Global Anchors (Intersection)
    # Find genes that appear in the top K of EVERY selected region
    num_regions = len(file_paths)
    gene_counts = combined_df['Gene_Symbol'].value_counts()
    global_anchor_genes = gene_counts[gene_counts == num_regions].index.tolist()
    
    intersection_df = combined_df[combined_df['Gene_Symbol'].isin(global_anchor_genes)].copy()
    intersection_filename = "Stable_Anchors_Intersection_AllRegions.csv"
    intersection_df.to_csv(intersection_filename, index=False)

    # Output Results
    print("\n" + "="*60)
    print("REGIONAL STABILITY SUMMARY")
    print("="*60)
    summary_df = pd.DataFrame(summary_list)
    print(summary_df.to_string(index=False))
    print("="*60)

    print(f"\nSaved combined file: {combined_filename}")
    print(f"Saved intersection file: {intersection_filename}")
    print(f"Found {len(global_anchor_genes)} genes appearing in all regions.")
    
    print("\n--- COMBINED ANCHORS PREVIEW (HEAD) ---")
    print(combined_df.head(15))

if __name__ == "__main__":
    main()