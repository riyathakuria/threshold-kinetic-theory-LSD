import subprocess
import pandas as pd
import os

# 1. Define your exact file paths
large_tsv_path = "/Users/riyathakuria/Documents/LSD_Insilico_Analysis_24:3:26/240526_JH_NP_28_LysoIP_Fcrls_Microglia_candidates.tsv"
lsd_csv_path = "/Users/riyathakuria/Documents/Zip files/LSD final insilico/Protein list/AA Unique protein.csv"

# Temporary pattern file and TSV output for the grep stream
pattern_file = "temp_lsd_patterns.txt"
temp_tsv_path = "/Users/riyathakuria/Documents/LSD_Insilico_Analysis_24:3:26/temp_mapped_output.tsv"

# NEW: Saving to an "Optimized" file to see the difference!
final_csv_path = "/Users/riyathakuria/Documents/LSD_Insilico_Analysis_24:3:26/Mapped_LSD_Microglia_Overlap_Optimized.csv"

# 2. Extract the LSD proteins to use as search patterns
print("1. Extracting LSD proteins for the search query...")
lsd_df = pd.read_csv(lsd_csv_path)

# Filter out the stray 'x' row so grep doesn't match random x's
protein_list = lsd_df[lsd_df['Gene'] != 'x']['Gene'].dropna().unique() 

with open(pattern_file, 'w') as f:
    for protein in protein_list:
        f.write(f"{protein}\n")

# 3. Execute grep to safely process the 20GB file
print("2. Initiating Case-Insensitive grep stream processing. Bypassing 8GB RAM limit...")
with open(temp_tsv_path, "w") as outfile:
    # Grab the header row from your 20GB file first
    subprocess.run(["head", "-n", "1", large_tsv_path], stdout=outfile)
    
    # THE FIX: Added '-i' for case-insensitive matching across Mouse/Human nomenclature
    # -i: ignore case (Matches HEXA with Hexa)
    # -F: fixed strings
    # -w: whole word (prevents APP from matching APPLE)
    subprocess.run(["grep", "-i", "-F", "-w", "-f", pattern_file, large_tsv_path], stdout=outfile)

# 4. Convert the shrunken, mapped data to CSV
print("3. Converting mapped data to a clean CSV format...")
try:
    # Load the small, filtered dataset into Pandas
    mapped_df = pd.read_csv(temp_tsv_path, sep='\t')
    
    # Save it out as a true CSV
    mapped_df.to_csv(final_csv_path, index=False)
    print(f"\nSuccess! Your fully optimized CSV is saved at: {final_csv_path}")
    print(f"Total overlapping targets found: {len(mapped_df) - 1}") # -1 for header
except Exception as e:
    print(f"Error during CSV conversion: {e}")

# 5. Clean up the temporary files
os.remove(pattern_file)
if os.path.exists(temp_tsv_path):
    os.remove(temp_tsv_path)