import subprocess

# 1. Define the path to your 20GB raw dataset
large_file_path = "/Users/riyathakuria/Documents/LSD_Insilico_Analysis_24:3:26/240526_JH_NP_28_LysoIP_Fcrls_Microglia_candidates.tsv"

# 2. The Shell Pipeline
# - awk -F'\t': Reads the file using Tabs as separators (change to ',' if it's a true CSV)
# - '{print $10}': Extracts the 10th column (the 'Genes' column in your dataset structure)
# - tail -n +2: Skips the header row
# - sort -u: Alphabetizes and keeps only UNIQUE entries
# - wc -l: Counts the final number of lines/proteins
shell_command = f"awk -F'\\t' '{{print $10}}' '{large_file_path}' | tail -n +2 | sort -u | wc -l"

print("Scanning the 20GB Brain Lysosome Atlas...")
print("(This may take a minute as it streams through the hard drive)")

try:
    # Execute the shell stream directly from Python
    result = subprocess.run(shell_command, shell=True, text=True, capture_output=True)
    
    # Extract and clean the output number
    unique_count = result.stdout.strip()
    
    print(f"\nSuccess! Total unique proteins in the Brain Lysosome Atlas: {unique_count}")

except Exception as e:
    print(f"\nError executing command: {e}")