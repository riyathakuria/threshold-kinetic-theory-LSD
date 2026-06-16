import pandas as pd
import numpy as np
import os

# --- 1. CONFIGURATION & FILE PATHS ---
DATA_DIR = "/Users/riyathakuria/Documents/LSD insilico/Raw data 2"  # Update path if necessary
INPUT_CSV = os.path.join(DATA_DIR, "Normalized_Dynamic_Biomarkers.csv")
OUTPUT_CSV = os.path.join(DATA_DIR, "Ranked_Enrichment_Scores.csv")

# Safeguard to prevent unstable, inflated ratios from very sparse background expression
PSEUDOCOUNT = 0.01 

# Defensible minimum abundance threshold (Replaces the arbitrary 0.2% percentile)
# Acts as a practical noise floor to drop uninformative transcriptomic "dust"
MIN_ABUNDANCE_THRESHOLD = 0.05 

print("Phase 1: Loading NormFinder-Anchored Dynamic Biomarkers...")
df = pd.read_csv(INPUT_CSV)

# Extract 'Region' from the 'Sample_ID' column (e.g., "1 yrs_Cerebellum" -> "Cerebellum")
# This assumes Sample_ID format is "Age_Region"
df['Region'] = df['Sample_ID'].apply(lambda x: x.split('_')[1] if '_' in str(x) else 'Unknown')

# --- 2. THE TECHNICAL NOISE FLOOR FILTER ---
print("\nPhase 2: Applying Minimum Abundance Filter...")
# This is NOT a biological cutoff. It is a technical noise floor used to 
# reduce sparse background expression before ranking. We use a practical 
# abundance threshold rather than a percentile to ensure it doesn't collapse to zero.
print(f"-> Technical Noise Floor established at normalized expression: {MIN_ABUNDANCE_THRESHOLD}")

# Filter out sparse background expression to clean the data before ranking
initial_rows = len(df)
df_clean = df[df['Normalized_Expression'] >= MIN_ABUNDANCE_THRESHOLD].copy()
filtered_rows = len(df_clean)
print(f"-> Filtered out {initial_rows - filtered_rows} sparse background artifacts.")

# --- 3. CALCULATE LOCAL AND GLOBAL MEANS ---
print("\nPhase 3: Calculating Local and Global Expression Means...")
# Calculate the Local Mean (Specific Gene in a Specific Region during a Specific Stage)
local_means = df_clean.groupby(['Region', 'Stage', 'Gene'])['Normalized_Expression'].mean().reset_index()
local_means.rename(columns={'Normalized_Expression': 'Local_Mean_Expression'}, inplace=True)

# Calculate the Global Mean (Specific Gene averaged across the entire brain and lifespan)
global_means = df_clean.groupby('Gene')['Normalized_Expression'].mean().reset_index()
global_means.rename(columns={'Normalized_Expression': 'Global_Mean_Expression'}, inplace=True)

# --- 4. CALCULATE ENRICHMENT SCORES & RANK ---
print("\nPhase 4: Computing Fold-Enrichment Scores with Pseudocount Safeguard...")
# Merge local and global means
enrichment_df = pd.merge(local_means, global_means, on='Gene')

# Enrichment Score Formula: (Local Mean + Pseudocount) / (Global Mean + Pseudocount)
# The pseudocount prevents massive ratio inflation for genes with near-zero global means.
enrichment_df['Enrichment_Score'] = (enrichment_df['Local_Mean_Expression'] + PSEUDOCOUNT) / (enrichment_df['Global_Mean_Expression'] + PSEUDOCOUNT)

# Sort and Rank the data
# Enrichment Score acts as the primary discovery metric. 
# Local Mean Expression is utilized secondarily as a tie-breaker.
ranked_df = enrichment_df.sort_values(
    by=['Region', 'Stage', 'Enrichment_Score', 'Local_Mean_Expression'], 
    ascending=[True, True, False, False]
)

# --- 5. EXPORT FINAL RANKINGS ---
print("\nPhase 5: Exporting Ranked Matrix...")
ranked_df.to_csv(OUTPUT_CSV, index=False)
print(f"✅ Enrichment Ranking successfully saved to:\n   {OUTPUT_CSV}")

# --- PREVIEW THE TOP RESULTS ---
print("\n=== Quick Preview: Top 3 Enriched Genes per Region (Postnatal) ===")
preview_df = ranked_df[ranked_df['Stage'] == 'Postnatal'].groupby('Region').head(3)
print(preview_df[['Region', 'Gene', 'Enrichment_Score', 'Local_Mean_Expression']].to_string(index=False))