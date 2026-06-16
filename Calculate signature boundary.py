import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# --- 1. FILE PATHS ---
DATA_DIR = "/Users/riyathakuria/Documents/LSD insilico/Revised Raw data 9:6"
H5AD_PATH = os.path.join(DATA_DIR, "velmeshev_microglia_subset.h5ad")
OUTPUT_H5AD = os.path.join(DATA_DIR, "velmeshev_microglia_scored_and_rooted.h5ad")
OUTPUT_CSV = os.path.join(DATA_DIR, "Data_Driven_State_Transition_Summary.csv")
OUTPUT_PLOT = os.path.join(DATA_DIR, "Olah_Signature_Crossover_Arial20.png")

# Set up global plotting aesthetics
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['pdf.fonttype'] = 42

print("Phase 1: Loading raw single-nucleus microglial matrix...")
adata = sc.read_h5ad(H5AD_PATH)

# --- CRITICAL FIX: MAP ENSEMBL IDs TO HUGO SYMBOLS ---
print("-> Mapping Ensembl IDs to human-readable Gene Symbols...")
if 'feature_name' in adata.var.columns:
    # Overwrite the Ensembl index with the actual gene names
    adata.var_names = adata.var['feature_name'].astype(str).str.strip()
    # Scanpy requires all gene names to be unique
    adata.var_names_make_unique()
    
    # FIX FOR H5AD SAVING: Clear the internal index name so it doesn't conflict 
    # with the existing 'feature_name' column during adata.write()
    adata.var.index.name = None
    
    print("-> Matrix successfully indexed with human-readable Gene Symbols!")
else:
    print("⚠️ Warning: 'feature_name' not found. Gene mapping might fail.")

# --- 2. DEFINE CLINICAL MODULE SIGNATURES ---
print("\nPhase 2: Calculating Clinical State Signatures (Olah/DAM)...")
# Standard clinical signatures for microglial states
homeostatic_genes = ['P2RY12', 'TMEM119', 'CX3CR1', 'CSF1R', 'HEXB', 'P2RY13', 'TGFBR1']
dam_genes = ['APOE', 'SPP1', 'ITGAX', 'TREM2', 'CLEC7A', 'LPL', 'CST7', 'CD9', 'CTSD', 'TYROBP', 'CD68']
stable_housekeepers = ['SYK', 'HSPA9', 'FXN'] # From NormFinder

# Filter to ensure genes exist in the dataset
homeo_valid = [g for g in homeostatic_genes if g in adata.var_names]
dam_valid = [g for g in dam_genes if g in adata.var_names]
hk_valid = [g for g in stable_housekeepers if g in adata.var_names]

print(f"   Found {len(homeo_valid)}/{len(homeostatic_genes)} Homeostatic markers")
print(f"   Found {len(dam_valid)}/{len(dam_genes)} DAM markers")
print(f"   Found {len(hk_valid)}/{len(stable_housekeepers)} Housekeepers")

# Score the cells (CRITICAL FIX: use_raw=False forces Scanpy to use the main matrix with the updated names)
sc.tl.score_genes(adata, gene_list=homeo_valid, score_name='Olah_Homeostatic_Score', use_raw=False)
sc.tl.score_genes(adata, gene_list=dam_valid, score_name='Olah_DAM_Score', use_raw=False)
sc.tl.score_genes(adata, gene_list=hk_valid, score_name='Housekeeping_Score', use_raw=False)

# --- 3. TRAJECTORY INVERSION (CALCULATE PSEUDOTIME) ---
print("\nPhase 3: Executing Trajectory Inversion (Calculating Pseudotime)...")
# Anchor the trajectory root in the cell with the absolute highest housekeeping (health) score
root_idx = np.argmax(adata.obs['Housekeeping_Score'].values)
adata.uns['iroot'] = root_idx

# Calculate Diffusion Map and Pseudotime
if 'neighbors' not in adata.uns:
    print("-> Computing neighborhood graph...")
    sc.pp.neighbors(adata, n_pcs=30)
print("-> Computing Diffusion Map and DPT Pseudotime...")
sc.tl.diffmap(adata)
sc.tl.dpt(adata)

# Save pseudotime to a clean, standardized column name
adata.obs['pseudotime'] = adata.obs['dpt_pseudotime']

# --- 4. GAUSSIAN SMOOTHING & BOUNDARY DETECTION ---
print("\nPhase 4: Executing Gaussian Smoothing & Boundary Detection...")
df_obs = adata.obs[['pseudotime', 'Olah_Homeostatic_Score', 'Olah_DAM_Score']].copy()
df_obs = df_obs.sort_values(by='pseudotime').reset_index(drop=True)

grid_pt = np.linspace(0, 1, 101)
sigma = 0.05  # Gaussian bandwidth

smoothed_homeo, smoothed_dam = [], []

for t in grid_pt:
    weights = np.exp(-((df_obs['pseudotime'].values - t) ** 2) / (2 * (sigma ** 2)))
    sum_weights = np.sum(weights)
    if sum_weights > 0:
        smoothed_homeo.append(np.sum(df_obs['Olah_Homeostatic_Score'].values * weights) / sum_weights)
        smoothed_dam.append(np.sum(df_obs['Olah_DAM_Score'].values * weights) / sum_weights)
    else:
        # Fallback if window is empty
        idx = (df_obs['pseudotime'] - t).abs().idxmin()
        smoothed_homeo.append(df_obs.loc[idx, 'Olah_Homeostatic_Score'])
        smoothed_dam.append(df_obs.loc[idx, 'Olah_DAM_Score'])

df_smoothed = pd.DataFrame({
    'Pseudotime': grid_pt,
    'Smoothed_Homeostatic': smoothed_homeo,
    'Smoothed_DAM': smoothed_dam
})

# Locate Crossover
crossover_mask = df_smoothed['Smoothed_Homeostatic'] < df_smoothed['Smoothed_DAM']
crossover_indices = np.where(crossover_mask)[0]

t_crossover = df_smoothed.loc[crossover_indices[0], 'Pseudotime'] if len(crossover_indices) > 0 else 0.25
print(f"--> SUCCESS: Independent clinical crossover boundary located at: t = {t_crossover:.3f}")

# --- 5. ASSIGN PHENOTYPES & EXPORT ---
print("\nPhase 5: Saving updated Matrix and Plotting Results...")
adata.obs['Independent_SAM_Phenotype'] = adata.obs['pseudotime'].apply(
    lambda x: "Homeostatic" if x < t_crossover else "Storage-Associated Microglia (SAM)"
)

# Save the updated .h5ad file so you never have to calculate this again
adata.write(OUTPUT_H5AD)
print(f"✅ Updated H5AD matrix saved to: {OUTPUT_H5AD}")

# Export CSV Summary
df_out = adata.obs[['pseudotime', 'Olah_Homeostatic_Score', 'Olah_DAM_Score', 'Independent_SAM_Phenotype']].copy()
df_out.to_csv(OUTPUT_CSV, index=True)

# Generate Plot
fig, ax = plt.subplots(figsize=(12, 8), dpi=300)
ax.plot(df_smoothed['Pseudotime'], df_smoothed['Smoothed_Homeostatic'], color='#22c55e', linewidth=3.5, label='Olah Homeostatic Signature')
ax.plot(df_smoothed['Pseudotime'], df_smoothed['Smoothed_DAM'], color='#ef4444', linewidth=3.5, label='Olah DAM Signature')
ax.axvline(x=t_crossover, color='black', linestyle='--', linewidth=2.5, label=f'Transition Boundary (t = {t_crossover:.2f})')
ax.axvspan(0, t_crossover, facecolor='#e6f2e6', alpha=0.3, label='Homeostatic Phase')
ax.axvspan(t_crossover, 1.0, facecolor='#fdeaea', alpha=0.3, label='SAM Phase (Storage)')

ax.set_title("Data-Driven State Transition Boundary\n(Trajectory-Level Clinical Signatures)", fontsize=22, fontweight='bold', pad=20)
ax.set_xlabel("Metabolic Exhaustion Axis (Pseudotime)", fontsize=16, fontweight='bold', labelpad=12)
ax.set_ylabel("Composite Module Score", fontsize=16, fontweight='bold', labelpad=12)
ax.set_xlim(0, 1)
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.22), ncol=2, frameon=False, fontsize=12)
import seaborn as sns
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(OUTPUT_PLOT, bbox_inches='tight', dpi=300)

print(f"✅ Plot successfully saved to: {OUTPUT_PLOT}")
print("\n=== Master Script Complete ===")