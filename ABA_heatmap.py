import pandas as pd
import numpy as np
from scipy.stats import zscore
from scipy.cluster.hierarchy import linkage, leaves_list
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches

# 1. Define configurations
file_path = r"/Users/riyathakuria/Documents/LSD_Insilico_Analysis_24:3:26/Raw data/PPI_SpatioTemporal_Expression.csv"
output_plot = r"/Users/riyathakuria/Documents/LSD_Insilico_Analysis_24:3:26/Raw data/PPI_Spatiotemporal_Heatmap_ThickBars.png"

chronological_ages = [
    '8 pcw', '9 pcw', '12 pcw', '13 pcw', '16 pcw', '17 pcw', '19 pcw', '21 pcw', 
    '24 pcw', '25 pcw', '26 pcw', '35 pcw', '37 pcw', '4 mos', '10 mos', '1 yrs', 
    '2 yrs', '3 yrs', '4 yrs', '8 yrs', '11 yrs'
]

regions = ['Cerebellum', 'Hippocampus', 'Striatum', 'Visual Cortex']

region_colors = {
    'Cerebellum': '#d09292',
    'Hippocampus': '#a4b6dd',
    'Striatum': '#a2d0c0',
    'Visual Cortex': '#c37892'
}

# 2. Load and preprocess the data
df = pd.read_csv(file_path)
agg = df.groupby(['Gene', 'Lineage', 'Age'])['Expression'].mean().reset_index()
pivot = agg.pivot(index='Gene', columns=['Lineage', 'Age'], values='Expression')

# Reorder columns
sorted_cols = []
for r in regions:
    for a in chronological_ages:
        if (r, a) in pivot.columns:
            sorted_cols.append((r, a))
pivot = pivot[sorted_cols]

# Fill missing and Z-score
pivot = pivot.apply(lambda row: row.fillna(row.mean()), axis=1)
z_matrix_array = zscore(pivot.values, axis=1)
z_matrix = pd.DataFrame(z_matrix_array, index=pivot.index, columns=pivot.columns)

# 3. Hierarchical clustering
Z = linkage(z_matrix, method='ward')
idx = leaves_list(Z)
z_matrix_sorted = z_matrix.iloc[idx]

# 4. Color Maps
expr_cmap = plt.cm.RdBu_r 
age_cmap = plt.cm.viridis
age_mapping = {age: age_cmap(i / (len(chronological_ages)-1)) for i, age in enumerate(chronological_ages)}

# 5. Build the plot architecture
fig = plt.figure(figsize=(24, 20))
widths = [sum(c[0] == r for c in z_matrix.columns) for r in regions] + [4] 

# UPDATED: height_ratios increased from 0.05 to 0.3 for thicker bars
gs = fig.add_gridspec(3, 5, width_ratios=widths, height_ratios=[0.3, 0.3, 10], wspace=0.1, hspace=0.02)

# 6. Plot panels
for i, r in enumerate(regions):
    r_cols = [c for c in z_matrix_sorted.columns if c[0] == r]
    data = z_matrix_sorted[r_cols]
    
    # Top Row: Lineage Color Bar (THICKER)
    ax_lin = fig.add_subplot(gs[0, i])
    ax_lin.imshow([[1]*len(r_cols)], aspect='auto', 
                  cmap=LinearSegmentedColormap.from_list("", [region_colors[r], region_colors[r]]))
    ax_lin.axis('off')
    ax_lin.set_title(r, fontsize=16, fontweight='bold', pad=10)
    
    # Middle Row: Age Color Bar (THICKER)
    ax_age = fig.add_subplot(gs[1, i])
    age_color_matrix = [[age_mapping[c[1]] for c in r_cols]]
    ax_age.imshow(age_color_matrix, aspect='auto')
    ax_age.axis('off')
    
    # Bottom Row: Heatmap
    ax_hm = fig.add_subplot(gs[2, i])
    im = ax_hm.imshow(data, aspect='auto', cmap=expr_cmap, vmin=-2.5, vmax=2.5)
    ax_hm.set_xticks([])
    if i == 0:
        ax_hm.set_yticks(range(len(z_matrix_sorted)))
        ax_hm.set_yticklabels(z_matrix_sorted.index, fontsize=9)
    else:
        ax_hm.set_yticks([])

# 7. Legends
gs_leg = gs[2, 4].subgridspec(3, 1, height_ratios=[1, 4, 1], hspace=0.4)
ax_leg_region = fig.add_subplot(gs_leg[0])
ax_leg_region.axis('off')
patches = [mpatches.Patch(color=region_colors[r], label=r) for r in regions]
ax_leg_region.legend(handles=patches, title="Brain Region", loc='upper left', frameon=False)

ax_leg_age = fig.add_subplot(gs_leg[1])
ax_leg_age.axis('off')
cax_age = ax_leg_age.inset_axes([0.1, 0.05, 0.15, 0.9])
cb_age = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=len(chronological_ages)-1), cmap=age_cmap), cax=cax_age)
cb_age.set_ticks(range(len(chronological_ages)))
cb_age.set_ticklabels(chronological_ages, fontsize=11)

ax_leg_expr = fig.add_subplot(gs_leg[2])
ax_leg_expr.axis('off')
cax_expr = ax_leg_expr.inset_axes([0.0, 0.5, 0.9, 0.2])
fig.colorbar(im, cax=cax_expr, orientation='horizontal')

# 8. Save
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.show()