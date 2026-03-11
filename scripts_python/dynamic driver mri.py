import pandas as pd
import matplotlib.pyplot as plt
from nilearn import plotting
import tkinter as tk
from tkinter import filedialog
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import numpy as np

def select_file():
    root = tk.Tk()
    root.withdraw()
    path = filedialog.askopenfilename(title="Select Top3 Dynamic Drivers CSV")
    root.destroy()
    return path

def plot_top1_dynamic_genes(csv_path):
    # 1. Load and Filter for Top 1 Gene per Region
    df = pd.read_csv(csv_path)
    
    # Sort by Enrichment Score and keep only the top entry for each region
    df_top1 = df.sort_values('Enrichment_Score_E', ascending=False).drop_duplicates('Region').copy()
    
    # 2. Coordinates Mapping (Adjusted 'Visual' to match your CSV name)
    coords_map = {
        'Cerebellum': [0, -60, -20],
        'Hippocampus': [25, -20, -15],
        'Striatum': [20, 10, 0],
        'Visual': [15, -85, 5] 
    }
    df_top1['coords'] = df_top1['Region'].map(coords_map)
    coords_list = list(df_top1['coords'])

    # 3. Layout Configuration
    fig = plt.figure(figsize=(20, 10))
    gs = fig.add_gridspec(2, 2, width_ratios=[12, 8], height_ratios=[1, 1])
    
    colors_list = plt.cm.Set1(np.linspace(0, 1, len(df_top1)))
    custom_cmap = mcolors.ListedColormap(colors_list)
    scale_factor = 600 # Increased for better visibility
    sizes_dyn = df_top1['Enrichment_Score_E'] * scale_factor

    # 4. Plot Brain Markers
    ax_brain = fig.add_subplot(gs[:, 0])
    display = plotting.plot_markers(
        node_values=range(len(df_top1)),
        node_coords=coords_list,
        node_cmap=custom_cmap,
        node_size=sizes_dyn,
        display_mode='ortho',
        title="Top 1 Dynamic Gene Enrichment per Region",
        colorbar=False,
        figure=fig,
        axes=ax_brain
    )

    # 5. Legend 1: Gene Identity
    ax_leg1 = fig.add_subplot(gs[0, 1])
    ax_leg1.axis('off')
    
    id_elements = [
        Line2D([0], [0], marker='o', color='w', 
               label=f"{row['Gene_Symbol']} ({row['Region']})",
               markerfacecolor=colors_list[i], markersize=14)
        for i, row in df_top1.reset_index().iterrows()
    ]
    ax_leg1.legend(handles=id_elements, loc='center left', title="Top Dynamic Gene", 
                   fontsize=12, title_fontsize=14, labelspacing=1.8, frameon=False)

    # 6. Legend 2: Enrichment Scale
    ax_leg2 = fig.add_subplot(gs[1, 1])
    ax_leg2.axis('off')
    
    size_refs = [1.0, 1.2, 1.4]
    size_elements = [
        Line2D([0], [0], marker='o', color='w', label=f"E = {s}",
               markerfacecolor='gray', markersize=np.sqrt(s * scale_factor)) 
        for s in size_refs
    ]
    ax_leg2.legend(handles=size_elements, loc='center left', title="Enrichment Score (E)", 
                   fontsize=12, title_fontsize=14, labelspacing=2.0, frameon=False)

    plt.tight_layout()
    output_path = "Top1_Dynamic_Genes_Map.png"
    plt.savefig(output_path, dpi=300)
    print(f"Figure saved to: {output_path}")
    plt.show()

if __name__ == "__main__":
    path = select_file()
    if path:
        plot_top1_dynamic_genes(path)