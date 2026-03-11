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
    path = filedialog.askopenfilename(title="Select Region Metrics CSV")
    root.destroy()
    return path

def plot_gene_markers_clean_final(csv_path):
    # 1. Load and Filter
    df = pd.read_csv(csv_path)
    df = df[df['Region'] != 'Dorsal Thalamus'].copy().reset_index(drop=True)
    
    # 2. Coordinates
    coords_map = {
        'Cerebellum': [0, -60, -20],
        'Hippocampus': [25, -20, -15],
        'Striatum': [20, 10, 0],
        'Visual Cortex': [15, -85, 5]
    }
    df['coords'] = df['Region'].map(coords_map)
    coords_list = list(df['coords'])

    # ---------------------------------------------------------
    # LAYOUT: Use GridSpec to create a dedicated sidebar column
    # ---------------------------------------------------------
    fig = plt.figure(figsize=(24, 10))
    # 3 columns: Brain (width 15), Spacer (width 1), Legends (width 8)
    gs = fig.add_gridspec(2, 3, width_ratios=[15, 1, 8], height_ratios=[1, 1])
    
    # Define colors and sizes
    colors_list = plt.cm.tab10(np.linspace(0, 1, len(df)))
    custom_cmap = mcolors.ListedColormap(colors_list)
    scale_factor = 500
    sizes_dyn = df['Dynamic_Enrichment_E'] * scale_factor

    # 3. Plot Brain on the left (spanning both rows)
    ax_brain = fig.add_subplot(gs[:, 0])
    display = plotting.plot_markers(
        node_values=range(len(df)),
        node_coords=coords_list,
        node_cmap=custom_cmap,
        node_size=sizes_dyn,
        display_mode='ortho',
        title="Dynamic Genes: Regional Enrichment Distribution",
        colorbar=False,
        figure=fig,
        axes=ax_brain
    )

    # 4. Use the RIGHT-MOST column for Legends
    # Identity Legend (Top Right)
    ax_leg1 = fig.add_subplot(gs[0, 2])
    ax_leg1.axis('off') # Hide the axis lines
    
    id_elements = [
        Line2D([0], [0], marker='o', color='w', label=f"{row['Top_Dynamic_Gene']} ({row['Region']})",
               markerfacecolor=colors_list[i], markersize=14)
        for i, row in df.iterrows()
    ]
    ax_leg1.legend(handles=id_elements, loc='upper left', title="Gene Identity", 
                   fontsize=12, title_fontsize=14, labelspacing=2.0, frameon=False)

    # Scale Legend (Bottom Right)
    ax_leg2 = fig.add_subplot(gs[1, 2])
    ax_leg2.axis('off')
    
    size_refs = [1.5, 2.0, 2.5]
    size_elements = [
        Line2D([0], [0], marker='o', color='w', label=f"E = {s}",
               markerfacecolor='gray', markersize=np.sqrt(s * scale_factor)) 
        for s in size_refs
    ]
    ax_leg2.legend(handles=size_elements, loc='upper left', title="Enrichment Scale", 
                   fontsize=12, title_fontsize=14, labelspacing=2.5, frameon=False)

    output_path = "Figure_2_Dynamic_NoOverlap.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Final publication-ready figure saved to: {output_path}")
    plt.show()

if __name__ == "__main__":
    path = select_file()
    if path:
        plot_gene_markers_clean_final(path)