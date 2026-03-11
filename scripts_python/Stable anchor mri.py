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
    # Opens file dialog to select the Stable_Anchors_Intersection_AllRegions.csv file
    path = filedialog.askopenfilename(title="Select Stable Anchors CSV")
    root.destroy()
    return path

def plot_psap_markers(csv_path):
    # 1. Load and Filter specifically for PSAP
    df = pd.read_csv(csv_path)
    # Filter for the specific gene 'PSAP'
    df_psap = df[df['Gene_Symbol'] == 'PSAP'].copy()
    
    # 2. Coordinates Mapping
    # Matching the CSV region names: Cerebellum, Hippocampus, Striatum, Visual
    coords_map = {
        'Cerebellum': [0, -60, -20],
        'Hippocampus': [25, -20, -15],
        'Striatum': [20, 10, 0],
        'Visual': [15, -85, 5]
    }
    df_psap['coords'] = df_psap['Region'].map(coords_map)
    coords_list = list(df_psap['coords'])

    # 3. Layout: Brain Visualization and Legends
    fig = plt.figure(figsize=(20, 10))
    gs = fig.add_gridspec(2, 2, width_ratios=[12, 8], height_ratios=[1, 1])
    
    # Using a distinct color palette for the 4 regions
    colors_list = plt.cm.Dark2(np.linspace(0, 1, len(df_psap)))
    custom_cmap = mcolors.ListedColormap(colors_list)
    
    # Scale for marker size (Enrichment Score E)
    scale_factor = 600 
    sizes_psap = df_psap['Enrichment_Score_E'] * scale_factor

    # 4. Plot Brain Markers (Ortho view)
    ax_brain = fig.add_subplot(gs[:, 0])
    display = plotting.plot_markers(
        node_values=range(len(df_psap)),
        node_coords=coords_list,
        node_cmap=custom_cmap,
        node_size=sizes_psap,
        display_mode='ortho',
        title="Stable Anchor Distribution: PSAP",
        colorbar=False,
        figure=fig,
        axes=ax_brain
    )

    # 5. Top Right: Region/Gene Legend
    ax_leg1 = fig.add_subplot(gs[0, 1])
    ax_leg1.axis('off')
    
    id_elements = [
        Line2D([0], [0], marker='o', color='w', 
               label=f"PSAP ({row['Region']})",
               markerfacecolor=colors_list[i], markersize=14)
        for i, row in df_psap.reset_index().iterrows()
    ]
    ax_leg1.legend(handles=id_elements, loc='center left', title="Gene & Region", 
                   fontsize=12, title_fontsize=14, labelspacing=1.8, frameon=False)

    # 6. Bottom Right: Scale Legend (Enrichment Score E)
    ax_leg2 = fig.add_subplot(gs[1, 1])
    ax_leg2.axis('off')
    
    size_refs = [0.95, 1.0, 1.05] # Range based on PSAP data (approx 0.97 to 1.03)
    size_elements = [
        Line2D([0], [0], marker='o', color='w', label=f"E = {s}",
               markerfacecolor='gray', markersize=np.sqrt(s * scale_factor)) 
        for s in size_refs
    ]
    ax_leg2.legend(handles=size_elements, loc='center left', title="Enrichment Score (E)", 
                   fontsize=12, title_fontsize=14, labelspacing=2.0, frameon=False)

    plt.tight_layout()
    output_path = "PSAP_Stable_Anchors_Map.png"
    plt.savefig(output_path, dpi=300)
    print(f"Publication-ready figure saved to: {output_path}")
    plt.show()

if __name__ == "__main__":
    path = select_file()
    if path:
        plot_psap_markers(path)