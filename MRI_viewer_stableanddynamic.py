import pandas as pd
import matplotlib.pyplot as plt
from nilearn import plotting
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import numpy as np

def generate_maps():
    print("Loading datasets...")
    path_dynamic = "/Users/riyathakuria/Documents/LSD_Insilico_Analysis_24:3:26/Raw data/Top_Dynamic_Drivers_Per_Region.csv"
    path_stable = "/Users/riyathakuria/Documents/LSD_Insilico_Analysis_24:3:26/Raw data/Top_Stable_Anchors_Per_Region.csv"
    
    df_dyn = pd.read_csv(path_dynamic)
    df_sta = pd.read_csv(path_stable)

    # 1. Filter for the Top 1 gene per region
    # For dynamic: Sort by Enrichment_Score
    df_dyn_top1 = df_dyn.sort_values('Enrichment_Score', ascending=False).drop_duplicates('Lineage').copy()
    
    # For stable: Filter by Rank 1.0
    df_sta_top1 = df_sta[df_sta['Rank'] == 1.0].drop_duplicates('Lineage').copy()

    # 2. Coordinates Mapping
    coords_map = {
        'Cerebellum': [0, -60, -20],
        'Hippocampus': [25, -20, -15],
        'Striatum': [20, 10, 0],
        'Visual Cortex': [15, -85, 5],
        'Visual': [15, -85, 5] # Fallback
    }
    
    df_dyn_top1['coords'] = df_dyn_top1['Lineage'].map(coords_map)
    df_sta_top1['coords'] = df_sta_top1['Lineage'].map(coords_map)

    def plot_working_map(df, value_col, title_str, scale_factor, legend_title, output_filename):
        # 3. Layout Configuration - INCREASED CANVAS SIZE (24, 12)
        fig = plt.figure(figsize=(24, 12), facecolor='white')
        gs = fig.add_gridspec(2, 2, width_ratios=[12, 8], height_ratios=[1, 1])
        
        colors_list = plt.cm.Set1(np.linspace(0, 1, len(df)))
        custom_cmap = mcolors.ListedColormap(colors_list)
        
        coords_list = list(df['coords'])
        node_sizes = df[value_col] * scale_factor

        # 4. Plot Brain Markers
        ax_brain = fig.add_subplot(gs[:, 0])
        ax_brain.set_facecolor('white')
        
        display = plotting.plot_markers(
            node_values=range(len(df)),
            node_coords=coords_list,
            node_cmap=custom_cmap,
            node_size=node_sizes,
            display_mode='ortho',
            title=title_str,
            colorbar=False,
            figure=fig,
            axes=ax_brain
        )

        # 5. Legend 1: Gene Identity
        ax_leg1 = fig.add_subplot(gs[0, 1])
        ax_leg1.axis('off')
        
        id_elements = [
            Line2D([0], [0], marker='o', color='w', 
                   label=f"{row['Gene']} ({row['Lineage']})",
                   markerfacecolor=colors_list[i], markersize=14)
            for i, row in df.reset_index().iterrows()
        ]
        ax_leg1.legend(handles=id_elements, loc='center left', title="Region's Top Gene", 
                       fontsize=14, title_fontsize=16, labelspacing=1.8, frameon=False)

        # 6. Legend 2: Size Scale
        ax_leg2 = fig.add_subplot(gs[1, 1])
        ax_leg2.axis('off')
        
        val_min, val_max = df[value_col].min(), df[value_col].max()
        val_mid = (val_min + val_max) / 2
        size_refs = [val_max, val_mid, val_min]
        
        size_elements = [
            Line2D([0], [0], marker='o', color='w', label=f"{s:.2f}",
                   markerfacecolor='gray', markersize=np.sqrt(s * scale_factor)) 
            for s in size_refs
        ]
        
        # INCREASED LABELSPACING to 3.5 to prevent bubble overlap
        ax_leg2.legend(handles=size_elements, loc='center left', title=legend_title, 
                       fontsize=14, title_fontsize=16, labelspacing=3.5, borderpad=1.5, frameon=False)

        plt.tight_layout()
        
        # Save as High-Res PNG
        plt.savefig(output_filename, dpi=600, bbox_inches='tight', facecolor='white')
        print(f"Figure saved successfully to: {output_filename}")
        plt.close(fig)

    # Dynamic Drivers Plot
    plot_working_map(
        df=df_dyn_top1, 
        value_col='Enrichment_Score', 
        title_str="Top 1 Dynamic Driver Per Region", 
        scale_factor=200, 
        legend_title="Enrichment Score", 
        output_filename="Top1_Dynamic_Drivers_Map.png"
    )

    # Stable Anchors Plot - SCALE FACTOR REDUCED TO 2
    plot_working_map(
        df=df_sta_top1, 
        value_col='Mean_Expression', 
        title_str="Top 1 Stable Anchor Per Region", 
        scale_factor=2, # <--- Change this up or down slightly if the bubbles need fine-tuning
        legend_title="Mean Expression", 
        output_filename="Top1_Stable_Anchors_Map.png"
    )

if __name__ == "__main__":
    generate_maps()