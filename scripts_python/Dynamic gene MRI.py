import pandas as pd
import numpy as np
import nilearn.plotting as plotting
from nilearn import datasets, image
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import tkinter as tk
from tkinter import filedialog

# 1. Interactive File Selection
root = tk.Tk()
root.withdraw()
print("Please select your 'Final_10_Targets_Summary copy.csv' file...")
file_path = filedialog.askopenfilename(title="Select Dynamic Genes CSV")

if not file_path:
    print("No file selected. Exiting.")
else:
    df = pd.read_csv(file_path)
    
    # 2. Setup Atlases (Standard MNI152)
    sub_atlas = datasets.fetch_atlas_harvard_oxford('sub-maxprob-thr25-2mm')
    cor_atlas = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm')
    sub_img = image.load_img(sub_atlas.maps)
    cor_img = image.load_img(cor_atlas.maps)

    region_info = {
        "Cerebellum": {"img": sub_img, "labels": [3, 14]},
        "Hippocampus": {"img": sub_img, "labels": [9, 19]},
        "MD_Thalamus": {"img": sub_img, "labels": [4, 15]},
        "Striatum": {"img": sub_img, "labels": [11, 21]},
        "Visual_Cortex": {"img": cor_img, "labels": [33]}
    }

    # 3. Generate 5 Separate Maps with Unified Scale (vmax=2.5)
    for _, row in df.iterrows():
        gene = row['Dynamic_Gene']
        val = row['Dynamic_Enrichment']
        reg = row['Region']
        
        # Create individual gene volume
        vol_data = np.zeros(sub_img.shape)
        if reg in region_info:
            info = region_info[reg]
            atlas_data = info['img'].get_fdata()
            for label_idx in info['labels']:
                vol_data[atlas_data == label_idx] = val
                
        gene_img = image.new_img_like(sub_img, vol_data)
        
        # Focus slices on the specific target region
        coords = plotting.find_xyz_cut_coords(gene_img)
        
        # Generate Figure with Tier 1 Title and Scaling
        display = plotting.plot_stat_map(
            gene_img, 
            title=f"Anatomical Enrichment of {gene} in {reg}",
            display_mode='ortho', 
            cut_coords=coords,
            cmap='magma', 
            colorbar=True,
            vmax=2.5,            # UNIFIED SCALE: Crucial for comparison
            annotate=True, 
            draw_cross=False, 
            black_bg=True,
            threshold=0.01
        )
        
        # Format Scale: Changes scientific notation to clean decimals (e.g., 2.5)
        display._cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
        # Save Outputs
        display.savefig(f"Figure_Dynamic_{gene}_Unified.png", dpi=300)
        display.savefig(f"Figure_Dynamic_{gene}_Unified.pdf")
        plt.close()
        
        print(f"Generated Figure for {gene}")

    print("\nAll dynamic gene maps are ready with unified scaling.")