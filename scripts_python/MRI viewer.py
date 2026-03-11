import pandas as pd
import numpy as np
import nilearn.plotting as plotting
from nilearn import datasets, image
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt

# 1. Interactive File Selection
root = tk.Tk()
root.withdraw()
print("Please select your 'Final_10_Targets_Summary copy.csv' file...")
file_path = filedialog.askopenfilename(title="Select your CSV file")

if not file_path:
    print("No file selected. Exiting.")
else:
    df = pd.read_csv(file_path)
    
    # 2. Fetch the standard Harvard-Oxford Atlases
    sub_atlas = datasets.fetch_atlas_harvard_oxford('sub-maxprob-thr25-2mm')
    cor_atlas = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm')
    
    sub_img = image.load_img(sub_atlas.maps)
    cor_img = image.load_img(cor_atlas.maps)

    # Professional Region Mapping (L/R combined for balance)
    region_mapping = {
        "Cerebellum": {"img": sub_img, "labels": [3, 14]},
        "Hippocampus": {"img": sub_img, "labels": [9, 19]},
        "MD_Thalamus": {"img": sub_img, "labels": [4, 15]},
        "Striatum": {"img": sub_img, "labels": [11, 21]},
        "Visual_Cortex": {"img": cor_img, "labels": [33]}
    }

    for _, row in df.iterrows():
        gene = row['Dynamic_Gene']
        val = row['Dynamic_Enrichment']
        reg = row['Region']
        
        if reg in region_mapping:
            mapping = region_mapping[reg]
            atlas_data = mapping['img'].get_fdata()
            
            # Create the 3D enrichment volume
            res_data = np.zeros(atlas_data.shape)
            for label in mapping['labels']:
                res_data[atlas_data == label] = val
            
            res_img = image.new_img_like(mapping['img'], res_data)

            # 3. Find ROI center for perfect slice framing
            # This ensures the 'ortho' view is centered exactly on the region
            coords = plotting.find_xyz_cut_coords(res_img)

            # 4. Generate High-Res Plot (Magma Colormap)
            fig_name = f"{gene}_{reg}_Tier1_Plot"
            
            # Use 'plot_stat_map' for the professional enrichment look
            display = plotting.plot_stat_map(
                res_img,
                title=f"{gene} Enrichment: {reg}",
                display_mode='ortho',
                cut_coords=coords,
                cmap='magma',        # High-impact standard
                colorbar=True,
                annotate=True,
                draw_cross=False,
                black_bg=True,
                threshold=0.001      # Removes noise
            )
            
            # Save both PNG (for draft) and PDF (for publication)
            display.savefig(f"{fig_name}.png", dpi=300)
            display.savefig(f"{fig_name}.pdf")
            plt.close()

            # 5. Save NIfTI for FSLeyes manual polish
            res_img.to_filename(f"{gene}_volume.nii.gz")
            print(f"Generated files for: {gene}")

    print("\nAll files are ready in your folder.")