import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tkinter import Tk, filedialog
import sys

def generate_aba_custom_plots():
    # 1. Pop-up to select file
    root = Tk()
    root.withdraw()
    root.attributes('-topmost', True)
    file_path = filedialog.askopenfilename(
        title="Select your Long_Format_ABA_Data CSV file",
        filetypes=[("CSV files", "*.csv")]
    )
    root.destroy()

    if not file_path:
        print("No file selected. Exiting.")
        sys.exit()

    # 2. Load Data
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error loading file: {e}")
        return

    # 3. Define Stage Order and Custom Hex Colors
    # Bottom to top order for stacking
    stage_color_map = {
        'Early Prenatal (8-12 pcw)': '#191F31',
        'Mid Prenatal (13-24 pcw)': '#365359',
        'Late Prenatal (25-38 pcw)': '#5F707F',
        'Infancy (0-12 mos)': '#CD853F',
        'Early Childhood (1-5 yrs)': '#DAA520',
        'Late Childhood (6-11 yrs)': '#8B4513'
    }
    
    stage_order = list(stage_color_map.keys())
    colors = [stage_color_map[stage] for stage in stage_order]

    # Ensure categorical ordering
    df['Developmental_Stage'] = pd.Categorical(df['Developmental_Stage'], categories=stage_order, ordered=True)

    # 4. Data Transformation: 100% Normalization
    gene_totals = df.groupby('Gene_Symbol')['Expression_Value'].transform('sum')
    df['Relative_Contribution'] = df['Expression_Value'] / gene_totals

    # 5. Plotting Function
    def plot_stacked_bars(data, classification_type):
        # Sort genes alphabetically
        sorted_genes = sorted(data['Gene_Symbol'].unique())
        
        # Pivot for plotting
        pivot_df = data.pivot_table(
            index='Gene_Symbol', 
            columns='Developmental_Stage', 
            values='Relative_Contribution'
        ).reindex(sorted_genes)

        # Plotting
        fig, ax = plt.subplots(figsize=(16, 9), dpi=300)
        
        pivot_df.plot(
            kind='bar', 
            stacked=True, 
            ax=ax, 
            color=colors, 
            width=0.8,
            edgecolor='white',
            linewidth=0.3
        )

        # Aesthetics
        ax.set_title(f'100% Stacked Bar Chart: {classification_type} Genes', fontsize=20, fontweight='bold', pad=25)
        ax.set_ylabel('Relative Expression Contribution', fontsize=14, fontweight='bold')
        ax.set_xlabel('Gene Symbol', fontsize=14, fontweight='bold')
        ax.set_ylim(0, 1.0)
        
        # Grid and Spines
        ax.yaxis.set_major_locator(plt.FixedLocator([0.25, 0.50, 0.75, 1.0]))
        ax.yaxis.grid(True, linestyle='--', alpha=0.3, color='grey')
        ax.set_axisbelow(True)
        sns.despine(ax=ax, top=True, right=True)

        # X-Axis Labels Rotation
        rotation = 90 if classification_type == "Stable" else 45
        plt.xticks(rotation=rotation, ha='right', fontsize=10)
        plt.yticks(fontsize=12)

        # Legend
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            reversed(handles), reversed(labels), # Reverse to match stacking order (top label = top segment)
            title="Developmental Stage", 
            bbox_to_anchor=(1.02, 1), 
            loc='upper left', 
            fontsize=11, 
            title_fontsize=13,
            frameon=False
        )

        plt.tight_layout()
        filename = f"{classification_type}_Genes_Custom_Colors.png"
        plt.savefig(filename, bbox_inches='tight')
        print(f"Saved: {filename}")
        plt.show()

    # 6. Generate Separated Plots
    dynamic_genes = df[df['Classification'] == 'Dynamic']
    stable_genes = df[df['Classification'] == 'Stable']

    if not dynamic_genes.empty:
        plot_stacked_bars(dynamic_genes, "Dynamic")
    
    if not stable_genes.empty:
        plot_stacked_bars(stable_genes, "Stable")

if __name__ == "__main__":
    generate_aba_custom_plots()