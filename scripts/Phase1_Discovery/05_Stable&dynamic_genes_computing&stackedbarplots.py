import pandas as pd
import numpy as np
from scipy.stats import chisquare
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Configuration and File Path
ROOT = Path(__file__).resolve().parents[2]

# 2. Load Data
df_raw = pd.read_csv(file_path)

# 3. Pre-processing: Map ages and average across brain regions
def map_age_to_stage(age_str):
    if 'pcw' in age_str:
        return 'Prenatal'
    elif 'mos' in age_str:
        return 'Postnatal'
    elif 'yrs' in age_str:
        years = int(age_str.split(' ')[0])
        return 'Early Childhood' if years <= 4 else 'Late Childhood'
    return 'Unknown'

df_raw['Age_Category'] = df_raw['Age'].apply(map_age_to_stage)

# Average expression per Gene and Stage (removing Lineage)
df = df_raw.groupby(['Gene', 'Age_Category'])['Expression'].mean().reset_index()
df.columns = ['Gene', 'Age_Category', 'Mean']

# 4. Metric Calculation
gene_totals = df.groupby('Gene')['Mean'].sum().reset_index()
gene_totals.columns = ['Gene', 'Total_Mean']
df = df.merge(gene_totals, on='Gene')

# Percentage of total expression for each stage
df['Stage_Percent'] = (df['Mean'] / df['Total_Mean']) * 100
df['Abs_Dev'] = (df['Stage_Percent'] - 25).abs()

# Inferential Statistical Analysis (Chi-square)
def calc_metrics(group):
    # Skip genes with zero total expression to avoid division errors
    total = group['Mean'].sum()
    if total == 0:
        return pd.Series({'Avg_Dev': 0, 'p_value': 1.0})
    
    avg_dev = group['Abs_Dev'].mean()
    observed = group['Mean'].values
    # Expected: Total sum split equally across 4 stages (25% each)
    expected = [total / 4] * 4
    _, p_val = chisquare(observed, f_exp=expected)
    return pd.Series({'Avg_Dev': avg_dev, 'p_value': p_val})

metrics_df = df.groupby('Gene').apply(calc_metrics).reset_index()

# 5. Apply Dual-Criteria Classification
def classify_genes(row):
    # Stable: Close to 25% across stages (<=6% dev) AND p > 0.05 (no significant flux)
    if row['Avg_Dev'] <= 6 and row['p_value'] > 0.05:
        return 'Stable'
    else:
        return 'Dynamic'

metrics_df['Stability_Status'] = metrics_df.apply(classify_genes, axis=1)

# Save the final classification table
metrics_df.to_csv('Gene_Stability_Classification.csv', index=False)

# 6. Visualization
plot_data = df.merge(metrics_df[['Gene', 'Stability_Status']], on='Gene')
plot_data['Relative_Contribution'] = plot_data['Stage_Percent'] / 100

stage_order = ['Prenatal', 'Postnatal', 'Early Childhood', 'Late Childhood']
stage_colors = {
    'Prenatal': '#191F31', 
    'Postnatal': '#CD853F', 
    'Early Childhood': '#DAA520', 
    'Late Childhood': '#8B4513'
}
colors = [stage_colors[s] for s in stage_order]
plot_data['Age_Category'] = pd.Categorical(plot_data['Age_Category'], categories=stage_order, ordered=True)

def generate_plot(status, filename):
    subset = plot_data[plot_data['Stability_Status'] == status]
    if subset.empty:
        print(f"No genes classified as {status}.")
        return
    
    # Sort X-axis by Total Mean Expression
    gene_sort = subset.groupby('Gene')['Total_Mean'].first().sort_values().index
    pivot_df = subset.pivot_table(
        index='Gene', 
        columns='Age_Category', 
        values='Relative_Contribution', 
        aggfunc='sum'
    ).reindex(gene_sort)
    
    fig_w = max(16, len(pivot_df) * 0.4)
    fig, ax = plt.subplots(figsize=(fig_w, 10), dpi=300)
    
    pivot_df.plot(kind='bar', stacked=True, ax=ax, color=colors, width=0.8, edgecolor='white', linewidth=0.3)
    
    ax.set_title(f'100% Stacked Bar Chart: {status} Genes (Dev <= 6% & p > 0.05)', fontsize=22, fontweight='bold', pad=30)
    ax.set_ylabel('Proportional Contribution', fontsize=14, fontweight='bold')
    ax.set_xlabel('Gene Symbol (Sorted by Total Mean Expression)', fontsize=14, fontweight='bold')
    ax.set_ylim(0, 1.0)
    
    ax.yaxis.set_major_locator(plt.FixedLocator([0.25, 0.50, 0.75, 1.0]))
    ax.yaxis.grid(True, linestyle='--', alpha=0.3, color='grey')
    sns.despine(ax=ax, top=True, right=True)
    
    plt.xticks(rotation=90, ha='center', fontsize=9)
    ax.legend(reversed(ax.get_legend_handles_labels()[0]), reversed(ax.get_legend_handles_labels()[1]), 
              title="Developmental Stage", bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False)
    
    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight')
    plt.close()

# Run the plot generation
generate_plot('Stable', 'Stable_Genes_Development.png')
generate_plot('Dynamic', 'Dynamic_Genes_Development.png')

print("Analysis complete. Classification saved and plots generated.")
