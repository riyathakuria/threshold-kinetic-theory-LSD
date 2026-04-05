import pandas as pd

# Define File Paths
expression_path = r"/Users/riyathakuria/Documents/LSD_Insilico_Analysis_24:3:26/Raw data/PPI_SpatioTemporal_Expression.csv"
stability_path = r"/Users/riyathakuria/Documents/LSD_Insilico_Analysis_24:3:26/Raw data/Gene_Stability_Classification.csv"

# 1. Load Datasets
df_expr = pd.read_csv(expression_path)
df_stab = pd.read_csv(stability_path)

# 2. Chronological Stage Mapping
# This ensures that temporal bins are consistent across both files
def map_developmental_stage(age_str):
    if 'pcw' in age_str: return 'Prenatal'
    if 'mos' in age_str: return 'Postnatal'
    if 'yrs' in age_str:
        years = int(age_str.split(' ')[0])
        return 'Early Childhood' if years <= 4 else 'Late Childhood'
    return 'Unknown'

df_expr['Stage'] = df_expr['Age'].apply(map_developmental_stage)

# 3. Mathematical Normalization (Aggregation)
# Calculates the Mean Expression for every Gene in each Brain Region (Lineage) per Stage
# This provides the stable baseline required for further anchor and driver computations
agg_expr = df_expr.groupby(['Gene', 'Lineage', 'Stage'])['Expression'].mean().reset_index()

# 4. Pivot Stages into Columns
# Reformatting the data to have one row per Gene-Region pair for downstream processing
pivot_expr = agg_expr.pivot_table(
    index=['Gene', 'Lineage'], 
    columns='Stage', 
    values='Expression'
).reset_index()

# Reorder columns chronologically
ordered_stages = ['Prenatal', 'Postnatal', 'Early Childhood', 'Late Childhood']
pivot_expr = pivot_expr[['Gene', 'Lineage'] + ordered_stages]

# 5. Relational Merge with Global Stability Metrics
# This connects the regional profile to the global Avg_Dev, p_value, and Stability_Status
master_df = pivot_expr.merge(df_stab, on='Gene', how='inner')

# 6. Save the Master CSV
master_df.to_csv('Master_Regional_Stability_Expression.csv', index=False)

print("Master File successfully generated with regional mapping.")
print(master_df.head())