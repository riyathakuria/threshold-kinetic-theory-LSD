import pandas as pd

# Define the local file path
file_path = r'/Users/riyathakuria/Documents/LSD_Insilico_Analysis_24:3:26/Raw data/PPI_SpatioTemporal_Expression.csv'

# Load the dataset
df = pd.read_csv(file_path)

# Function to map Age strings to developmental stages
def map_age_to_stage(age_str):
    if 'pcw' in age_str:
        return 'Prenatal'
    elif 'mos' in age_str:
        return 'Postnatal'
    elif 'yrs' in age_str:
        # Extract the numeric year
        years = int(age_str.split(' ')[0])
        if years <= 4:
            return 'Early Childhood'
        else:
            return 'Late Childhood'
    return 'Unknown'

# 1. Apply the stage mapping
df['Stage'] = df['Age'].apply(map_age_to_stage)

# 2. Group by Gene and Stage only (removing Lineage) and calculate average expression
summary_df = df.groupby(['Gene', 'Stage'])['Expression'].mean().reset_index()

# 3. Create the summarized table (Pivot)
pivot_df = summary_df.pivot_table(
    index='Gene', 
    columns='Stage', 
    values='Expression'
).reset_index()

# 4. Reorder columns for chronological progression
ordered_stages = ['Prenatal', 'Postnatal', 'Early Childhood', 'Late Childhood']
existing_stages = [stage for stage in ordered_stages if stage in pivot_df.columns]
final_columns = ['Gene'] + existing_stages
pivot_df = pivot_df[final_columns]

# Save the summarized result to a new CSV file
pivot_df.to_csv('summarized_expression_by_stage_no_lineage.csv', index=False)

print("Processing complete. Summarized table (no lineages) saved.")
print(pivot_df.head())