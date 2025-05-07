import pandas as pd

# Load the CSV file
file_path = "/Users/alexlawson/Masters-Data-Final/Morphoglia/PVN_all/Morphoglia/MorphoGlia/Data/Morphology_HDBSCAN_30_0.1_150_5.csv"
df = pd.read_csv(file_path)

# Define the remapping
cluster_mapping = {
    0: 0,
    4: 1,
    3: 2,
    1: 3,
    2: 4
}

# Apply the mapping to the 'Clusters' column
df['Clusters'] = df['Clusters'].map(cluster_mapping)

# Save to a new CSV file
output_path = "/Users/alexlawson/Masters-Data-Final/Morphoglia/PVN_all/Morphoglia/MorphoGlia/Data/Morphology_HDBSCAN_30_0.1_150_5_reformatted.csv"
df.to_csv(output_path, index=False)

print(f"Reformatted clusters saved to {output_path}")
