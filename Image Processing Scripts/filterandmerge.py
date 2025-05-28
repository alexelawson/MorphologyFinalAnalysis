import pandas as pd

# Example: load or define your three dataframes
df1 = pd.read_csv('/Users/alexlawson/Masters-Data-Final/Morphoglia/PVN_all/Morphology.csv')
df2 = pd.read_csv('/Users/alexlawson/Desktop/mr196.csv')

# Step 2: Check column headers match
if list(df1.columns) == list(df2.columns):
    # Step 3: Merge the DataFrames
    merged_df = pd.concat([df1, df2], ignore_index=True)
    print("DataFrames merged successfully.")

    # Step 4: Save to CSV
    merged_df.to_csv("/Users/alexlawson/Masters-Data-Final/Morphoglia/PVN_all/merged_filtered_data.csv", index=False)
    print("Saved merged DataFrame to 'merged_filtered_data.csv'.")
else:
    print("Error: DataFrames do not have the same column headers.")
    print("df1 columns:", df1.columns)
    print("df2 columns:", df2.columns)
