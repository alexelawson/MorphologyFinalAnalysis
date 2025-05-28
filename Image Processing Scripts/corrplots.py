import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency
from matplotlib.colors import LinearSegmentedColormap

# ---------------------------
# SETUP PATHS AND LOAD DATA
# ---------------------------
# Define paths (adjust as needed)
data_path = "/Users/alexlawson/Masters-Data-Final/Morphoglia/PVN_all/Morphoglia/MorphoGlia/Data/Morphology_HDBSCAN_30_0.1_150_5_reformatted.csv"
plot_path = "/Users/alexlawson/Masters-Data-Final/Morphoglia/PVN_all/Morphoglia/MorphoGlia/Data/stats"

# Create the output directory if it doesn't exist
if not os.path.exists(plot_path):
    os.makedirs(plot_path)

# Load the CSV data
data = pd.read_csv(data_path)

# Filter out noise: remove rows where Clusters == -1
data = data[data['Clusters'] != -1]

# ---------------------------
# SPLIT DATA BY TISSUE & PREP LABELS
# ---------------------------
# Subset by tissue
data_M = data[data['tissue'] == "M"].copy()
data_F = data[data['tissue'] == "F"].copy()

# Create "Cluster_Labels" as a categorical column (if needed)
data_M['Cluster_Labels'] = data_M['Clusters'].astype('category')
data_F['Cluster_Labels'] = data_F['Clusters'].astype('category')

# Create "categories" by replacing underscores in treatment names with spaces
data_M['categories'] = data_M['treatment'].str.replace('_', ' ')
data_F['categories'] = data_F['treatment'].str.replace('_', ' ')

# ---------------------------
# CREATE CONTINGENCY TABLES & CHI-SQUARE TESTS
# ---------------------------
# Create contingency tables: rows = treatment, columns = Clusters
contingency_table_M = pd.crosstab(data_M['treatment'], data_M['Clusters'])
contingency_table_F = pd.crosstab(data_F['treatment'], data_F['Clusters'])

print("Contingency Table - Males:")
print(contingency_table_M)
print("\nContingency Table - Females:")
print(contingency_table_F)

# Chi-square test for Males
chi2_M, p_M, dof_M, expected_M = chi2_contingency(contingency_table_M)
# Compute standardized residuals for males: (observed - expected) / sqrt(expected)
residuals_M = (contingency_table_M.values - expected_M) / np.sqrt(expected_M)
residuals_M_df = pd.DataFrame(residuals_M, index=contingency_table_M.index, columns=contingency_table_M.columns)
print("\nChi-square residuals for Males:")
print(residuals_M_df)

# Chi-square test for Females
chi2_F, p_F, dof_F, expected_F = chi2_contingency(contingency_table_F)
# Compute standardized residuals for females
residuals_F = (contingency_table_F.values - expected_F) / np.sqrt(expected_F)
residuals_F_df = pd.DataFrame(residuals_F, index=contingency_table_F.index, columns=contingency_table_F.columns)
print("\nChi-square residuals for Females:")
print(residuals_F_df)

# ---------------------------
# PLOT: BALLOON PLOT FOR MALES
# ---------------------------
plt.figure(figsize=(8, 6))
# Get x (Clusters) and y (treatment) labels from the contingency table
x_labels = contingency_table_M.columns.astype(str)
y_labels = contingency_table_M.index.astype(str)

# Create a meshgrid for plotting
x, y = np.meshgrid(np.arange(len(x_labels)), np.arange(len(y_labels)))

# Flatten arrays for plotting scatter
x_flat = x.flatten()
y_flat = y.flatten()
counts = contingency_table_M.values.flatten()

# Set a scaling factor for bubble sizes (adjust for aesthetics)
scale = 200
sizes = counts * scale

# Plot the bubble plot: each bubble size reflects the count
plt.scatter(x_flat, y_flat, s=sizes, color='lightblue', alpha=0.6, edgecolors='black')

# Annotate each bubble with its count
for (i, j, count) in zip(x_flat, y_flat, counts):
    plt.text(i, j, str(count), ha='center', va='center', fontsize=12)

# Set axes ticks and labels
plt.xticks(np.arange(len(x_labels)), x_labels)
plt.yticks(np.arange(len(y_labels)), y_labels)
plt.xlabel("Clusters")
plt.ylabel("Treatment")
plt.title("Cluster Frequencies - Males")
plt.gca().invert_yaxis()  # Invert y-axis to mirror typical balloonplot orientation
plt.tight_layout()

# Save the balloon plot
balloon_plot_file = os.path.join(plot_path, "Cluster_Frequencies-M-reformatted-python.png")
plt.savefig(balloon_plot_file, dpi=300)
plt.close()

# ---------------------------
# PLOT: CORRELATION PLOT (HEATMAP) FOR FEMALES
# ---------------------------
# Define universal scale for the heatmap
vmin, vmax = -1.67, 1.67

# Create a custom color palette from your specified colors
colors = [
    "#2166AC", "#92C5DE",  # negative residuals (pale blue)
    "#F0F0F0", "#F7F7F7", "#FFFFFF", "#F7F7F7",  # widened white zone
    "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"  # positive residuals (strong red)
]
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=100)

plt.figure(figsize=(10, 9))
# Plot the heatmap with a fixed color range to highlight the low variation in females
sns.heatmap(residuals_F_df, cmap=custom_cmap, annot=True, fmt=".2f",
            vmin=vmin, vmax=vmax, linewidths=0.5, linecolor='black',
            cbar_kws={'shrink': 0.8})

plt.title("Chi-Square Residuals - Females\n(Universal Scale: -1.67 to 1.67)")
plt.tight_layout()

# Save the heatmap
heatmap_file = os.path.join(plot_path, "Chi_Square_Correlogram-F-reformatted-python.png")
plt.savefig(heatmap_file, dpi=300)
plt.close()

print(f"Plots saved:\n- Balloon Plot: {balloon_plot_file}\n- Heatmap: {heatmap_file}")
