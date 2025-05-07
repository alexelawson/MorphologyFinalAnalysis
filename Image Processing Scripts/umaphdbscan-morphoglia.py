import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

def rgb_string_to_tuple(rgb_string):
    """Convert '(R,G,B)' string to reversed normalized RGB tuple (0–1 scale)."""
    rgb = tuple(int(x.strip()) for x in rgb_string.strip('()').split(','))
    return tuple(x / 255 for x in reversed(rgb))

def plot_umap_with_clusters(
    input_path="your/hardcoded/input/umap_data.csv",
    output_path="your/hardcoded/output/umap_plot.png"
):
    # Read CSV
    df = pd.read_csv(input_path)
    df = df[df['Clusters'] != -1]

    # Convert RGB string to matplotlib-compatible color
    df['Color'] = df['Cluster_Color'].apply(rgb_string_to_tuple)

    # Create plot
    fig, ax = plt.subplots(figsize=(9, 8))

    for cluster_id in sorted(df['Clusters'].unique()):
        cluster_data = df[df['Clusters'] == cluster_id]
        ax.scatter(
            cluster_data['UMAP_1'],
            cluster_data['UMAP_2'],
            label=f"Cluster {cluster_id}",
            c=[cluster_data['Color'].iloc[0]],
            s=30,
            edgecolors='grey',   # ← Outline color
            linewidths=0.25,
        )

    # Create custom legend
    handles = [
        mpatches.Patch(
            facecolor=rgb_string_to_tuple(df[df['Clusters'] == cluster]['Cluster_Color'].iloc[0]),
            edgecolor='grey',
            linewidth=0.25,
            label=f'Cluster {cluster}'
        )
        for cluster in sorted(df['Clusters'].unique())
    ]
    ax.legend(
        handles=handles, 
        loc='upper right', 
        title='Clusters',
        frameon=True)

    ax.set_title("UMAP Colored by HDBSCAN Clusters")
    ax.set_xlabel("UMAP_1")
    ax.set_ylabel("UMAP_2")
    plt.tight_layout()

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Save figure
    plt.savefig(output_path, dpi=300)
    print(f"UMAP plot saved to: {output_path}")

    plt.close()

# Example usage
plot_umap_with_clusters(
    input_path="/Users/alexlawson/Masters-Data-Final/Morphoglia/PVN_all/Morphoglia/MorphoGlia/Data/Morphology_HDBSCAN_30_0.1_150_5-v2.csv",
    output_path="/Users/alexlawson/Masters-Data-Final/Morphoglia/PVN_all/Morphoglia/MorphoGlia/Data/Plots/umap_hdbscan_clusters-reformatted02.png"
)
