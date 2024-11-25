import pandas as pd
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage
from matplotlib.colors import LinearSegmentedColormap

def get_custom_sample_names(original_names):
    """
    Prompt the user to enter custom names for samples.
    """
    print("Original sample names detected:")
    for name in original_names:
        print(f" - {name}")
    
    print("\nEnter custom names for each sample (press Enter to keep the original name):")
    name_mapping = {}
    for name in original_names:
        custom_name = input(f"Custom name for '{name}': ").strip()
        name_mapping[name] = custom_name if custom_name else name
    
    print("\nSample names have been updated:")
    for old, new in name_mapping.items():
        print(f" - {old} → {new}")
    
    return name_mapping

def load_bracken_files(input_folder, rank="G", top_n=15):
    """
    Load and combine Bracken report files from a specified folder, filter by taxonomic level (rank), 
    and keep top N taxa.
    """
    # Get a list of all Bracken files in the input folder
    file_paths = glob.glob(os.path.join(input_folder, "*.txt"))
    dataframes = []

    # Extract original sample names
    sample_ids = [os.path.basename(file).replace("_bracken.txt", "") for file in file_paths]
    
    # Allow user to rename samples
    name_mapping = get_custom_sample_names(sample_ids)

    # Load each Bracken report file
    for file in file_paths:
        sample_id = os.path.basename(file).replace("_bracken.txt", "")  # Extract sample ID from file name
        custom_sample_id = name_mapping[sample_id]  # Map to custom name
        df = pd.read_csv(file, sep="\t", usecols=["name", "taxonomy_lvl", "fraction_total_reads"])
        df["sample_id"] = custom_sample_id  # Add sample ID column with custom name
        dataframes.append(df)

    # Combine all data into a single DataFrame
    combined_df = pd.concat(dataframes)

    # Ensure 'taxonomy_lvl' is being read as expected
    print(f"Unique values in 'taxonomy_lvl' column: {combined_df['taxonomy_lvl'].unique()}")

    # Filter for the specified taxonomic level (rank)
    filtered_df = combined_df[combined_df["taxonomy_lvl"] == rank]

    # Check if filtering worked by displaying a few records
    print(f"Filtered data for rank '{rank}':")
    print(filtered_df.head())

    # If there's no data after filtering, return an empty dataframe
    if filtered_df.empty:
        print(f"No data found for rank '{rank}'")
        return pd.DataFrame()

    # Find the top N taxa across all samples based on fraction_total_reads
    top_taxa = filtered_df.groupby("name")["fraction_total_reads"].sum().nlargest(top_n).index
    top_taxa_df = filtered_df[filtered_df["name"].isin(top_taxa)]

    # Pivot data for the heatmap
    heatmap_data = top_taxa_df.pivot(index="name", columns="sample_id", values="fraction_total_reads")
    heatmap_data.fillna(0, inplace=True)  # Replace NaNs with 0 for absent taxa

    return heatmap_data

def plot_clustergram(data, output_file="taxonomic_abundance_clustergram.png"):
    """
    Plot a clustergram of taxonomic abundance using Seaborn and save to file.
    """
    if data.empty:
        print("No data available to plot the clustergram.")
        return
    
    # Compute Bray-Curtis distances
    bray_curtis_dist = pdist(data.T, metric="braycurtis")
    linkage_matrix = linkage(bray_curtis_dist, method="average")

    # Set Seaborn style with faint blue background
    sns.set(style="whitegrid", palette="muted", rc={"axes.facecolor": "#E0F7FA"})  # Light background

    # Custom color scale based on your provided colors (Light grey, Ice, Pale teal, Dusty teal, Dark aqua, Deep teal)
    colors = ['#E3F0E0', '#AED6C9', '#7CB7AE', '#539895', '#30797A', '#0E5960']
    custom_cmap = LinearSegmentedColormap.from_list("custom_mint", colors)

    # Adjust figure size based on the data dimensions (making cells square)
    num_rows, num_cols = data.shape
    figsize = (num_cols, num_rows)  # Ensure 1:1 aspect ratio for square cells

    # Create the clustergram with square cells
    clustergrid = sns.clustermap(
        data,
        row_cluster=True,
        col_cluster=True,
        method="average",
        metric="braycurtis",
        cmap=custom_cmap,  # Applying the custom continuous color map
        linewidths=0.5,
        figsize=figsize,  # Adjusted size for square cells
        square=True  # Ensure square cells
    )

    # Save the clustergram as a PNG
    plt.title("Taxonomic Abundance Clustergram")
    plt.savefig(output_file, dpi=300)  # Saving with high resolution
    print(f"Clustergram saved to {output_file}")
    plt.show()

if __name__ == "__main__":
    # Define input folder and output file
    input_folder = "/home/viroicbas/scriptTeste/bracken_reports"  # Change this to your folder with Bracken files
    output_file = "taxonomic_abundance_clustergram.png"  # Change this to your desired output file path

    # Load data and generate clustergram
    heatmap_data = load_bracken_files(input_folder, rank="G", top_n=15)  # Change rank to "S" for species-level if needed
    plot_clustergram(heatmap_data, output_file)
