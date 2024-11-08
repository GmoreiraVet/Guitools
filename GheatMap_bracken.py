import pandas as pd
import glob
import plotly.express as px
import os

def load_bracken_files(input_folder, rank="G", top_n=20): #Change parameters here and at the bottom of the script, for reasons
    """
    Load and combine Bracken report files from a specified folder, filter by taxonomic level (rank), 
    and keep top N taxa.
    
    Args:
    - input_folder (str): Path to the folder containing Bracken report files.
    - rank (str): The taxonomic level to filter (e.g., 'G' for genus, 'S' for species).
    - top_n (int): Number of top taxa to retain based on overall abundance.
    
    Returns:
    - pd.DataFrame: A pivoted DataFrame ready for heatmap visualization.
    """
    # Get a list of all Bracken files in the input folder
    file_paths = glob.glob(os.path.join(input_folder, "*.txt"))
    dataframes = []

    # Load each Bracken report file
    for file in file_paths:
        sample_id = os.path.basename(file).replace("_bracken.txt", "")  # Extract sample ID from file name
        df = pd.read_csv(file, sep="\t", usecols=["name", "taxonomy_lvl", "fraction_total_reads"])
        df["sample_id"] = sample_id  # Add sample ID column
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

def plot_heatmap(data, output_file="taxonomic_abundance_heatmap.html"):
    """
    Plot a heatmap of taxonomic abundance using plotly and save to file.
    
    Args:
    - data (pd.DataFrame): The pivoted data ready for heatmap plotting.
    - output_file (str): Path to save the generated heatmap image.
    """
    if data.empty:
        print("No data available to plot the heatmap.")
        return
    
    fig = px.imshow(data, 
                    labels=dict(x="Sample ID", y="Taxon", color="Relative Abundance"),
                    title="Taxonomic Abundance Heatmap Across Samples",
                    color_continuous_scale='viridis')
    fig.update_layout(xaxis_tickangle=-45)
    fig.write_html(output_file)
    fig.show()
    print(f"Heatmap saved to {output_file}")

if __name__ == "__main__":
    # Define input folder and output file
    input_folder = "/home/viroicbas/scriptTeste/bracken_reports"  # Change this to your folder with Bracken files. Names must end with ...bracken_reports
    output_file = "taxonomic_abundance_heatmap.html"  # Change this to your desired output file path

    # Load data and generate heatmap
    heatmap_data = load_bracken_files(input_folder, rank="G", top_n=20)  # Change rank to "S" for species-level, change top in accordance with the top of the script
    plot_heatmap(heatmap_data, output_file)

