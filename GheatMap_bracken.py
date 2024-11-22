import pandas as pd
import glob
import plotly.express as px
import os

def get_custom_sample_names(original_names):
    """
    Prompt the user to enter custom names for samples.
    
    Args:
    - original_names (list): List of original sample names.
    
    Returns:
    - dict: A dictionary mapping original names to custom names.
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
                    labels=dict(x="", y="", color="Relative Abundance"),
                    title="Taxonomic Abundance Heatmap Across Samples",
                    color_continuous_scale='mint')
    fig.update_layout(xaxis_tickangle=-45)

    # Save as HTML
    fig.write_html(output_file)
    print(f"Heatmap saved to {output_file}")
    
    # Show the heatmap
    fig.show()

    # Ask if the user wants to save as PNG
    save_as_png = input("Would you like to save the heatmap as a PNG file? (yes/no): ").strip().lower()
    if save_as_png in ["yes", "y"]:
        png_file = output_file.replace(".html", ".png")
        fig.write_image(png_file)
        print(f"Heatmap saved as PNG to {png_file}")

if __name__ == "__main__":
    # Define input folder and output file
    input_folder = "/home/viroicbas/scriptTeste/bracken_reports"  # Change this to your folder with Bracken files
    output_file = "taxonomic_abundance_heatmap.html"  # Change this to your desired output file path

    # Load data and generate heatmap
    heatmap_data = load_bracken_files(input_folder, rank="G", top_n=15)  # Change rank to "S" for species-level if needed
    plot_heatmap(heatmap_data, output_file)
