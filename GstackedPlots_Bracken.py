import pandas as pd
import glob
import plotly.express as px
import os
import re

def natural_sort(file_paths):
    """
    Sort the files using a natural sorting algorithm by extracting the leading number from each file name.
    """
    def extract_number(filename):
        # Extract the number at the beginning of the file name using a regular expression.
        match = re.match(r"(\d+)_", filename)
        return int(match.group(1)) if match else float('inf')  # Return infinity if no number is found
    
    # Sort the files based on the extracted number from the file name
    return sorted(file_paths, key=lambda x: extract_number(os.path.basename(x)))

def load_bracken_files(input_folder, rank="G", top_n=15):
    """
    Load and combine Bracken report files from a specified folder, filter by taxonomic level (rank), 
    and keep top N taxa.
    
    Args:
    - input_folder (str): Path to the folder containing Bracken report files.
    - rank (str): The taxonomic level to filter (e.g., 'G' for genus, 'S' for species).
    - top_n (int): Number of top taxa to retain based on overall abundance.
    
    Returns:
    - pd.DataFrame: A DataFrame ready for stacked bar plot visualization.
    """
    file_paths = glob.glob(os.path.join(input_folder, "*.txt"))
    file_paths = natural_sort(file_paths)  # Sort files naturally by numeric order
    
    # Assign sample names automatically based on filenames
    sample_names = {os.path.basename(file).replace("_bracken.txt", ""): os.path.basename(file).replace("_bracken.txt", "") for file in file_paths}

    dataframes = []

    # Load each Bracken report file
    for file in file_paths:
        sample_id = os.path.basename(file).replace("_bracken.txt", "")
        
        # Use the automatically assigned sample name
        custom_sample_name = sample_names.get(sample_id, sample_id)
        
        # Load the Bracken file and assign the custom sample name
        df = pd.read_csv(file, sep="\t", usecols=["name", "taxonomy_lvl", "fraction_total_reads"])
        df["sample_id"] = custom_sample_name
        dataframes.append(df)

    # Combine all data into a single DataFrame
    combined_df = pd.concat(dataframes)

    # Filter for the specified taxonomic level (rank)
    filtered_df = combined_df[combined_df["taxonomy_lvl"] == rank]

    if filtered_df.empty:
        print(f"No data found for rank '{rank}'")
        return pd.DataFrame()

    # Find the top N taxa across all samples based on fraction_total_reads
    top_taxa = filtered_df.groupby("name")["fraction_total_reads"].sum().nlargest(top_n).index
    top_taxa_df = filtered_df[filtered_df["name"].isin(top_taxa)]
    
    # Calculate the "Other" category (sum of the remaining taxons)
    other_taxa_df = filtered_df[~filtered_df["name"].isin(top_taxa)]
    other_taxa_sum = other_taxa_df.groupby(["sample_id"])["fraction_total_reads"].sum().reset_index()
    other_taxa_sum["name"] = "Other"
    
    # Combine the top taxa with the "Other" category
    combined_df_with_other = pd.concat([top_taxa_df, other_taxa_sum])

    return combined_df_with_other

def plot_stacked_bar(data, output_file="taxonomic_abundance_stacked_bar.html"):
    """
    Plot a stacked bar chart of taxonomic abundance using plotly and save to file.
    
    Args:
    - data (pd.DataFrame): The data ready for stacked bar plotting.
    - output_file (str): Path to save the generated plot in HTML format.
    """
    if data.empty:
        print("No data available to plot.")
        return
    
    # Define a custom color sequence or use a predefined one
    color_palette = px.colors.qualitative.Pastel  # You can choose 'Vivid', 'Pastel', 'Dark2', etc.

    # Get unique taxa names
    taxa_names = data['name'].unique()
    
    # Dynamically assign colors by cycling through the chosen palette
    color_map = {taxa: color_palette[i % len(color_palette)] for i, taxa in enumerate(taxa_names)}
    
    # Force "Other" category to be grey
    color_map["Other"] = "gray"

    # Sort the x-axis labels alphabetically or numerically
    sorted_samples = sorted(data['sample_id'].unique(), key=lambda x: (int(re.match(r'(\d+)', x).group(1)), x))

    # Plot with the automatically mapped colors
    fig = px.bar(data, 
                 x="sample_id", 
                 y="fraction_total_reads", 
                 color="name",
                 color_discrete_map=color_map,
                 labels={"fraction_total_reads": "Relative Abundance", "sample_id": "Sample ID", "name": "Taxon"},
                 title="Taxonomic Abundance Stacked Bar Plot Across Samples")
    
    # Ensure the x-axis is ordered correctly
    fig.update_layout(barmode='stack', xaxis_tickangle=-45, xaxis=dict(categoryorder="array", categoryarray=sorted_samples))
    
    fig.write_html(output_file)
    fig.show()
    print(f"Stacked bar plot saved to {output_file}")

if __name__ == "__main__":
    input_folder = "/home/viroicbas/scriptTeste/bracken_reports"
    output_file = "taxonomic_abundance_stacked_bar.html"
    
    # Load data and generate stacked bar plot
    bar_data = load_bracken_files(input_folder, rank="G", top_n=15)
    plot_stacked_bar(bar_data, output_file)
