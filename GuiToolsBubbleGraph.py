import pandas as pd
import glob
import plotly.express as px
import os
import re

def natural_sort(file_paths):
    def extract_number(filename):
        match = re.match(r"(\d+)_", filename)
        return int(match.group(1)) if match else float('inf')
    return sorted(file_paths, key=lambda x: extract_number(os.path.basename(x)))

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
    Load and process Bracken files with an option to rename samples and combine results.
    """
    file_paths = glob.glob(os.path.join(input_folder, "*.txt"))
    file_paths = natural_sort(file_paths)
    sample_ids = [os.path.basename(file).replace("_bracken.txt", "") for file in file_paths]
    
    # Allow custom renaming of samples
    sample_names = get_custom_sample_names(sample_ids)

    dataframes = []
    for file in file_paths:
        sample_id = os.path.basename(file).replace("_bracken.txt", "")
        custom_sample_name = sample_names.get(sample_id, sample_id)
        df = pd.read_csv(file, sep="\t", usecols=["name", "taxonomy_lvl", "fraction_total_reads"])
        df["sample_id"] = custom_sample_name
        dataframes.append(df)

    combined_df = pd.concat(dataframes)
    filtered_df = combined_df[combined_df["taxonomy_lvl"] == rank]

    if filtered_df.empty:
        print(f"No data found for rank '{rank}'")
        return pd.DataFrame()

    top_taxa = filtered_df.groupby("name")["fraction_total_reads"].sum().nlargest(top_n).index
    top_taxa_df = filtered_df[filtered_df["name"].isin(top_taxa)]
    other_taxa_df = filtered_df[~filtered_df["name"].isin(top_taxa)]
    other_taxa_sum = other_taxa_df.groupby(["sample_id"])["fraction_total_reads"].sum().reset_index()
    other_taxa_sum["name"] = "Other"

    combined_df_with_other = pd.concat([top_taxa_df, other_taxa_sum])

    return combined_df_with_other

def plot_bubble_chart(data, output_file="taxonomic_abundance_bubble_chart.html"):
    """
    Generate an interactive bubble chart using Plotly.
    """
    if data.empty:
        print("No data available to plot.")
        return

    fig = px.scatter(data,
                     x="sample_id",
                     y="name",
                     size="fraction_total_reads",
                     color="name",
                     labels={"fraction_total_reads": "Relative Abundance", "sample_id": "Sample ID", "name": "Taxon"},
                     title="Taxonomic Abundance Bubble Chart Across Samples",
                     size_max=60)

    fig.update_layout(xaxis_tickangle=-45)
    fig.write_html(output_file)
    fig.show()
    print(f"Bubble chart saved to {output_file}")

if __name__ == "__main__":
    input_folder = "/home/viroicbas/scriptTeste/bracken_reports"
    output_file = "taxonomic_abundance_bubble_chart.html"

    bubble_data = load_bracken_files(input_folder, rank="G", top_n=15)
    plot_bubble_chart(bubble_data, output_file)
