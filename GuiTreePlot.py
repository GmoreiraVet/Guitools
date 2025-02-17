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
    file_paths = glob.glob(os.path.join(input_folder, "*.txt"))
    file_paths = natural_sort(file_paths)
    sample_ids = [os.path.basename(file).replace("_bracken.txt", "") for file in file_paths]
    name_mapping = get_custom_sample_names(sample_ids)
    dataframes = []
    for file in file_paths:
        sample_id = os.path.basename(file).replace("_bracken.txt", "")
        custom_sample_name = name_mapping[sample_id]
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

def plot_treemap(data, output_file="taxonomic_abundance_treemap.html"):
    if data.empty:
        print("No data available to plot.")
        return
    fig = px.treemap(data,
                     path=["sample_id", "name"],
                     values="fraction_total_reads",
                     color="name",
                     title="Taxonomic Abundance Treemap Across Samples",
                     color_discrete_sequence=px.colors.qualitative.Pastel)
    fig.write_html(output_file)
    fig.show()
    print(f"Treemap plot saved to {output_file}")

if __name__ == "__main__":
    input_folder = "/home/viroicbas/scriptTeste/bracken_reports"
    output_file = "taxonomic_abundance_treemap.html"
    treemap_data = load_bracken_files(input_folder, rank="G", top_n=15)
    plot_treemap(treemap_data, output_file)
