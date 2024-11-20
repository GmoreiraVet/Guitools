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

def load_bracken_files(input_folder, rank="G", top_n=15):
    file_paths = glob.glob(os.path.join(input_folder, "*.txt"))
    file_paths = natural_sort(file_paths)
    sample_names = {os.path.basename(file).replace("_bracken.txt", ""): os.path.basename(file).replace("_bracken.txt", "") for file in file_paths}

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

