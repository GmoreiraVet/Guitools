import pandas as pd
import glob
import plotly.graph_objects as go
import os
import numpy as np
from scipy.optimize import curve_fit

def load_bracken_files(input_folder):
    """
    Load Bracken report files from a specified folder and prepare data for rarefaction.

    Args:
    - input_folder (str): Path to the folder containing Bracken report files.

    Returns:
    - pd.DataFrame: Combined DataFrame with sample_id, taxonomy_id, and number of reads (new_est_reads).
    """
    file_paths = glob.glob(os.path.join(input_folder, "*.txt"))
    dataframes = []

    for file in file_paths:
        sample_id = os.path.basename(file).replace("_bracken.txt", "")
        df = pd.read_csv(file, sep="\t", usecols=["taxonomy_id", "new_est_reads"])
        df["sample_id"] = sample_id
        dataframes.append(df)

    combined_df = pd.concat(dataframes)
    return combined_df

def rarefaction_curve(data, sample_id):
    """
    Generate rarefaction curve for a single sample based on reads.

    Args:
    - data (pd.DataFrame): Filtered DataFrame for a specific sample.
    - sample_id (str): The sample ID to analyze.

    Returns:
    - dict: Dictionary of subsampling levels and corresponding unique OTUs (taxa).
    """
    sample_data = data[data["sample_id"] == sample_id]
    rarefaction_data = {}

    total_reads = sample_data["new_est_reads"].sum()
    for depth in range(10, total_reads, 10):
        subsample = sample_data.sample(n=depth, replace=True, weights=sample_data["new_est_reads"])
        unique_taxa = subsample["taxonomy_id"].nunique()
        rarefaction_data[depth] = unique_taxa

    return rarefaction_data

def log_model(x, a, b):
    """ Logarithmic model: y = a * log(x) + b """
    return a * np.log(x) + b

def plot_rarefaction_curves_html(rarefaction_results, output_file="rarefaction_curve_with_log_fit.html"):
    """
    Plot rarefaction curves for all samples as an interactive HTML.

    Args:
    - rarefaction_results (dict): Dictionary of sample_id and their rarefaction data.
    - output_file (str): Path to save the rarefaction curve HTML file.
    """
    fig = go.Figure()
    color_scale = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    color_idx = 0
    sample_names = {}
    max_otus = 0

    for sample_id, data in rarefaction_results.items():
        sample_name = input(f"Enter a name for sample '{sample_id}': ")
        sample_names[sample_id] = sample_name

        depths = list(data.keys())
        species_counts = list(data.values())

        popt, _ = curve_fit(log_model, depths, species_counts, maxfev=10000)
        smooth_depths = np.linspace(min(depths), max(depths), 500)
        smooth_counts = log_model(smooth_depths, *popt)

        color = color_scale[color_idx % len(color_scale)]
        fig.add_trace(go.Scatter(x=smooth_depths, y=smooth_counts, mode='lines', name=sample_name, line=dict(color=color)))

        std_counts = np.std(species_counts)

        fill_color = f"rgba{tuple(int(color[i:i+2], 16) for i in (1, 3, 5)) + (0.2,)}"
        fig.add_trace(go.Scatter(
            x=smooth_depths.tolist() + smooth_depths[::-1].tolist(),
            y=[count + std_counts for count in smooth_counts] + [count - std_counts for count in smooth_counts[::-1]],
            fill='toself', fillcolor=fill_color,
            line=dict(color='rgba(255, 255, 255, 0)'), name=f"SD Area: {sample_name}"
        ))

        max_otus = max(max_otus, max(species_counts))
        color_idx += 1

    plot_title = input("Enter a title for the rarefaction curve plot: ")

    fig.update_layout(
        title=plot_title,
        xaxis_title="Number of Reads Sampled",
        yaxis_title="Unique OTUs (Taxa)",
        legend_title="Samples",
        template="plotly_white",
        yaxis=dict(range=[0, max_otus * 1.2])
    )

    fig.write_html(output_file)
    fig.show()
    print(f"Rarefaction curve saved to {output_file}")

def print_success_ascii():
    """ Prints a bunny with a speech bubble saying 'SUCCESS'. """
    print("""
    
  SUCCESS !! 
   
   __     __
  /_/|   |\\_\\  
   |U|___|U|
   |       |
   | ,   , |
  (  = Y =  )
   |      |
  /|       |\\
  \\| |   | |/
 (_|_|___|_|_)
   '"'   '"'

------------------------------------------------

    """)

if __name__ == "__main__":
    input_folder = input("Enter the path to the folder containing Bracken files: ").strip()
    
    if not os.path.exists(input_folder):
        print("Error: The specified folder does not exist. Please check the path and try again.")
    else:
        output_file = "rarefaction_curve_with_log_fit.html"

        data = load_bracken_files(input_folder)
        rarefaction_results = {}

        for sample_id in data["sample_id"].unique():
            rarefaction_results[sample_id] = rarefaction_curve(data, sample_id)

        plot_rarefaction_curves_html(rarefaction_results, output_file)

        print_success_ascii()
