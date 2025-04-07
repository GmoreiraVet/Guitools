import pandas as pd
import glob
import plotly.graph_objects as go
import os
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import entropy

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
        # Load Bracken data, focusing on taxonomy_id and new_est_reads columns
        df = pd.read_csv(file, sep="\t", usecols=["taxonomy_id", "new_est_reads"])
        df["sample_id"] = sample_id
        dataframes.append(df)

    combined_df = pd.concat(dataframes)
    return combined_df

def shannon_index(taxa_counts):
    """
    Calculate the Shannon Index (Shannon-Wiener index) for a given list of taxa counts.

    Args:
    - taxa_counts (list): List of counts for each unique taxa in the subsample.

    Returns:
    - float: Shannon Index value.
    """
    total_count = sum(taxa_counts)
    proportions = [count / total_count for count in taxa_counts]
    return entropy(proportions, base=np.e)

def rarefaction_curve_with_shannon(data, sample_id):
    """
    Generate rarefaction curve for a single sample, including Shannon Index at each depth.

    Args:
    - data (pd.DataFrame): Filtered DataFrame for a specific sample.
    - sample_id (str): The sample ID to analyze.

    Returns:
    - dict: Dictionary of subsampling levels and corresponding Shannon Index values.
    """
    sample_data = data[data["sample_id"] == sample_id]
    rarefaction_data = []

    # Simulate subsampling at increasing depths starting from 1
    total_reads = sample_data["new_est_reads"].sum()
    for depth in range(1, total_reads + 1):  # Start at 1 read
        # Subsample by sampling "depth" reads (randomly selecting from the available reads)
        subsample = sample_data.sample(n=depth, replace=True, weights=sample_data["new_est_reads"])

        # Count the occurrences of each taxa in the subsample
        taxa_counts = subsample["taxonomy_id"].value_counts()

        # Calculate the Shannon Index for this subsample
        shannon_val = shannon_index(taxa_counts)

        # Store Shannon Index and depth for this subsample
        rarefaction_data.append((depth, shannon_val))

    return np.array(rarefaction_data)

def log_model(x, a, b):
    """
    Logarithmic model: y = a * log(x) + b
    """
    return a * np.log(x) + b

def plot_rarefaction_curves_html_with_shannon(rarefaction_results, output_file="rarefaction_curve_with_shannon.html"):
    fig = go.Figure()

    color_scale = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    color_idx = 0

    max_shannon = 0  # Track the maximum Shannon Index

    for sample_id, data in rarefaction_results.items():
        sample_name = input(f"Enter a name for sample '{sample_id}': ")

        depths = data[:, 0]  # Subsampling depths
        shannon_values = data[:, 1]  # Shannon Index values

        # Fit the rarefaction data to a logarithmic model
        popt, _ = curve_fit(log_model, depths, shannon_values, maxfev=10000)

        # Generate smoothed curve using the fitted logarithmic model
        smooth_depths = np.linspace(min(depths), max(depths), 500)
        smooth_shannon_values = log_model(smooth_depths, *popt)

        # Calculate standard deviation of Shannon values for the shaded area
        # Here, we approximate it by computing the standard deviation of the observed Shannon values
        std_dev = np.std(shannon_values)
        upper_bound = smooth_shannon_values + std_dev
        lower_bound = smooth_shannon_values - std_dev

        color = color_scale[color_idx % len(color_scale)]

        # Plot the smoothed curve
        fig.add_trace(go.Scatter(x=smooth_depths, y=smooth_shannon_values, mode='lines', name=sample_name,
                                 line=dict(color=color)))

        # Plot the shaded area representing the standard deviation
        fig.add_trace(go.Scatter(
            x=np.concatenate([smooth_depths, smooth_depths[::-1]]),
            y=np.concatenate([upper_bound, lower_bound[::-1]]),
            fill='toself',
            fillcolor=f"rgba{tuple(int(color[i:i+2], 16) for i in (1, 3, 5)) + (0.2,)}",  # Apply RGBA color
            line=dict(color='rgba(255, 255, 255, 0)'),  # Hide border line
            name=f"SD Area: {sample_name}"
        ))

        max_shannon = max(max_shannon, max(smooth_shannon_values))
        color_idx += 1

    plot_title = input("Enter a title for the rarefaction curve plot: ")

    fig.update_layout(
        title=plot_title,
        xaxis_title="Number of Reads Sampled",
        yaxis_title="Shannon Index",
        legend_title="Samples",
        template="plotly_white",
        yaxis=dict(range=[0, max_shannon * 1.2])
    )

    fig.write_html(output_file)
    fig.show()
    print(f"Rarefaction curve with Shannon Index and SD area saved to {output_file}")

def print_success_ascii():
    """
    Prints a bunny with a speech bubble saying 'SUCCESS'.
    """
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
    # Prompt for folder path
    input_folder = input("Enter the path to the folder containing Bracken report files: ")
    
    output_file = "rarefaction_curve_with_shannon.html"
    
    # Load Bracken files
    data = load_bracken_files(input_folder)
    rarefaction_results = {}

    # For each sample, calculate the rarefaction curve with Shannon Index
    for sample_id in data["sample_id"].unique():
        rarefaction_results[sample_id] = rarefaction_curve_with_shannon(data, sample_id)

    # Plot and save the rarefaction curves with Shannon Index as an HTML file
    plot_rarefaction_curves_html_with_shannon(rarefaction_results, output_file)

    # After plotting, print the success message with the ASCII bunny
    print_success_ascii()
