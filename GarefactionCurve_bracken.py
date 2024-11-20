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
        # Load Bracken data, focusing on taxonomy_id and new_est_reads columns
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

    # Simulate subsampling at increasing depths
    total_reads = sample_data["new_est_reads"].sum()
    for depth in range(10, total_reads, 10):
        # Subsample by sampling "depth" reads (randomly selecting from the available reads)
        subsample = sample_data.sample(n=depth, replace=True, weights=sample_data["new_est_reads"])
        unique_taxa = subsample["taxonomy_id"].nunique()  # Count unique OTUs
        rarefaction_data[depth] = unique_taxa

    return rarefaction_data

def log_model(x, a, b):
    """
    Logarithmic model: y = a * log(x) + b
    """
    return a * np.log(x) + b

def plot_rarefaction_curves_html(rarefaction_results, output_file="rarefaction_curve_with_log_fit.html"):
    """
    Plot rarefaction curves for all samples as an interactive HTML, with a smoothed logarithmic-like curve.

    Args:
    - rarefaction_results (dict): Dictionary of sample_id and their rarefaction data.
    - output_file (str): Path to save the rarefaction curve HTML file.
    """
    fig = go.Figure()

    # Define color scale (to assign distinct colors for each sample)
    color_scale = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    color_idx = 0  # Start index for colors
    sample_names = {}  # Dictionary to store user-defined names for each sample

    # Initialize the variable to track the max OTUs across all samples
    max_otus = 0

    # Plot individual rarefaction curves with smoothed log-like curve
    for sample_id, data in rarefaction_results.items():
        # Prompt the user to provide a name for each sample
        sample_name = input(f"Enter a name for sample '{sample_id}': ")
        sample_names[sample_id] = sample_name

        depths = list(data.keys())
        species_counts = list(data.values())

        # Fit the rarefaction data to a logarithmic model
        popt, _ = curve_fit(log_model, depths, species_counts, maxfev=10000)

        # Generate smoothed curve using the fitted logarithmic model
        smooth_depths = np.linspace(min(depths), max(depths), 500)
        smooth_counts = log_model(smooth_depths, *popt)

        # Pick color for the current sample's line and shaded area
        color = color_scale[color_idx % len(color_scale)]  # Loop through the color scale

        # Plot the smoothed logarithmic curve for each sample
        fig.add_trace(go.Scatter(x=smooth_depths, y=smooth_counts, mode='lines', name=sample_name,
                                 line=dict(color=color)))

        # Plot the shaded area representing standard deviation (approximated as range)
        std_counts = np.std(species_counts)  # Standard deviation (though this is an approximation)

        # Using RGBA format for fillcolor (with alpha transparency)
        fill_color = f"rgba{tuple(int(color[i:i+2], 16) for i in (1, 3, 5)) + (0.2,)}"  # Convert hex to RGBA

        fig.add_trace(go.Scatter(
            x=smooth_depths.tolist() + smooth_depths[::-1].tolist(),  # Reverse the depths for the shaded area
            y=[count + std_counts for count in smooth_counts] + [count - std_counts for count in smooth_counts[::-1]],  # Upper and lower bounds for the shaded area
            fill='toself',
            fillcolor=fill_color,  # Apply RGBA color
            line=dict(color='rgba(255, 255, 255, 0)'),  # Hide the border line
            name=f"SD Area: {sample_name}"
        ))

        # Track the maximum number of OTUs for dynamic Y-axis scaling
        max_otus = max(max_otus, max(species_counts))

        color_idx += 1  # Increment color index

    # Add a user prompt to input a custom title for the plot
    plot_title = input("Enter a title for the rarefaction curve plot: ")

    # Set the Y-axis range based on the maximum OTUs found
    fig.update_layout(
        title=plot_title,  # Use the user-defined title
        xaxis_title="Number of Reads Sampled",
        yaxis_title="Unique OTUs (Taxa)",
        legend_title="Samples",
        template="plotly_white",
        yaxis=dict(range=[0, max_otus * 1.2])  # Set Y-axis range based on max OTUs (with a 20% margin)
    )

    fig.write_html(output_file)
    fig.show()
    print(f"Rarefaction curve saved to {output_file}")

if __name__ == "__main__":
    input_folder = "/home/viroicbas/scriptTeste/bracken_reports"  # Update with your folder containing Bracken files
    output_file = "rarefaction_curve_with_log_fit.html"
    
    # Load Bracken files
    data = load_bracken_files(input_folder)
    rarefaction_results = {}

    # For each sample, calculate the rarefaction curve
    for sample_id in data["sample_id"].unique():
        rarefaction_results[sample_id] = rarefaction_curve(data, sample_id)

    # Plot and save the rarefaction curves as an HTML file with the logarithmic fit and standard deviation
    plot_rarefaction_curves_html(rarefaction_results, output_file)

