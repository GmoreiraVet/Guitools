import os
import glob
import pandas as pd
import numpy as np
import plotly.express as px
import math
from scipy.optimize import root_scalar

# Diversity calculations
def shannons_alpha(p):
    return -1 * sum(i * math.log(i) for i in p if i > 0)

def simpsons_index(D):
    return 1 - D

def fishers_alpha(counts):
    counts = counts[counts > 0]
    S = len(counts)  # Number of species
    N = counts.sum()  # Total number of individuals
    
    if N == 0 or S == 0:
        return 0
    
    # Solve for alpha in S = alpha * ln(1 + N/alpha)
    def equation(alpha):
        return alpha * math.log(1 + N/alpha) - S
    
    try:
        solution = root_scalar(equation, bracket=[1e-10, 1e6], method='brentq')
        if solution.converged:
            return solution.root
        else:
            return 0
    except ValueError:
        return 0

def calculate_diversities(counts):
    counts = counts[counts > 0]
    N = counts.sum()
    if N < 2:
        return {"Shannon Index": 0, "Simpson Index (1-D)": 0, "Fisher's Alpha": 0}
    
    p = (counts / N).tolist()
    D = (counts * (counts - 1)).sum() / (N * (N - 1))
    
    return {
        "Shannon Index": shannons_alpha(p),
        "Simpson Index (1-D)": simpsons_index(D),
        "Fisher's Alpha": fishers_alpha(counts)
    }

def load_bracken_files_with_custom_names(input_folder):
    file_paths = glob.glob(os.path.join(input_folder, "*.txt"))
    all_data = []

    print("\nDetected Bracken files:")
    sample_name_map = {}
    for file in file_paths:
        default_name = os.path.basename(file).replace("_bracken.txt", "")
        new_name = input(f"Enter a custom sample name for '{default_name}' (press Enter to keep default): ").strip()
        sample_name_map[file] = new_name if new_name else default_name

    for file in file_paths:
        sample_id = sample_name_map[file]
        try:
            df = pd.read_csv(file, sep="\t")
            if "new_est_reads" not in df.columns:
                continue
            counts = df["new_est_reads"]
            diversities = calculate_diversities(counts)
            for index, value in diversities.items():
                all_data.append({
                    "Sample": sample_id,
                    "Index": index,
                    "Value": value
                })
        except Exception as e:
            print(f"Error processing {file}: {e}")
            continue

    return pd.DataFrame(all_data)

def plot_diversity(df, title=None):
    fig = px.bar(
        df,
        x="Sample",
        y="Value",
        color="Index",
        title=title if title else "Alpha Diversity Indices per Sample",
        template="plotly_white",
        barmode="group",
        color_discrete_sequence=px.colors.qualitative.Pastel,
    )
    
    fig.update_traces(width=0.4)
    fig.update_layout(
        yaxis_title="Diversity Value",
        xaxis_title="Sample ID",
        legend_title="Diversity Index",
        plot_bgcolor='white'
    )
    fig.update_xaxes(tickangle=-45)
    return fig

if __name__ == "__main__":
    input_folder = input("Enter the path to the folder containing Bracken files: ").strip()
    df = load_bracken_files_with_custom_names(input_folder)

    if df.empty:
        print("No valid Bracken files found.")
    else:
        custom_title = input("\nEnter a custom title for the plot (press Enter for default): ").strip()
        fig = plot_diversity(df, title=custom_title)
        fig.show()
