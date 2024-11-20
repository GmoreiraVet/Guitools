#!/usr/bin/env python
import os
import sys
import argparse
import numpy as np
import pandas as pd
from skbio.stats.ordination import pcoa
from skbio.stats.distance import DistanceMatrix
import plotly.express as px

def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', '--input-files', required=True, dest='in_files', nargs='+',
                        help='Input files (one per community) for which to compare for Bray-Curtis dissimilarity metrics')
    parser.add_argument('--type', required=False, default='single', dest='filetype', choices=['single', 'simple', 'bracken', 'kreport', 'kreport2', 'krona'],
                        help='Type of input file[s]: single, simple, bracken, kreport, kreport2, krona. See docs for details')
    parser.add_argument('--cols', '--columns', dest='cols', required=False, default='1,2',
                        help='Specify category/counts separated by single comma: cat,counts (1 = first col)')
    parser.add_argument('--level', '-l', dest='lvl', required=False, default='all', choices=['all', 'S', 'G', 'F', 'O'],
                        help='For Kraken or Krona files, taxonomy level for which to compare samples. Default: all')
    args = parser.parse_args()

    # Step 1: Read in Samples from Files
    i2totals = {}
    i2counts = {}
    i2names = {}
    num_samples = 0
    num_categories = 0
    categ_col, count_col = map(int, args.cols.split(','))
    categ_col -= 1
    count_col -= 1

    for f in args.in_files:
        if not os.path.isfile(f):
            sys.stderr.write(f"File {f} not found\n")
            exit(1)

        i_file = open(f, 'r')
        i2totals[num_samples] = 0
        i2counts[num_samples] = {}

        # Ask user for a custom name for each sample
        sample_name = input(f"Enter a name for the sample '{f}': ")
        i2names[num_samples] = sample_name

        for line in i_file:
            l_vals = line.strip().split("\t")
            if len(l_vals) == 0 or (not l_vals[count_col].isdigit()) or l_vals[0] == '#':
                continue

            curr_categ = l_vals[categ_col]
            count = int(l_vals[count_col])
            i2totals[num_samples] += count
            if curr_categ not in i2counts[num_samples]:
                i2counts[num_samples][curr_categ] = 0
            i2counts[num_samples][curr_categ] += count
            num_categories += 1
        i_file.close()
        num_samples += 1

    # Step 2: Calculate Bray-Curtis Dissimilarities
    bc = np.zeros((num_samples, num_samples))
    for i in range(0, num_samples):
        i_tot = i2totals[i]
        for j in range(i + 1, num_samples):
            j_tot = i2totals[j]
            C_ij = 0.0
            for cat in i2counts[i]:
                if cat in i2counts[j]:
                    C_ij += min(i2counts[i][cat], i2counts[j][cat])
            # Calculate Bray-Curtis dissimilarity
            bc_ij = 1.0 - ((2.0 * C_ij) / float(i_tot + j_tot))
            bc[i][j] = bc_ij
            bc[j][i] = bc_ij

    # Step 3: Perform PCoA
    dist_matrix = DistanceMatrix(bc.tolist(), list(range(num_samples)))
    pcoa_results = pcoa(dist_matrix)

    # Extract PCoA results
    pcoa_df = pd.DataFrame({
        'Sample': [i2names[i] for i in range(num_samples)],  # Use custom sample names
        'PC1': pcoa_results.samples['PC1'],
        'PC2': pcoa_results.samples['PC2'],
        'PC3': pcoa_results.samples['PC3'],
        'PC4': pcoa_results.samples['PC4'],
        'PC5': pcoa_results.samples['PC5'],
    })

    # Step 4: Add jitter to PC1 and PC2 values
    jitter_strength = 0.01  # Adjust this value for more or less jitter
    pcoa_df['PC1'] += np.random.uniform(-jitter_strength, jitter_strength, size=len(pcoa_df))
    pcoa_df['PC2'] += np.random.uniform(-jitter_strength, jitter_strength, size=len(pcoa_df))

    # Step 5: Calculate variance explained for each PC
    explained_variance = pcoa_results.proportion_explained * 100  # Multiply by 100 for percentage
    pc1_variance = explained_variance[0]
    pc2_variance = explained_variance[1]

    # Define the axis limits for zooming (smaller scale by reducing padding)
    x_range = [pcoa_df['PC1'].min() - 0.05, pcoa_df['PC1'].max() + 0.05]  # Reduced padding
    y_range = [pcoa_df['PC2'].min() - 0.05, pcoa_df['PC2'].max() + 0.05]  # Reduced padding

    # Plot using Plotly with text positioned below dots
    fig = px.scatter(pcoa_df, x='PC1', y='PC2', text='Sample', color='Sample', 
                     title=f'Principal Coordinate Analysis (PCoA) of Bray-Curtis Dissimilarity\nPC1: {pc1_variance:.2f}% variance, PC2: {pc2_variance:.2f}% variance',
                     labels={'PC1': 'PC1', 'PC2': 'PC2'})

    fig.update_traces(marker=dict(size=12, line=dict(width=2, color='DarkSlateGrey')),
                      textposition='bottom center')  # Offset the text position below the points
    fig.update_layout(
        showlegend=False, 
        margin=dict(t=50, b=50, l=50, r=50),
        xaxis_title=f'PC1 ({pc1_variance:.2f}%)',
        yaxis_title=f'PC2 ({pc2_variance:.2f}%)',
        xaxis=dict(range=x_range),  # Set x-axis range (smaller scale)
        yaxis=dict(range=y_range),   # Set y-axis range (smaller scale)
        width=1000,   # Slightly larger plot width
        height=800    # Slightly larger plot height
    )

    # Show the plot
    fig.show()

if __name__ == "__main__":
    main()

