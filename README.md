# Guitools

**Guitools** is a Python-based toolkit designed to process Bracken output files and generate comprehensive taxonomic abundance visualizations. It supports various analyses, including heatmaps, dendrograms, stacked bar plots, beta diversity (PCoA), bubble graphs, and rarefaction curve analyses. This makes it ideal for microbial ecology and metagenomics studies.

## Features
- **Heatmaps**: Taxonomic abundance across samples.
- **Dendrograms**: Sample clustering based on Bray-Curtis dissimilarity.
- **Stacked Bar Plots**: Visualize relative abundance across taxonomic groups.
- **Beta Diversity (PCoA)**: Bray-Curtis dissimilarity for visualizing compositional differences.
- **Bubble Graphs**: Compare taxa abundances.
- **Rarefaction Curves**: Assess alpha diversity based on sequencing depth.

## Installation

### Prerequisites
- Python 3.8 or higher
- Required libraries: `pandas`, `glob`, `os`, `plotly`, `matplotlib`, `seaborn`, `scipy`, `skbio`

Install dependencies:

```bash
pip install pandas plotly matplotlib seaborn scipy scikit-bio
```
Clone the Repository
```bash
git clone https://github.com/GmoreiraVet/Guitools.git
cd Guitools
```
Usage
Preparing Input Files
Input Folder: Place Bracken report files in the designated directory. Default path:
```python
input_folder = "/home/viroicbas/scriptTeste/bracken_reports/"
```
Modify input_folder in the scripts if needed.
File Naming: Ensure Bracken report files follow the naming format SAMPLE_bracken.txt

##Outputs
Each script saves visualizations as HTML or PNG files in the working directory, ready for publication.

# :books:<span style="color: green;"> Citation :books:

If you use Guitools in your research, please consider citing the software. You can generate the citation directly from GitHub by using the "Cite this repository" button located in the repository's header (next to the "Fork" and "Star" buttons). This will generate a citation formatted in the recommended style for GitHub repositories.

Thank you for your support!




