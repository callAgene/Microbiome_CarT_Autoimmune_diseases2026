# Custom Analysis Scripts for: [Microbiome_recapitulation_afte_-CAR-T_therapy_Autoimmune_diseases]


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview
This repository contains the collection of custom R and Python scripts used for the data analysis and figure generation in our manuscript. 

Rather than a single automated pipeline, these are standalone scripts categorized by omics layer:
1. **scRNA-seq**: Scripts for single-cell clustering and visualization (`scanpy`).
2. **Metagenomics**: PathwayDF.R Scripts for differential abundance and GSEA (`limma`, `clusterProfiler`).
3. **Metabolomics**: metabolism.R Scripts for data preprocessing, scaling, and PCA.
4. **16S rRNA Microbiome**: lefseFor16Sfull.R Scripts for biomarker discovery (`LEfSe`).


A detailed description of the statistical methods applied in these scripts can be found in the **Methods** section of the main manuscript.

---

## 1. System Requirements

The scripts are intended to be run interactively in standard R and Python environments.

### Software Requirements
* **Operating Systems:** macOS (tested on Sonoma 14.2), Linux (Ubuntu 22.04 LTS), or Windows 11.
* **R Environment (v4.2.0+):** `dplyr`, `ggplot2`, `limma`, `clusterProfiler`, `lefser`, `SummarizedExperiment`.
* **Python Environment (v3.9+):** `scanpy`, `anndata`, `pandas`, `numpy`, `matplotlib`, `seaborn`.

---

## 2. Installation Guide

For the R scripts, open R/RStudio and run:
```R
# Install core and Bioconductor dependencies
install.packages(c("dplyr", "ggplot2", "openxlsx", "factoextra"))
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "clusterProfiler", "lefser", "SummarizedExperiment"))
```


For the Python single-cell scripts, we recommend using conda:
```BASH
conda create -n sc_env python=3.9
conda activate sc_env
pip install scanpy anndata pandas numpy matplotlib seaborn
```

## 3. Demo (scRNA-seq Analysis)
To fulfill the code reproducibility requirements, we have provided a 1% randomly downsampled subset of our single-cell RNA-seq (scRNA-seq) data. This small, anonymized dataset demonstrates the execution of our core preprocessing, dimensionality reduction, and clustering scripts while remaining lightweight enough to run quickly on a standard computer.

### Instructions to Run
1. Download this repository and set it as your working directory.
2. Ensure you have activated the required Python environment (see Section 2).
   ```bash
   conda activate sc_env
   ```
3. Locate the demo script and data in the ./Demo_scRNAseq/ folder.

4. Run the script via your command line:
    ```Bash
    python ./Demo/scRNAseq_demo.py
    ```
5. Expected Output
The script will successfully execute the Scanpy workflow and output a .h5ad processed object along with a UMAP visualization plot (umap_demo_clusters.pdf) in the same directory.

(Note: Because this demo uses only 1% of the total cells, the resulting UMAP will appear sparse and the cluster resolutions will differ from the main figures in the manuscript. The purpose of this demo is solely to validate the functionality of the code infrastructure.)
(Note: The full scRNA-seq analysis script requires a machine with at least 32GB of RAM and may take over an hour to compute on the full dataset.)

 6. Expected run time: < 2 minutes on a normal desktop computer.


   
