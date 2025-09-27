# Replication of Pulmonary Arterial Hypertension (PAH) Gene Expression Analysis

This repository contains an R script for a comprehensive bioinformatics analysis of Pulmonary Arterial Hypertension (PAH) using publicly available gene expression data from the Gene Expression Omnibus (GEO). The script replicates the study **Bioinformatics analysis of the GEO database for the identification of novel biomarkers and potential targeted drugs for pulmonary hypertension** by Lu et al. (2024) (https://doi.org/10.36922/gpd025080018)

The analysis integrates two GEO datasets (GSE117261 and GSE113439), corrects for batch effects, and performs several downstream analyses including differential gene expression, weighted gene co-expression network analysis (WGCNA), immune cell infiltration, and candidate drug screening. Also, a second script performs scRNA seq analysis using the hub genes discovered from the functional analyses and the 3 PAH samples from a different dataset, GSE169471.

## Workflow Overview

The analysis pipeline is structured as follows:

1.  **Data Acquisition and Pre-processing**:
    *   Downloads two PAH-related datasets (GSE117261 and GSE113439) from GEO.
    *   Harmonizes the datasets by finding common genes and standardizing phenotype information ("PAH" vs. "Normal").
    *   Merges the datasets and annotates microarray probes to official gene symbols.
    *   Visualizes batch effects using Principal Component Analysis (PCA).

2.  **Batch Effect Correction**:
    *   Applies the ComBat algorithm from the `sva` package to correct for technical differences between the two datasets.
    *   Visualizes the data post-correction via PCA to confirm the removal of the batch effect.

3.  **Differential Gene Expression (DEG) Analysis**:
    *   Uses the `limma` package to identify genes that are differentially expressed between PAH and normal samples.
    *   Visualizes the DEGs using a volcano plot and a heatmap of the most significant genes.

4.  **Weighted Gene Co-expression Network Analysis (WGCNA)**:
    *   Identifies modules of co-expressed genes from the top 25% most variable genes.
    *   Correlates these modules with the clinical trait (PAH vs. Normal).
    *   Identifies the most significant positive and negative correlated modules.
    *   Finds high-confidence "Module DEGs" (MDEGs) by intersecting the significant WGCNA modules with the list of DEGs.

5.  **Functional and Pathway Analysis**:
    *   Performs Gene Ontology (GO) and KEGG pathway enrichment analysis on the final MDEG list to understand their biological functions.
    *   Visualizes enriched terms and pathways.

6.  **Protein-Protein Interaction (PPI) Network and Hub Gene Identification**:
    *   Constructs a PPI network of the MDEGs using the STRING database.
    *   Identifies the top 10 hub genes based on network centrality (degree).
    *   Visualizes the PPI network and the hub genes.

7.  **Immune Infiltration Analysis**:
    *   Uses the CIBERSORT algorithm to estimate the relative proportions of 22 immune cell types from the bulk gene expression data.
    *   Visualizes immune cell proportions and compares them between PAH and normal samples.

8.  **Candidate Drug Screening**:
    *   Queries the Comparative Toxicogenomic Database (CTD) using an API call.
    *   Screens for chemical compounds that are known to affect the expression of the identified hub genes.
    *   Generates a ranked list of potential drug candidates for PAH.

9. **Single-Cell RNA-seq Analysis (scRNA-seq)**:

   * Downloads and loads the GSE169471 dataset, performs QC to filter out low-quality cells, and applies standard normalization and scaling procedures using the Seurat package.  
   * Utilizes PCA and UMAP for dimensionality reduction and visualization, followed by automated cell type annotation for each cluster using the SingleR package against a reference atlas.  
   * Validates the findings from the bulk RNA-seq analysis by creating feature plots and violin plots to visualize the expression of key hub genes across different identified cell types.  


## Prerequisites

To run this script, you need R and the following R packages installed. You can install them from CRAN or Bioconductor.

**CRAN Packages:**
```R
install.packages(c("dplyr", "tidyr", "ggplot2", "ggrepel", "VennDiagram", "ggvenn", "conflicted", "httr", "jsonlite", "pheatmap", "corrplot", "vioplot", "ggpubr"))
```

**Bioconductor Packages:**
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GEOquery", "Biobase", "affy", "limma", "sva", "ComplexHeatmap", "circlize", "clusterProfiler", "enrichplot", "org.Hs.eg.db", "hugene10sttranscriptcluster.db", "STRINGdb", "igraph", "ggraph", "WGCNA"))
```

## How to Run the Script

1.  **Clone the repository or save the script**: Download the R script (`PAH_analysis.Rmd`) to your local machine.
2.  **Set Working Directory**: Open the R script in RStudio or your preferred R environment and set the working directory to the script's location.
3.  **Execute the Script**: Run the script from top to bottom. The script is designed to be executed sequentially.

**Note on External Tools**:
*   **CIBERSORT**: The script is set up to load pre-computed CIBERSORT results from a CSV file. If you wish to run CIBERSORT yourself, you will need the `CIBERSORT.R` script and the `LM22.txt` signature matrix file, which can be obtained from the [CIBERSORT website](https://cibersort.stanford.edu/).
*   **NetworkAnalyst**: The regulatory network construction step is intended to be performed on the [NetworkAnalyst web server](https://www.networkanalyst.ca/), not directly in R.

## Output Files

The script will generate several plots and data files in subdirectories (`/plots/`, `/geo_files/`, etc.), including:

*   **PCA Plots**: `PCA_Before.png`, `PCA_After.png`
*   **DEG Visualizations**: `Volcano Plot of DEGs.png`, `Heatmap of Top Differentially Expressed Genes.png`
*   **WGCNA Plots**: `sample_outlier_dendrogram.png`, `power_selection_plots.png`, `Cluster Dendrogram.png`, `module_trait_heatmap.png`
*   **Functional Analysis Plots**: `GOplot.png`, `KEGG_dotplot2.png`
*   **PPI Network Plots**: `PPI_Network_Enhanced.png`, `Top10_Hub_Genes_Barplot.png`
*   **Data Files**: `Candidate_Drug_Interactions.csv` and intermediate files for CIBERSORT.
