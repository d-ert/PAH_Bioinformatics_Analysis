# Replication of Pulmonary Arterial Hypertension (PAH) GEO Data Analysis

## Project Overview

This R script performs a comprehensive bioinformatic analysis to identify key genes, biological pathways, and potential therapeutic targets for Pulmonary Arterial Hypertension (PAH). The analysis integrates data from two distinct GEO datasets (GSE117261 and GSE113439), harmonizes them, and applies a multi-step workflow including differential gene expression, weighted gene co-expression network analysis (WGCNA), functional enrichment, immune infiltration analysis, and candidate drug screening.

## Analysis Workflow

The script follows a structured pipeline to process and analyze the gene expression data:

1.  **Data Acquisition and Pre-processing:**
    *   Downloads two PAH-related gene expression datasets (GSE117261 and GSE113439) from the Gene Expression Omnibus (GEO).
    *   Extracts expression and phenotype data.
    *   Identifies common genes between the two datasets and standardizes the phenotype data ("PAH" vs. "Normal").

2.  **Batch Effect Correction:**
    *   Merges the two datasets into a single expression matrix.
    *   Visualizes the batch effect between the two datasets using Principal Component Analysis (PCA).
    *   Applies the ComBat algorithm from the `sva` package to correct for batch effects while preserving the biological differences between PAH and normal samples.
    *   Visualizes the data post-correction via PCA to confirm the removal of the batch effect.

3.  **Differential Gene Expression (DEG) Analysis:**
    *   Uses the `limma` package to identify Differentially Expressed Genes (DEGs) between PAH and normal samples.
    *   Visualizes the DEGs using a volcano plot and a heatmap.

4.  **Weighted Gene Co-expression Network Analysis (WGCNA):**
    *   Constructs a gene co-expression network to identify modules of highly correlated genes.
    *   Selects an appropriate soft-thresholding power to achieve a scale-free network topology.
    *   Identifies gene modules and correlates them with the clinical trait (PAH).
    *   Identifies the modules most significantly associated (both positively and negatively) with PAH.

5.  **Hub Gene Identification and Functional Analysis:**
    *   Identifies high-confidence "Module DEGs" (MDEGs) by finding the intersection between the DEGs and the genes within the most significant WGCNA modules.
    *   Performs Gene Ontology (GO) and KEGG pathway enrichment analysis on the final MDEG list to understand their biological functions.
    *   Constructs a Protein-Protein Interaction (PPI) network using the STRING database to identify the top 10 hub genes based on network centrality.

6.  **Immune Infiltration Analysis:**
    *   Uses the CIBERSORT algorithm to estimate the relative proportions of 22 immune cell types from the bulk gene expression data.
    *   Visualizes the immune cell composition and compares the proportions between PAH and normal groups.

7.  **Candidate Drug Screening:**
    *   Queries the Comparative Toxicogenomic Database (CTD) with the identified hub genes.
    *   Screens for chemical compounds that are known to affect the expression of these genes, identifying potential therapeutic candidates for PAH.

## How to Run the Script

1.  **Install Dependencies:** Ensure all the required R packages are installed. You can install them using the following command in R:
    ```R
    install.packages(c("dplyr", "GEOquery", "Biobase", "affy", "limma", "sva", "tidyr", "ggplot2", "ComplexHeatmap", "circlize", "clusterProfiler", "enrichplot", "ggrepel", "VennDiagram", "org.Hs.eg.db", "hugene10sttranscriptcluster.db", "STRINGdb", "igraph", "ggraph", "ggvenn", "WGCNA", "httr", "jsonlite", "corrplot", "vioplot", "ggpubr"))
    ```
2.  **Set Working Directory:** Set your R working directory to the location where the script is saved.
3.  **Run Script:** Execute the R script in your R environment (e.g., RStudio). The script will automatically create directories for plots and downloaded data.

## Dependencies

This analysis relies on the following R packages:

*   **Data Handling & Manipulation:** `dplyr`, `tidyr`
*   **Bioinformatics & Genomics:** `GEOquery`, `Biobase`, `affy`, `limma`, `sva`, `clusterProfiler`, `org.Hs.eg.db`, `hugene10sttranscriptcluster.db`, `WGCNA`
*   **Visualization:** `ggplot2`, `ComplexHeatmap`, `circlize`, `enrichplot`, `ggrepel`, `VennDiagram`, `ggraph`, `ggvenn`, `corrplot`, `vioplot`, `ggpubr`
*   **Network Analysis:** `STRINGdb`, `igraph`
*   **API/Web:** `httr`, `jsonlite`
*   **Utilities:** `conflicted`

## Output Files

The script will generate several output files, saved in dedicated subdirectories (`plots/`, `cibersort_results/plots/`):

*   **PCA Plots:** Before and after batch effect correction.
*   **DEG Visualizations:** Volcano plot and heatmap of top DEGs.
*   **WGCNA Plots:** Sample dendrogram, scale-free topology plots, and module-trait relationship heatmap.
*   **Venn Diagram:** Showing the intersection of DEGs and significant WGCNA module genes.
*   **MDEG Heatmap:** Heatmap of the final high-confidence genes.
*   **Enrichment Plots:** Bar plots for GO and KEGG enrichment results.
*   **PPI Network Graphs:** Network visualization of key PAH genes and a bar chart of top hub genes.
*   **Immune Infiltration Plots:** Heatmap, violin plots, and correlation matrix of immune cell subtypes.
*   **Candidate Drugs CSV:** A file named `Candidate_Drug_Interactions.csv` containing a list of potential drugs and their interactions with the hub genes.
