# scRNA-Data-analysis-
Standard workflow of scRNA data analysis using R
-------------------------------------

This project contains a single R script that performs end-to-end single-cell RNA sequencing (scRNA-seq) analysis using the Seurat package in R.
The analysis covers preprocessing, dimensionality reduction (PCA, UMAP), clustering, and visualization.

-------------------------------------------------------
Project Contents
-------------------------------------------------------
scrna_analysis.R   - Complete Seurat workflow (from data loading to visualization)
README.txt         - Project documentation
Dataset	used        - Non small cell lung cancer sample data - 20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5
-------------------------------------------------------
Requirements
-------------------------------------------------------
R version: 4.2 or higher

R Packages:
  - Seurat
  - SeuratDisk
  - tidyverse
  - reticulate

Python package (via reticulate):
  - umap-learn

To install required packages in R:
  install.packages(c("Seurat", "SeuratDisk", "tidyverse", "reticulate"))
  reticulate::py_install("umap-learn")

-------------------------------------------------------
How to Run
-------------------------------------------------------
1. Open the script in R or RStudio.
2. Run the following command in R:
   source("scrna_analysis.R")

   OR, if using an R Markdown file:
   rmarkdown::render("nsclc_scrna_analysis.Rmd")

3. The script performs:
   - Data preprocessing
   - PCA and clustering
   - UMAP dimensionality reduction
   - Visualization of cell clusters

-------------------------------------------------------
Standard Workflow of scRNA-seq Analysis
-------------------------------------------------------
1. Upload scRNA data  
2. Convert to Seurat object  
3. Quality Check  
   3.1. Measure % mitochondrial content 
   3.2. Visualization:  
        - Violin plots to inspect feature distributions  
        - Scatter plots to examine feature correlations  
4. Filtering  
   - Remove low-quality cells   
   - Typical thresholds:  
        * min.features > 200  
        * min.cells > 3  
        * percent.mt < 5  
5. Normalization  
6. Identify highly variable genes  
7. Scale data  
8. Dimensionality reduction (PCA)  
9. Clustering  
10. Non-linear dimensionality reduction (UMAP or t-SNE)  
11. Visualization of clusters  

-------------------------------------------------------
Author
-------------------------------------------------------
Arpudhamary V.
