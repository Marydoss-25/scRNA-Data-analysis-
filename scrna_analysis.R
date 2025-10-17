setwd("../scrna/Analysis/")

library(Seurat)
library(SeuratDisk)
library(tidyverse)

#Read10x data
nsclc_raw <- Read10X_h5("../Analysis/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
nsclc_cts <- nsclc_raw$`Gene Expression`

#Createseuratobject
nsclc.seurat <- CreateSeuratObject(counts = nsclc_cts, project = "NSCLC" , min.cells = 3, min.features = 200)

View(nsclc.seurat@meta.data)
range(nsclc.seurat$nCount_RNA)
range(nsclc.seurat$nFeature_RNA)

#QC - %MT 
nsclc.seurat[["percent.MT"]] <- PercentageFeatureSet(nsclc.seurat,pattern = "^MT-")
VlnPlot(nsclc.seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.MT"))
FeatureScatter(nsclc.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")

#Filtering

nsclc.seurat <- subset(nsclc.seurat, subset = nFeature_RNA >200 & nFeature_RNA < 2500 & percent.MT < 5)

#Normalizing

nsclc.seurat <- NormalizeData(nsclc.seurat, normalization.method = "LogNormalize" , scale.factor = 10000)

#Filtering High variable genes

nsclc.seurat <- FindVariableFeatures(nsclc.seurat, selection.method = "vst", nfeatures = 2000)

#selecting top 10 hgh variable features

top10 <- head(VariableFeatures(nsclc.seurat), 10)

#ploting variable features
plot1 <- VariableFeaturePlot(nsclc.seurat)
LabelPoints(plot1, points = top10, repel = TRUE)
  
#Scaling the data
allgenes <- rownames(nsclc.seurat)
nsclc.seurat <- ScaleData(nsclc.seurat, features = allgenes )

#Dimensional reduction
nsclc.seurat <- RunPCA(nsclc.seurat, features = VariableFeatures(nsclc.seurat))
print(nsclc.seurat[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat, dims = 1 , cells = 500, balanced = TRUE)

#Determine the dimesnionality reduction
ElbowPlot(nsclc.seurat)

#Clustering
nsclc.seurat <- FindNeighbors(nsclc.seurat, dims = 1:15)
nsclc.seurat <- FindClusters(nsclc.seurat, resolution = c(0.1,0.3,0.5,0.7,0.9,1))
DimPlot(nsclc.seurat, group.by = "RNA_snn_res.0.1", label = TRUE)
  
#Setting the identity
Idents(nsclc.seurat)
Idents(nsclc.seurat) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat)

#UMAP visualization
nsclc.seurat <- RunUMAP(nsclc.seurat, dims = 1:15)
DimPlot(nsclc.seurat, reduction = "umap", label = TRUE)

