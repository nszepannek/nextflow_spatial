#!/usr/bin/env Rscript

# Load packages:
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(zellkonverter)
library(SingleCellExperiment)
library(org.Mm.eg.db)
library(SingleR)
library(reticulate)

# Load input path: takes spaceranger output as input
args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]

# Load Seuratobject: takes spaceranger output and returns seurat object (contains spot-level expression and image of tissue)
seurat_obj <- Load10X_Spatial(data.dir = data_dir)

# Create directory for result plots and csv
dir.create("plots_qc", showWarnings = FALSE)
dir.create("csv", showWarnings = FALSE)


### Data Preprocessing on spot by gene expression data
# Normalization: variance in sequencing depth (differences in cell density across the tissue)
# Violin-Plot: Number of RNA-molecules (UMIs) per spot (nCount_Spatial: Number of counts (genes x molecules) per spot) - distribution
plot1 <- VlnPlot(seurat_obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
# Show spatial sequencing depth: feature nCount_Spatial is shown spatially
plot2 <- SpatialFeaturePlot(seurat_obj, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2) # Plot is only generated here!

# Save QC Plots together
jpeg("./plots_qc/qc_spatial_depth.jpeg", width = 2000, height = 1000, res = 200)
wrap_plots(plot1, plot2)
dev.off()

# Variance in molecular counts (technical as well as dependent on tissue anatomy)
# --> Standardnormalization might be problematic!, sctransfrom as alternativ
# more sturdy than Normalize + Scale

# builds regularized negative binomial models of gene expression in order to account for technical artifacts while preserving biological variance
# normalizes the data, detects high-variance features, and stores the data in the SCT assay.

# Important for spatial data which variance in sequencing depth as well as spot to spot variability
# Needed for PCA, Clustering, DEG, as these need consistent and scaled data

seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE) # verbose: no output in console

# Dimensionality reduction, clustering and visualization
seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:30)

# Save the seurat object for further analysis with 
saveRDS(seurat_obj, file = "./seurat_obj_with_umap.rds")

# Extract UMAP coordinates for every spot (with cluster identity), save it
umap_df <- as.data.frame(Embeddings(seurat_obj, "umap"))
umap_df$cluster <- Idents(seurat_obj)
write.csv(umap_df, file = "umap_coords_with_clusters.csv", quote = FALSE)
