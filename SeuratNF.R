## Setting wdir?

# Load packages:
library(Seurat)
## library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(optparse)

# Load input path
#data_dir <- "~/Masterarbeit/nextflow_spatial/work/6b/4375573436ed2285cb486c118deb0a/testrun/outs"

option_list <- list(
  make_option("--sample_id", type = "character", help = "Sample ID"),
  make_option("--input_dir", type = "character", help = "Spaceranger outs folder"),
  make_option("--output_dir", type = "character", help = "Output directory for results")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Load Seuratobject: Takes SpaceRanger output and returns seurat object (Contains spot-level expression and image of tissue)
#seurat_obj <- Load10X_Spatial(data.dir = data_dir)
seurat_obj <- Load10X_Spatial(data.dir = opt$input_dir)

# Create directory for result plots
dir.create("plots", showWarnings = FALSE)

### Data Preprocessing on spot by gene expression data
# Normalization: variance in sequencing depth (differences in cell density across the tissue)
# Violin-Plot: Number of RNA-molecules (UMIs) per spot (nCount_Spatial: Number of counts (genes x molecules) per spot)
plot1 <- VlnPlot(seurat_obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
# Show spatial sequencing depth
plot2 <- SpatialFeaturePlot(seurat_obj, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2) # Plot is only generated here!

# Save QC Plots
jpeg("./plots/qc_spatial_depth.jpeg", width = 2000, height = 1000, res = 200)
wrap_plots(plot1, plot2)
dev.off()

# Result: variance in molecular counts (technical as well as dependent on tissue anatomy)
# --> Standardnormalization might be problematic!, sctransfrom as alternativ

# builds regularized negative binomial models of gene expression in order to account for technical artifacts while preserving biological variance
# normalizes the data, detects high-variance features, and stores the data in the SCT assay.

seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE)

install.packages('BiocManager') # for faster estimation
BiocManager::install('glmGamPoi')

# Gene expression visualisation: overlay molecular data on top of tissue histology
SpatialFeaturePlot(seurat_obj, features = c("Hpca", "Ttr"))

# Dimensionality reduction, clustering and visualization
seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:30)

# DimPlot: clustering in UMAP space
jpeg("./plots/clustering_UMAP.jpeg", width = 2000, height = 1000, res = 200)
clustering_UMAP <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
print(clustering_UMAP)
dev.off()


# SpatialDimPlot: Clustering overlaid on the image (label!)
jpeg("./plots/clustering_image_overlay.jpeg", width = 2000, height = 1000, res = 200)
clustering_image_overlay <- SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3)
print(clustering_image_overlay)
dev.off()

# Distinguish location of individual clusters (takes a few seconds)
jpeg("./plots/individual_clusters_on_image.jpeg", width = 2000, height = 1000, res = 200)
# Distinguish location of individual clusters (takes a few seconds)
individual_clusters_on_image <- SpatialDimPlot(seurat_obj, cells.highlight = CellsByIdentities(object = seurat_obj, idents = c(2, 1, 4, 3,
                                                                                           5, 8)), facet.highlight = TRUE, ncol = 3)
print(individual_clusters_on_image)
dev.off()

## Show ALL clusters
# Define clustes
# Idents(seurat_obj) <- "seurat_clusters"
# all_clusters <- levels(Idents(seurat_obj))

# jpeg("./plots/individual_clusters_on_image.jpeg", width = 2000, height = 1000, res = 200)
#individual_clusters_on_image <- SpatialDimPlot(
#   seurat_obj,
#   cells.highlight = CellsByIdentities(object = seurat_obj, idents = all_clusters),
#   facet.highlight = TRUE,
#   ncol = 3
# )
# print(individual_clusters_on_image)
# dev.off()

# Identification of Spatially Variable Features 
# presto package for faster implementation of Wilcoxon Rank Sum Test, automatically downloaded

# Differential expression based on pre-annotated anatomical regions (Either from unsupervised clustering or prior knowledge)
# Compares spots of two clusters (eg 5 and 6) --> de_markes is a dataframe, which contains p-Value,, Log2 Foldchange, % cells,...
# Goal: finding genes which are expressed more in Cluster 5 or 6
de_markers <- FindMarkers(seurat_obj, ident.1 = 5, ident.2 = 6)

# Plotting where these genes are expressedSpatialFeaturePlot() zeichnet, wo im Gewebe diese Gene exprimiert werden.
# features = rownames(de_markers)[1:3] takes first three genesnimmt die Top 3 Gene (also die ersten drei Gene im Ergebnis von FindMarkers()).
# alpha = c(0.1, 1) for transparency: low expression: 0.1, highly expressed: opaque
# ncol = 3 (3 Plots next to each other)
jpeg("./plots/Spat_Var_Features_preannotated.jpeg", width = 2000, height = 1000, res = 200)
Spat_Var_Features_preannotated <- SpatialFeaturePlot(object = seurat_obj, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
print(Spat_Var_Features_preannotated)
dev.off()


# Alternative without pre-annotation:
# Automatic installation: install.packages('Rfast2') for efficient implementation of Morans I calculation
# Moran’s I: finding genes with non-coincidental expression (clustered, patterns)
# Identify spatially relevant genes without knowing tissue-regions
seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, assay = "SCT", features = VariableFeatures(seurat_obj)[1:1000],
                                       selection.method = "moransi")

# Official Seuratcode doesnt work here, would be:
#top.features <- head(SpatiallyVariableFeatures(brain_x, selection.method = "moransi"), 6)
#SpatialFeaturePlot(brain_x, features = top.features, ncol = 3, alpha = c(0.1, 1))
# Gives back: Fehler in h(simpleError(msg, call)) : 
# Fehler bei der Auswertung des Argumentes 'x' bei der Methodenauswahl für Funktion 'head': kann nicht xtfrm Dataframe

# Workaround (results same as seurat-vignette)
# Get Metadata and genes as column
svf <- seurat_obj@assays$SCT@meta.features
svf$gene <- rownames(seurat_obj@assays$SCT@meta.features)

# IFlter spatial variable genes
svf <- svf[!is.na(svf$moransi.spatially.variable.rank), ]

# Order top genes
svf <- svf[order(svf$moransi.spatially.variable.rank), ]
top.features <- head(svf$gene, 6)

# Visualisation
jpeg("./plots/Spatial_variable_features_top6.jpeg", width = 2000, height = 1000, res = 200)
Spatial_variable_features_top6 <- SpatialFeaturePlot(seurat_obj, features = top.features, ncol = 3, alpha = c(0.1, 1))
print(Spatial_variable_features_top6)
dev.off()

# Additional: interactive plotting, subset of anatomical regions, wokrking with mulitple slices
# --> Bioinformatik/Masterarbeit/Seurat/Seurat.R


                                                                                                                                                                                                                                            


