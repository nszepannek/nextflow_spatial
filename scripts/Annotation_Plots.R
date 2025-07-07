#!/usr/bin/env Rscript

# Load packages:
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(zellkonverter)
library(pheatmap)
library(SingleCellExperiment)
library(org.Mm.eg.db)
library(SingleR)
library(viridis)
library(reticulate)

args <- commandArgs(trailingOnly = TRUE)
seurat_obj_path <- args[1]
annot_ref <- args[2]

# Create directorys for plots and csv
dir.create("plots_annotation", showWarnings = FALSE)
if (!dir.exists("csv")) dir.create("csv")

# Read in annotated reference
annot_ref_h5ad <- sub("\\.rds$", ".h5ad", annot_ref)

if (file.exists(annot_ref_h5ad)) {
  ref_h5ad <- readH5AD(annot_ref_h5ad)

  # Adapt Ensemble ID and Gene ID
  symbols <- mapIds(
   org.Mm.eg.db,
   keys = rownames(ref_h5ad),
   column = "SYMBOL",
   keytype = "ENSEMBL",
   multiVals = "first"
   )

  # Rename
  rownames(ref_h5ad) <- symbols
  saveRDS(ref_h5ad, annot_ref)
}

ref_rds <- readRDS(annot_ref)

# Read in seurat-object with UMAP
seurat_obj <- readRDS(seurat_obj_path)

seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE)

#### Convert seurat object to single cell experiment
class(as.SingleCellExperiment(seurat_obj))
seurat_obj.sce <- as.SingleCellExperiment(seurat_obj)
seurat_obj.sce

logcounts(ref_rds) <- assay(ref_rds, "X")


# Cellannotation
pred <- SingleR(
  test = seurat_obj.sce,
  ref = ref_rds,
  assay.type.test = "counts",  # von brain.sce
  assay.type.ref  = "X",          # von ref_rds
  labels = ref_rds$cell_type     # Zelltyp-Spalte
)

# SingleR
jpeg("./plots_annotation/pred.jpeg", width = 2000, height = 1000, res = 200)
plotScoreHeatmap(pred)
print(pred)
dev.off()

# Add annotation to Seurat object
seurat_obj[["celltype"]] <- pred$labels

## Annotation specific plots
jpeg("./plots_annotation/dim_plot_by_celltype.jpeg", width = 2000, height = 1000, res = 200)
dim_plot_by_celltype <- DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", label = FALSE) +
ggtitle("DimPlot of celltypes")
print(dim_plot_by_celltype)
dev.off()

# SpatialDimPlot: Clustering overlaid on the image (label!)
jpeg("./plots_annotation/spatialdim_by_celltype.jpeg", width = 2000, height = 1000, res = 200)
spatialdim_by_celltype <- SpatialDimPlot(seurat_obj, label = FALSE, group.by = "celltype")+
  ggtitle("Location of celltypes")
print(spatialdim_by_celltype)
dev.off()

# Extract top cell types
main_cells <- sort(table(pred$labels), decreasing = TRUE)
top_celltypes <- names(head(main_cells, 10))

# Save top cell types
write.csv(as.data.frame(main_cells), "./csv/main_cells.csv", row.names = TRUE)

# Extra csv for loupe, optional
predforloupe <- subset(pred, select = "labels")
write.table(predforloupe, file="./csv/predforloupe.csv", sep=",")
