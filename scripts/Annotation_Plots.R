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


dir.create("plots_annotation", showWarnings = FALSE)
if (!dir.exists("csv")) dir.create("csv")


# Read in seurat-object with UMAP
seurat_obj <- readRDS("seurat_obj_with_umap.rds")

seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE)



#### Convert seurat object to single cell experiment
class(as.SingleCellExperiment(seurat_obj))
seurat_obj.sce <- as.SingleCellExperiment(seurat_obj)
seurat_obj.sce

# Read in reference
ref_h5ad <- readRDS("ref_h5ad.rds")

# ref_h5ad <- readH5AD("f87d516e-83fa-4ca4-a37a-ed1ab7ff2199.h5ad")

# Adapt Ensemble ID and Gene ID
symbols <- mapIds(
  org.Mm.eg.db,
  keys = rownames(ref_h5ad),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Umbenennen
rownames(ref_h5ad) <- symbols
saveRDS(ref_h5ad, file = "ref_h5ad.rds")

logcounts(ref_h5ad) <- assay(ref_h5ad, "X")


# Cellannotation
pred <- SingleR(
  test = seurat_obj.sce,
  ref = ref_h5ad,
  assay.type.test = "counts",  # von brain.sce
  assay.type.ref  = "X",          # von ref_h5ad
  labels = ref_h5ad$cell_type     # Zelltyp-Spalte
)

# SingleR aufrufen
#pred <- SingleR(test = seurat_obj.sce, ref = ref_h5ad, labels = ref$label.main)
# pred <- SingleR(test = seurat_obj.sce, ref = ref_h5ad, labels = ref_h5ad$cell_type)
# Ergebnis z.B. als DataFrame
# Ergebnis z.B. als DataFrame
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

# Top Zelltypen extrahieren (nach Häufigkeit)
main_cells <- sort(table(pred$labels), decreasing = TRUE)
top_celltypes <- names(head(main_cells, 10))

# Speichere Zelltypen-Häufigkeit
write.csv(as.data.frame(main_cells), "./csv/main_cells.csv", row.names = TRUE)

predforloupe <- subset(pred, select = "labels")
write.table(predforloupe, file="./csv/predforloupe.csv", sep=",")


