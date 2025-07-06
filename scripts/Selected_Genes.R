#!/usr/bin/env Rscript

# Load packages:
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Output of first R script (seurat_obj_with_umap.rds) is the input given in nextflow
args <- commandArgs(trailingOnly = TRUE)
seurat_path <- args[1]
genes <- args[2]

genes <- unlist(strsplit(genes, split = ","))

# Load seurat object
seurat_obj <- readRDS(seurat_path)

dir.create("plots_selected_genes", showWarnings = FALSE)

# Get all available features
all_genes <- rownames(seurat_obj)

# Loop through provided genes
for (gene in genes) {
  if (!(gene %in% all_genes)) {
    message(paste0("Gene not found in Seurat object: ", gene))
    next
  }

  # VlnPlot
  jpeg(paste0("./plots_selected_genes/VlnPlot_", gene, ".jpeg"), width = 1000, height = 1000, res = 200)
  p1 <- VlnPlot(seurat_obj, features = gene) +
    ggtitle(paste("VlnPlot –", gene))
  print(p1)
  dev.off()

  # RidgePlot
  jpeg(paste0("./plots_selected_genes/RidgePlot_", gene, ".jpeg"), width = 1000, height = 1000, res = 200)
  p2 <- RidgePlot(seurat_obj, features = gene) +
    ggtitle(paste("RidgePlot –", gene))
  print(p2)
  dev.off()

  # FeaturePlot
  jpeg(paste0("./plots_selected_genes/FeaturePlot_", gene, ".jpeg"), width = 1000, height = 1000, res = 200)
  p3 <- FeaturePlot(seurat_obj, features = gene) +
    ggtitle(paste("FeaturePlot –", gene))
  print(p3)
  dev.off()
}