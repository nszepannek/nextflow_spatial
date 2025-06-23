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
seurat_obj <- args[1]
gene_input <- args[2]

seurat_obj <- readRDS(seurat_path)

# Gene-IDs aus Kommandozeile parsen
gene_ids <- if (gene_input != "") unlist(strsplit(gene_input, ",")) else character(0)

if (length(gene_ids) == 0) {
  message("No genes provided; skipping plotting.")
  file.create("not_found_genes.txt")
  quit(save = "no")
}

# Prüfen, welche Gene im Datensatz vorhanden sind
available_genes <- rownames(seurat_obj)
genes_found <- gene_ids[gene_ids %in% available_genes]
genes_not_found <- setdiff(gene_ids, genes_found)

# Nicht gefundene Gene speichern
writeLines(genes_not_found, "not_found_genes.txt")

if (length(genes_found) == 0) {
  message("None of the provided genes were found in the dataset.")
  quit(save = "no")
}

# Plots erzeugen
for (gene in genes_found) {
  p1 <- VlnPlot(seurat_obj, features = gene) + ggtitle(paste("Violin Plot:", gene))
  p2 <- SpatialFeaturePlot(seurat_obj, features = gene) + ggtitle(paste("Spatial:", gene))
  
  ggsave(paste0("violin_", gene, ".pdf"), p1, width = 6, height = 4)
  ggsave(paste0("spatial_", gene, ".pdf"), p2, width = 6, height = 5)
}
