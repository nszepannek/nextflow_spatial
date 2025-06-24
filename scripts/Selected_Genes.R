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
seurat_path <- args[1]
gene_input <- args[2]

seurat_obj <- readRDS(seurat_path)

# Gene-IDs von Kommandozeile einlesen (gene_input sind ENSEMBL-IDs als Komma-getrennter String)
gene_ids <- if (gene_input != "") unlist(strsplit(gene_input, ",")) else character(0)

if (length(gene_ids) == 0) {
  message("No genes provided; skipping plotting.")
  file.create("not_found_genes.txt")
  quit(save = "no")
}

# ENSEMBL-IDs → Gensymbole mappen
gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = gene_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Verfügbare Gene im Seurat-Objekt (Rownames = Gensymbole)
available_genes <- rownames(seurat_obj)

# Welche gemappten Symbole sind im Seurat-Objekt vorhanden?
genes_found <- gene_symbols[!is.na(gene_symbols) & gene_symbols %in% available_genes]

# Nicht gefundene Gene (also entweder Mapping fehlgeschlagen oder nicht im Seurat-Objekt)
genes_not_found <- gene_ids[is.na(gene_symbols) | !(gene_symbols %in% available_genes)]

# Speichern der nicht gefundenen Gene
writeLines(genes_not_found, "not_found_genes.txt")

# Abbruch, wenn keine gültigen Gene gefunden wurden
if (length(genes_found) == 0) {
  message("None of the provided genes were found in the dataset.")
  quit(save = "no")
}

# Plots erzeugen
for (gene_symbol in genes_found) {
  p1 <- VlnPlot(seurat_obj, features = gene_symbol) + ggtitle(paste("Violin Plot:", gene_symbol))
  p2 <- SpatialFeaturePlot(seurat_obj, features = gene_symbol) + ggtitle(paste("Spatial:", gene_symbol))
  
  ggsave(paste0("violin_", gene_symbol, ".pdf"), p1, width = 6, height = 4)
  ggsave(paste0("spatial_", gene_symbol, ".pdf"), p2, width = 6, height = 5)
}

