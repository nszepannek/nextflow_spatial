#!/usr/bin/env Rscript

# Load packages:
library(Seurat)
library(zellkonverter)

# https://datasets.cellxgene.cziscience.com/f87d516e-83fa-4ca4-a37a-ed1ab7ff2199.h5ad

ref_h5ad <- readH5AD("f87d516e-83fa-4ca4-a37a-ed1ab7ff2199.h5ad")

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