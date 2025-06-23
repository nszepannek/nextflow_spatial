#!/usr/bin/env Rscript

# Load packages:
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Output of first R script (seurat_obj_with_umap.rds) is the input given in nextflow
args <- commandArgs(trailingOnly = TRUE)
seurat_path <- args[1]

# Load seurat object
seurat_obj <- readRDS(seurat_path)

# Create directories for results (plots and csv)
dir.create("plots", showWarnings = FALSE)
dir.create("csv", showWarnings = FALSE)

# DimPlot: clustering in UMAP space
jpeg("./plots/dim_plot_by_umap.jpeg", width = 2000, height = 1000, res = 200)
dim_plot_by_umap <- DimPlot(seurat_obj, reduction = "umap", label = FALSE) +
  ggtitle("DimPlot of clusters")
print(dim_plot_by_umap)
dev.off()

# SpatialDimPlot: Clustering overlaid on the image (label!)
jpeg("./plots/spatialdim_by_umap.jpeg", width = 2000, height = 1000, res = 200)
spatialdim_by_umap <- SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3)+
  ggtitle("Location of clusters")
print(spatialdim_by_umap)
dev.off()

# Distinguish location of top 6 clusters (takes a few seconds)
jpeg("./plots/spatial_dim_top_6_cluster.jpeg", width = 2000, height = 1000, res = 200)
spatial_dim_top_6_cluster <- SpatialDimPlot(seurat_obj, cells.highlight = CellsByIdentities(object = seurat_obj, idents = c(0, 1, 2, 3,
                                                                                                                           4, 5)), facet.highlight = TRUE, ncol = 3) +
  ggtitle("Location of top 6 clusters")
print(spatial_dim_top_6_cluster)
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
Spatial_variable_features_top6 <- SpatialFeaturePlot(seurat_obj, features = top.features, ncol = 3, alpha = c(0.1, 1))+
  ggtitle("Location of top 6 variable features")
print(Spatial_variable_features_top6)
dev.off()

# Schleife über alle Gene
for (gene in top.features) {
  
  # --- VlnPlot ---
  jpeg(paste0("./plots/VlnPlot_", gene, ".jpeg"), width = 1000, height = 1000, res = 200)
  p1 <- VlnPlot(seurat_obj, features = gene) +
    ggtitle(paste("VlnPlot –", gene))
  print(p1)
  dev.off()
  
  # --- RidgePlot ---
  jpeg(paste0("./plots/RidgePlot_", gene, ".jpeg"), width = 1000, height = 1000, res = 200)
  p2 <- RidgePlot(seurat_obj, features = gene) +
    ggtitle(paste("RidgePlot –", gene))
  print(p2)
  dev.off()
  
  # --- FeaturePlot ---
  jpeg(paste0("./plots/FeaturePlot_", gene, ".jpeg"), width = 1000, height = 1000, res = 200)
  p3 <- FeaturePlot(seurat_obj, features = gene) +
    ggtitle(paste("FeaturePlot –", gene))
  print(p3)
  dev.off()
}

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

top_50genes <- head(VariableFeatures(seurat_obj), 50)

jpeg("./plots/DoHeatmap.jpeg", width = 1500, height = 1000, res = 200)
print(
  DoHeatmap(seurat_obj, features = top_50genes, size = 4, angle = 90) + NoLegend()
)
dev.off()

# Aktiviere Cluster als Identität
Idents(seurat_obj) <- "seurat_clusters"  # oder "celltype", falls annotiert

# Alle Cluster durchgehen
cluster_ids <- names(sort(table(Idents(seurat_obj)), decreasing = TRUE))
for (cluster in cluster_ids) {
  # Marker berechnen
  markers <- FindMarkers(seurat_obj, ident.1 = cluster)
  
  # Rownames als Spalte speichern
  markers$gene <- rownames(markers)
  
  # Top 10 Gene fürs Label
  top10 <- markers %>% arrange(p_val_adj) %>% head(10)
  
  # Volcano Plot
  p <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_text(data = top10, aes(label = gene), vjust = 1.5, size = 3) +
    labs(
      title = paste("Volcano Plot – Cluster", cluster),
      x = "Log2 Fold Change",
      y = "-log10 adjusted p-value"
    ) +
    theme_classic(base_size = 14)
  
  # Plot speichern
  plot_filename <- paste0("./plots/volcano_cluster_", cluster, ".jpg")
  ggsave(plot_filename, plot = p, width = 8, height = 6, dpi = 300)
  
  # CSV speichern
  csv_filename <- paste0("./csv/DEG_cluster_", cluster, ".csv")
  write.csv(markers, file = csv_filename, row.names = FALSE)
}


# Additional: interactive plotting, subset of anatomical regions, wokrking with mulitple slices
# --> Bioinformatik/Masterarbeit/Seurat/Seurat.R





