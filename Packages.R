# Installiere BiocManager, falls noch nicht vorhanden
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# CRAN-Pakete installieren
cran_packages <- c(
  "Seurat", "ggplot2", "patchwork", "dplyr",
  "pheatmap", "viridis", "reticulate", "hdf5r", "Rfast2"
)
install.packages(cran_packages, repos = "https://cloud.r-project.org")

# Bioconductor-Pakete installieren
bioc_packages <- c(
  "zellkonverter", "SingleCellExperiment",
  "org.Mm.eg.db", "SingleR", "glmGamPoi"
)
BiocManager::install(bioc_packages, ask = FALSE)
