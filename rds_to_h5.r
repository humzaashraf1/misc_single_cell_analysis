#!/usr/bin/env Rscript

# -------------------------
# Setup and Package Loading
# -------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(remotes)
  library(SeuratDisk)
})

cat("✅ Libraries loaded.\n")

# -------------------------
# Input RDS File (you edit this)
# -------------------------
rds_path <- "/path/to/object.seurat.RDS"
base <- sub("\\.rds$", "", basename(rds_path))
h5s_path <- paste0(base, ".h5Seurat")
h5ad_path <- paste0(base, ".h5ad")

cat("📦 Converting:\n")
cat("RDS       → ", rds_path, "\n")
cat("H5Seurat  → ", h5s_path, "\n")
cat("H5AD      → ", h5ad_path, "\n")

# -------------------------
# Load RDS
# -------------------------
obj <- readRDS(rds_path)

# Seurat v5 fix: ensure RNA assay is correct class
if ("RNA" %in% names(obj@assays)) {
  obj[["RNA"]] <- as(obj[["RNA"]], "Assay")
}
obj@meta.data[] <- lapply(obj@meta.data, function(x) if (is.factor(x)) as.character(x) else x)

# -------------------------
# Convert to H5Seurat and then H5AD
# -------------------------
SaveH5Seurat(obj, filename = h5s_path, overwrite = TRUE)
Convert(h5s_path, dest = "h5ad", overwrite = TRUE)

cat("✅ Conversion complete: ", h5ad_path, "\n")
