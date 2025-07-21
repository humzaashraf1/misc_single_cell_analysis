# Load required libraries
lapply(c("dplyr", "Seurat", "HGNChelper"), library, character.only = TRUE)

# Load ScType wrapper
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R")

# === Step 1: Load your Seurat object ===
# Replace with your actual .RDS file
seurat_obj <- readRDS("/path/to/object.seurat.RDS")

# Check for 'seurat_clusters'
if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  stop("The Seurat object does not contain 'seurat_clusters'.")
}

# === Step 2: Run ScType with full marker database ===
seurat_obj <- run_sctype(
  seurat_obj,
  assay = "RNA",
  scaled = TRUE,
  known_tissue_type = "Brain",
  custom_marker_file = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",
  name = "sctype_classification",
  plot = TRUE
)

# === Print cluster to cell type mapping ===
cat("Cluster to Cell Type Mapping:\n")
cluster_mapping <- seurat_obj@meta.data %>%
  dplyr::select(seurat_clusters, sctype_classification) %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(Major_Cell_Type = names(sort(table(sctype_classification), decreasing = TRUE)[1]))  # Most common label per cluster

print(cluster_mapping)

write.csv(cluster_mapping, "cluster_to_celltype_mapping.csv", row.names = FALSE)

# === Step 3: Save results ===
saveRDS(seurat_obj, file = "seurat_sctype_annotated_full.rds")
write.csv(seurat_obj@meta.data, "cell_metadata_with_sctype_full.csv")

