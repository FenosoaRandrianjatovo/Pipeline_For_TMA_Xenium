library(Seurat)
library(future)
library(ggplot2)
library(arrow)
library(glmGamPoi)
library(patchwork)
library(RColorBrewer)
library(patchwork)

setwd("/home/fenosoa/projects/def-salehlab-ab/fenosoa/code_source_for_Pipeline/data_object/")

path <- "/home/fenosoa/projects/def-salehlab-ab/fenosoa/code_source_for_Pipeline/data_object/UMAP_v2_xenium.obj.RData"
temp_env <- new.env()
print("xenium.obj :   Loading begin")
load(path, envir = temp_env)

xenium.obj <- temp_env$xenium.obj


groups <- levels(xenium.obj)

print("Nomber of Groups")
print(length(groups))

n <- length(groups)


print("# 1. Retrieve the SCT data matrix")
counts_mat <- GetAssayData(xenium.obj, assay = "SCT", slot = "data")

print("# 2. Select all cell barcodes (no downsampling)")
selected_cells <- colnames(counts_mat)

print("# 3. Subset the matrix to all cells (genes stay the same)")
counts_sub <- counts_mat[, selected_cells]

print("# 4. Transpose & coerce to a dense data.frame")
counts_df <- as.data.frame(t(as.matrix(counts_sub)))

print("# 5. Add the corresponding cluster labels for every cell")
counts_df$cluster <- xenium.obj@meta.data[selected_cells, "seurat_clusters"]

print("# 6. Write out to CSV file")
write.csv(
  counts_df,
  file = "TMA_counts_with_clusters_SCT_data_AllCells.csv",
  row.names = TRUE
)
