library(Seurat)
library(future)
library(ggplot2)
library(arrow)
library(glmGamPoi)
library(patchwork)
library(RColorBrewer)
library(patchwork)


cl <- parallelly::makeClusterPSOCK(25, outfile="")
plan(cluster, workers = cl)
# options(future.globals.maxSize = 700 * 1024^3) 

# remotes::install_version(package = 'Seurat', version = package_version('5.2.0'))

setwd("/home/fenosoa/projects/def-salehlab-ab/fenosoa/code_source_for_Pipeline/")

print("==================================================================================")
print(sessionInfo())
print("==================================================================================")


# path <- "/home/fenosoa/projects/def-salehlab-ab/TMA_Xenium/output-XETG00325__0051618__TMA__20250124__223259" 


# print("# Load the Xenium data")
# xenium.obj <- LoadXenium(path, fov = "fov")

# print("# remove cells with 0 counts")
# xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
# print("==================================================================================")
# print("Save xenium.obj as a R Data file")
# save(xenium.obj, file = "/home/fenosoa/projects/def-salehlab-ab/fenosoa/code_source_for_Pipeline/data_object/positive_count_xenium.obj.RData")
# print("Save xenium.obj as a R Data file is done")

print("==================================================================================")

print("# Load the Xenium data from positive_count_xenium.obj.RData")

path <- "/home/fenosoa/projects/def-salehlab-ab/fenosoa/code_source_for_Pipeline/data_object/positive_count_xenium.obj.RData"
temp_env <- new.env()
print("xenium.obj :   Loading begin")
load(path, envir = temp_env)

xenium.obj <- temp_env$xenium.obj

print("==================================================================================")

genes <- rownames(xenium.obj)

print("Save Genes as a text file")
# Save as a text file
write.table(genes, file = "genes_TMA.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


print("==================================================================================")
# plot0 <- ImageDimPlot(xenium.obj, fov = "fov", molecules = c("RETREG1",   "RETREG3", "PIMREG", "FOXRED1"), nmols = 20000)
# ggsave("/home/fenosoa/projects/def-salehlab-ab/fenosoa/code_source_for_Pipeline/ImageDimPlot_TMA_RETREG1_FOXRED1_v1.png", plot = plot0, width = 20, height = 15, dpi = 300)
print("==================================================================================")

print("SCTransform is running")
# xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
options(future.globals.maxSize = 480 * 1024^3)  # Set the limit to 400 GB

xenium.obj <- SCTransform(
  xenium.obj,
  assay = "Xenium",
  layer = "counts"      
)


#> dim(xenium.obj@assays$SCT@counts)
#[1]   248 36553

print("RunPCA is running")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))

print("RunUMAP is running")
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)

print("FindNeighbors is running")
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)


print("FindClusters is running")
xenium.obj <- FindClusters(xenium.obj, resolution = 0.3)

print("SCTransform, FindClusters and RunUMAP are done")
print("==================================================================================")
print("Save the xenium.obj with UMAP and SCT slot to disk under the name xenium.obj....And we can load it back")
save(xenium.obj, file = "/home/fenosoa/projects/def-salehlab-ab/fenosoa/code_source_for_Pipeline/data_object/UMAP_v2_xenium.obj.RData")

print("Saving is DONE")
print("==================================================================================")


if (!dir.exists("img")) {
  dir.create("img")
}
plot4 <- DimPlot(xenium.obj)
ggsave("img/UMAP_Plot.png", plot = plot4, width = 20, height = 15, dpi = 300)


plot1 <- FeaturePlot(xenium.obj, features = c("TREM1", "CD163"))
ggsave("img/FeaturePlot_of_MregTREM1__CD163.png", plot = plot1, 
       width = 80, height = 50, dpi = 300, limitsize = FALSE)


plot2 <- ImageDimPlot(xenium.obj, cols = "polychrome", size = 0.75)
ggsave("img/ImageDimPlot_R1.png", plot = plot2, width = 20, height = 15, dpi = 300)

ggsave("img/ImageDimPlot_VS_UMAP_R1.png", plot = plot2+ plot4, width = 30, height = 20, dpi = 300)

p <- ImageDimPlot(xenium.obj, fov = "fov", cols = "orange", 
                    cells = WhichCells(xenium.obj, idents = "3"))

ggsave("img/ImageDimPlot_for_3_R1.png", plot = p, width = 30, height = 20, dpi = 300)


groups <- levels(xenium.obj)

print("Nomber of Groups")
print(length(groups))

n <- length(groups)
# We have 12 Cluster

print("==================================================================================")


if (!dir.exists("ident")) {
  dir.create("ident")
}

# Loop through ident values from 0 to 11
for (ident in 0:(n-1)) {
  # Generate the plot
  p <- ImageDimPlot(xenium.obj, fov = "fov", cols = "red", 
                    cells = WhichCells(xenium.obj, idents = ident))
  
  # Define the filename dynamically
  filename <- paste0("ident/ImageDimPlot_for_ident", ident, ".png")
  
  # Save the plot
  ggsave(filename, plot = p, width = 30, height = 20, dpi = 300)
}


p1 <- ImageFeaturePlot(xenium.obj, features =c("TREM1", "CD163"), mols.cols = c("red", "blue"))
p2 <- ImageDimPlot(xenium.obj, molecules = c("TREM1", "CD163"), nmols = 10000, alpha = 0.3, mols.cols = c("red", "blue"))

ggsave("img/ImageDimPlot_for_melecule_MREG_TREM1__CD163.png", plot = p1 + p2, width = 30, height = 20, dpi = 300)



p1 <- ImageFeaturePlot(xenium.obj, features =c("FoxP3", "CD3E"))
p2 <- ImageDimPlot(xenium.obj, molecules = c("FoxP3", "CD3E"), nmols = 10000, alpha = 0.3, mols.cols = c("red", "blue"))

ggsave("img/ImageDimPlot_for_melecule_TREG_FoxP3__CD3E.png", plot = p1 + p2, width = 30, height = 20, dpi = 300)



p1 <- ImageFeaturePlot(xenium.obj, features =c("EPCAM", "CD24", "CD47"))
ggsave("img/test.png", plot = p1 , width = 30, height = 20, dpi = 300)
p2 <- ImageDimPlot(xenium.obj, molecules = c("EPCAM", "CD24", "CD47"), nmols = 10000, alpha = 0.3, mols.cols = c("red", "blue", "green"))

ggsave("img/ImageDimPlot_for_melecule_Cancer_stem_cells__EPCAM__CD24__CD47.png", plot = p1 + p2, width = 30, height = 20, dpi = 300)


p1 <- ImageFeaturePlot(xenium.obj, features =c("CD19", "CD200", "CD22"))
ggsave("img/test.png", plot = p1 , width = 30, height = 20, dpi = 300)
p2 <- ImageDimPlot(xenium.obj, molecules = c("CD19", "CD20", "CD22"), nmols = 10000, alpha = 0.3, mols.cols = c("red", "blue", "green"))

ggsave("img/ImageDimPlot_for_melecule_Bcells__CD19__CD20__CD22.png", plot = p1 + p2, width = 30, height = 20, dpi = 300)


# Define the features of interest
features <- c("CD19", "CD200", "CD22","EPCAM", "CD47","FOXP3", "CD3E","TREM1", "CD163")

# features <- c("FOXP3", "CD3E","TREM1", "CD163")


# Loop over each feature and save the corresponding plots
for (feature in features) {
  # Generate FeaturePlot and ImageFeaturePlot for the feature
  feature_plot <- FeaturePlot(xenium.obj, features = feature)
  image_feature_plot <- ImageFeaturePlot(xenium.obj, features = feature)
  
  # Combine the plots
  combined_plot <- image_feature_plot + feature_plot
  
  # Save the plot
  file_name <- paste0("img/feature_plot_VS_ImageFeaturePlot_", feature, ".png")
  ggsave(file_name, plot = combined_plot, width = 30, height = 20, dpi = 300)
  
  print(paste("Saved:", file_name))
}

if (!dir.exists("VlnPLot")) {
  dir.create("VlnPLot")
}

for (feature in features) {
      # Generate violin plot
      te <- VlnPlot(xenium.obj, features = feature, pt.size = 0.1)

      # Define filename dynamically
      filename <- paste0("VlnPLot/VlnPLot_gene_", feature, ".png")
      
      # Save the plot
      ggsave(filename, plot = te, width = 20, height = 15, dpi = 300)
      print(paste("Saved:", file_name))

  }




# Create the "markers" folder
if (!dir.exists("markers")) {
  dir.create("markers")
}


markers_list <- list()

# Loop through ident values from 0 to 11
for (ident in 0:(n-1)) {
  # Find marker genes for the current ident
  markers <- FindMarkers(xenium.obj, ident.1 = as.character(ident))
  
  # Store the result in a dynamically named variable
  assign(paste0("markers.", ident), markers, envir = .GlobalEnv)
  
  # Also store it in a list for easy future reference
  markers_list[[as.character(ident)]] <- markers
  
  # Define the filename dynamically
  filename <- paste0("markers/markers.", ident, ".csv")
  
  # Save the markers data as a CSV file
  print(paste("Saved:", filename))
  # write.csv(markers, file = filename, row.names = TRUE)
}





# Create the "markers" folder if it does not exist
if (!dir.exists("Feature")) {
  dir.create("Feature")
}



pa <- FeaturePlot(xenium.obj, features = rownames(markers.0)[1:6])
pb <- ImageFeaturePlot(xenium.obj, features = rownames(markers.0)[1:6])
ggsave("Feature/FeaturePlot_ImageFeaturePlot_6_markers.0.png", plot = pb, width = 50, height = 40, dpi = 300)



# Define marker indices (0 to 11)
marker_indices <- 0:(n-1)

# Loop through marker objects from markers.0 to markers.11
for (i in 0:(n-1)) {
  # Construct the variable name dynamically
  marker_var <- get(paste0("markers.", i), envir = .GlobalEnv)
  
  # Ensure the object exists and has at least 6 rows
  if (!is.null(marker_var) && nrow(marker_var) >= 6) {
    
    # Loop through the first six selected genes
    for (j in 1:6) {
      selected_gene <- rownames(marker_var)[j]
      
      # Generate FeaturePlot and ImageFeaturePlot for the current gene
      feature_plot <- FeaturePlot(xenium.obj, features = selected_gene)
      image_feature_plot <- ImageFeaturePlot(xenium.obj, features = selected_gene)
      
      # Define the filename dynamically
      filename <- paste0("Feature/FeaturePlot_ImageFeaturePlot_", selected_gene, "_ident", i, ".png")
      
      # Save the ImageFeaturePlot
      ggsave(filename, plot = image_feature_plot + feature_plot, width = 50, height = 40, dpi = 300, limitsize = FALSE)
      
      print(paste("Saved plot:", filename))
    }
  } else {
    print(paste("Skipping markers.", i, " - Not enough genes or NULL"))
  }
}



# Plot a marker gene’s expression: VlnPlot
if (!dir.exists("VlnPLot")) {
  dir.create("VlnPLot")
}
# Loop through markers.0 to markers.11
for (i in 0:(n-1)) {
  marker_var <- get(paste0("markers.", i))  # Get marker dataset dynamically
  
  if (!is.null(marker_var) && nrow(marker_var) >= 6) {
    for (j in 1:6) {  # Loop through the first 6 genes
      gene_name <- rownames(marker_var)[j]

      # Generate violin plot
      te <- VlnPlot(xenium.obj, features = gene_name, pt.size = 0.1)

      # Define filename dynamically
      filename <- paste0("VlnPLot/markers.", i, "_gene_", gene_name, "_VlnPlot.png")
      
      # Save the plot
      ggsave(filename, plot = te, width = 20, height = 15, dpi = 300)

      print(paste("Saved plot:", filename))
    }
  } else {
    print(paste("Skipping markers.", i, " - Not enough genes or NULL"))
  }
}

if (!dir.exists("VlnPlot_Feature_Count")) {
  dir.create("VlnPlot_Feature_Count")
}

test <- VlnPlot(xenium.obj, features = "nFeature_Xenium", ncol = 1, pt.size = 0.1)
ggsave("VlnPlot_Feature_Count/VlnPlot_nFeature_Xenium_after_Clustering.png", plot = test, width = 20, height = 15, dpi = 300)

test1 <- VlnPlot(xenium.obj, features = "nCount_Xenium", ncol = 1, pt.size = 0.1)
ggsave("VlnPlot_Feature_Count/VlnPlot_nCount_Xenium_after_Clustering.png", plot = test1, width = 20, height = 15, dpi = 300)


before <- VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), group.by = "orig.ident")
ggsave("VlnPlot_Feature_Count/VlnPlot_nFeature_Xenium_nCount_Xenium_before_Clustering.png", plot = before, width = 20, height = 15, dpi = 300)



# Run FindAllMarkers once to get discriminating genes for every cluster
all_markers <- FindAllMarkers(
  xenium.obj,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

# Split by cluster and write out one CSV per cluster
markers_by_cluster <- split(all_markers, all_markers$cluster)

for (cluster in names(markers_by_cluster)) {
  filename <- paste0("markers/markers.", cluster, ".csv")
  write.csv(markers_by_cluster[[cluster]], file = filename, row.names = TRUE)
  message("Saved: ", filename)
}



