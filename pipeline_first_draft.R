library(Seurat)
library(future)
library(ggplot2)
library(arrow)
library(glmGamPoi)
library(patchwork)
library(RColorBrewer)
library(patchwork)

plan("multisession", workers = 25)
options(future.globals.maxSize = 25 * 1024^3) 

setwd("/home/fenosoa/projects/def-salehlab-ab/fenosoa/code_source_for_Pipeline/")

path <- "/home/fenosoa/projects/def-salehlab-ab/TMA_Xenium/output-XETG00325__0051618__TMA__20250124__223259" 


# Load the Xenium data
xenium.obj <- LoadXenium(path, fov = "fov")

# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

print("==================================================================================")
genes <- rownames(xenium.obj)

print("Save Genes as a text file")
# Save as a text file
write.table(genes, file = "genes_TMA.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

print("==================================================================================")
plot0 <- ImageDimPlot(xenium.obj, fov = "fov", molecules = c("RETREG1",   "RETREG3", "PIMREG", "FOXRED1"), nmols = 20000)
ggsave("/home/fenosoa/projects/def-salehlab-ab/fenosoa/code_source_for_Pipeline/ImageDimPlot_TMA_RETREG1.png", plot = plot0, width = 20, height = 15, dpi = 300)
print("==================================================================================")
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
#> dim(xenium.obj@assays$SCT@counts)
#[1]   248 36553

xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.3)



print("==================================================================================")
print("Save the xenium.obj with UMAP and SCT slot to disk under the name xenium.obj....And we can load it back")
save(xenium.obj, file = "/home/fenosoa/projects/def-salehlab-ab/fenosoa/code_source_for_Pipeline/data_object/UMAP_xenium.obj.RData")

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


# To run FindMarkers we need presto
# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')
# Let's find all the Markers of the 12 clusters.


# Create the "markers" folder if it does not exist
if (!dir.exists("markers")) {
  dir.create("markers")
}

# Initialize an empty list to store marker genes
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
  # write.csv(markers, file = filename, row.names = TRUE)
}

# We can zoom into a region of tissue, creating a new field of view. 
# For example, we can zoom into a region that contains the hippocampus. 
# Once zoomed-in, we can set DefaultBoundary() to show cell segmentations. 
# You can also ‘simplify’ the cell segmentations, 
# reducing the number of edges in each polygon to speed up plotting.



# Now, markers.0 to markers.11 are available as variables
# You can also access them using markers_list[["0"]], markers_list[["1"]], etc.



# Create the "markers" folder if it does not exist
if (!dir.exists("Feature")) {
  dir.create("Feature")
}


# c("ATP1A1", "ATP7B", "CD5", "CD6", "CD7", "CD8A")
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





# Problem with Cropping after doin the following command install.packages("sf", dependencies = TRUE)
# DONE (spatstat)
# ERROR: dependencies ‘sf’, ‘units’ are not available for package ‘stars’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/stars’
# ERROR: dependencies ‘leaflet’, ‘raster’, ‘sf’ are not available for package ‘leafem’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/leafem’
# ERROR: dependencies ‘leaflet’, ‘sf’ are not available for package ‘leafgl’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/leafgl’
# ERROR: dependency ‘leaflet’ is not available for package ‘leaflegend’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/leaflegend’
# ERROR: dependency ‘leaflet’ is not available for package ‘leafsync’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/leafsync’
# ERROR: dependencies ‘sf’, ‘lwgeom’, ‘stars’, ‘units’ are not available for package ‘tmaptools’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/tmaptools’
# ERROR: dependencies ‘leafem’, ‘leaflet’, ‘leafpop’, ‘raster’, ‘satellite’, ‘sf’ are not available for package ‘mapview’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/mapview’
# ERROR: dependencies ‘leafem’, ‘leafgl’, ‘leaflegend’, ‘leaflet’, ‘leafsync’, ‘sf’, ‘stars’, ‘tmaptools’, ‘units’ are not available for package ‘tmap’
# * removing ‘/home/fenosoa/R/x86_64-pc-linux-gnu-library/4.4/tmap’

# The downloaded source packages are in
# 	‘/tmp/RtmpDLzjBp/downloaded_packages’
# There were 18 warnings (use warnings() to see them)


# To run crop we need The package "sf" is required to overlay spatial information on the image.
# create a Crop
# cropped.coords <- Crop(xenium.obj[["fov"]], x = c(1750, 3000), y = c(3750, 5250), coords = "plot")
# # set a new field of view (fov)
# xenium.obj[["tumor"]] <- cropped.coords

# List of genes you want to plot


# If needed, set it to the correct assay, for example:
# DefaultAssay(xenium.obj) <- "Xenium"  

# genes_to_plot <- c("KHK", "CD4", "CYBA", "TMPRSS6")

# # Get the feature names from the default assay
# assay_features <- rownames(xenium.obj@assays[[DefaultAssay(xenium.obj)]])

# # Check which genes are present
# present_genes <- genes_to_plot[genes_to_plot %in% assay_features]
# missing_genes <- genes_to_plot[!genes_to_plot %in% assay_features]

# print(paste("Present genes:", paste(present_genes, collapse = ", ")))
# print(paste("Missing genes:", paste(missing_genes, collapse = ", ")))

# xenium.obj_up <- UpdateSeuratObject(xenium.obj)
# plot2 <- Seurat::SpatialPlot(xenium.obj_up, features = "KHK", pt.size = 0.1)
# ggsave("SpatialPlot.png", plot = plot2, width = 20, height = 15, dpi = 300)

# View(Seurat::SpatialPlot)
