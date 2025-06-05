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

if (!dir.exists("img")) {
  dir.create("img")
}

# Create the "markers" folder if it does not exist
if (!dir.exists("markers")) {
  dir.create("markers")
}

# Initialize an empty list to store marker genes
markers_list <- list()

# Suppose 'n' is the number of idents (e.g., length of levels or some predefined integer)
# Make sure 'n' is defined before this loop:
# n <- length(levels(xenium.obj))  # or however you determine the number of idents

for (ident in 0:(n - 1)) {
  # Find marker genes for the current ident (as character)
  markers <- FindMarkers(xenium.obj, ident.1 = as.character(ident))
  
  # Store the result in a globally named variable
  assign(paste0("markers.", ident), markers, envir = .GlobalEnv)
  
  # Also store it in a list for easy future reference
  markers_list[[as.character(ident)]] <- markers
  
  # Define the filename dynamically
  filename <- paste0("markers/markers.", ident, ".csv")
  
  # Save the markers data as a CSV file
  write.csv(markers, file = filename, row.names = TRUE)
  message("Saved: ", filename)
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

# Create the "Feature" folder if it does not exist
if (!dir.exists("Feature")) {
  dir.create("Feature")
}

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


