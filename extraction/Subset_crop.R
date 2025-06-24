library(Seurat)
library(future)
library(ggplot2)
library(arrow)
library(glmGamPoi)
library(patchwork)
library(RColorBrewer)
library(patchwork)
plan(sequential)


setwd("/Users/fenosoa/Desktop/TMA_test/")

path <- "/Users/fenosoa/Desktop/TMA_test/positive_count_xenium.obj.RData"

temp_env <- new.env()
print("xenium.obj :   Loading begin")
load(path, envir = temp_env)

xenium.obj <- temp_env$xenium.obj



centroids <- GetTissueCoordinates(xenium.obj[["fov"]], coords = "centroids")


img_name <- names(xenium.obj@images)[1]

test <- xenium.obj@images[[img_name]]

coords_df <- xenium.obj@images[[img_name]]@molecules$molecules@.Data

meta <- xenium.obj@meta.data

meta2 <- cbind(meta, coords_df[rownames(meta), ])


# after extracting or adding x/y into meta.data
summary(xenium.obj$x)
summary(xenium.obj$y)

centroids <- GetTissueCoordinates(xenium.obj[["fov"]], coords = "tissue")


     <- plot(centroids$x, centroids$y,
     pch = 20, cex = 0.5,
     xlab = "x", ylab = "y",
     main = "All cell centroids")




options(future.globals.maxSize = 400 * 1024^3)

xrange <- sort(c(1621.154,  163.7940))
yrange <- sort(c(2305.799, 1198.2369))

cropped.coords <- Crop(
  xenium.obj[["fov"]],
  x      = xrange,
  y      = yrange,
  coords = "tissue"
)

#cropped.coords <- Crop(xenium.obj[["fov"]], x = c(1750, 3000), y = c(3750, 5250), coords = "plot")


cropped.coords <- Crop(xenium.obj[["fov"]], x = c(50, 2500), y = c(45, 4000), coords = "plot")
#FOV1 [526-C]
xenium.obj[["fov1"]] <- cropped.coords


p1 <- ImageDimPlot(xenium.obj, fov = "fov1", axes = TRUE, size = 0.7, border.color = "white", cols = "polychrome",
                   coord.fixed = FALSE)



p1

hippo_centroids <- GetTissueCoordinates(
  xenium.obj[["fov1"]],
  coords = "tissue"
)

hippo_cells <- hippo_centroids$cell

hippo_seurat <- subset(
  x = xenium.obj,
  cells = hippo_cells
)


hippo_seurat <- SCTransform(hippo_seurat, assay = "Xenium")
hippo_seurat <- RunPCA(   hippo_seurat, npcs = 30, features = rownames(hippo_seurat) )
hippo_seurat <- RunUMAP(  hippo_seurat, dims = 1:30  )
hippo_seurat <- FindNeighbors(hippo_seurat, reduction = "pca", dims = 1:30)
hippo_seurat <- FindClusters( hippo_seurat, resolution = 0.3 )

q1 <- DimPlot(hippo_seurat)

q2<- ImageDimPlot(hippo_seurat,alpha = 1, cols = "polychrome")

p1 <- ImageFeaturePlot(hippo_seurat, features =c("FoxP3", "CD3E"),cols = "polychrome")
p2 <- ImageDimPlot(hippo_seurat, features = c("FoxP3", "CD3E"), cols = "polychrome", nmols = 10000, alpha = 0.3, mols.cols = c("red", "blue"))

ggsave("TEST_UMAP_for_FOV1.png", plot = q1, bg='white')



















