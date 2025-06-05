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



