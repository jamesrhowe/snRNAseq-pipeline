## ----setup---------------------------------------------------------------------------------------
require(Seurat)
require(plotly)
require(tidyverse)

source("../misc.R")


## ----1-merge_specified_libs, warning = FALSE, message = FALSE------------------------------------
knitr::purl("MergeNormalize.Rmd")

args <- commandArgs(TRUE)
path <- args[1]
id <- args[2]

file_names <- list.files(path = "../../datasets/preprocessed/", pattern = "*.rds", full.names = TRUE)
for (i in 1:length(file_names)){
  array[[i]] <- readRDS(file_names[i])
  
}
array <- Reduce(function(...) merge(..., all=T), array)

paste("Total", id, "cells:", length(colnames(array)))
count_all(paste(id, "cell"))

plot_array <- cbind.data.frame(array$nCount_RNA, 
                               array$nFeature_RNA,
                               array$orig.ident,
                               colnames(array))
colnames(plot_array) <- c("UMIs", "Features", "Batch", "Cell_ID")

plot_ly(data = plot_array, x = ~UMIs, y = ~Features, color = ~Batch, text = ~Cell_ID) %>%
  add_markers(marker = list(size = 2, sizemode = "area")) %>%
  layout(title = paste("Batch Depth:", id))


## ----2-gene_filter, warning = FALSE, message = FALSE---------------------------------------------
plot_array <- cbind.data.frame(Matrix::rowMeans(array@assays$RNA@counts),
                               Matrix::rowSums(array@assays$RNA@counts != 0),
                               ifelse(Matrix::rowSums(array@assays$RNA@counts != 0) > 5, 
                                      "<5 Cells Expressing", ">5 Cells Expressing"),
                               rownames(array))
colnames(plot_array) <- c("Mean", "Total_cells", "Filter", "Gene_ID")

array <- subset(array, features = which(Matrix::rowSums(array) > 5))

count_all(paste(id, "cells after gene filter"))

plot_ly(data = plot_array) %>%
  add_markers(x = ~Mean, y = ~Total_cells, color = ~Filter, text = ~Gene_ID, 
              marker = list(size = 2, sizemode = "area")) %>%
  layout(title = paste("5 Cell Expression Filter:", id),
         xaxis = list(type = "log"), yaxis = list(type = "log"))


## ----3-sctransform, warning = FALSE, message = TRUE----------------------------------------------
# needs conserve.memory = TRUE, or else it crashes R session on laptop
array <- SCTransform(array, variable.features.n = 5000, 
                     conserve.memory = TRUE, verbose = FALSE)

plot_array <- cbind.data.frame(array@assays$SCT@meta.features$sct.gmean,
                               array@assays$SCT@meta.features$sct.residual_variance,
                               ifelse(array@assays$SCT@meta.features$sct.variable, "Variable", "Non-variable"),
                               rownames(array))
colnames(plot_array) <- c("Geometric_Mean", "Residual_Variance", "Variable_Feature", "Gene_ID")

plot_ly(data = plot_array) %>%
  add_markers(x = ~Geometric_Mean, y = ~Residual_Variance, color = ~Variable_Feature, text = ~Gene_ID, 
              marker = list(size = 2, sizemode = "area")) %>%
  layout(title = paste("Variable Feature Plot:", id), xaxis = list(type = "log"))


## ----4-pca, warning = FALSE, message = TRUE------------------------------------------------------
array <- RunPCA(array, verbose = FALSE)

plot_array <- cbind.data.frame(array@reductions$pca@cell.embeddings[,1:3], 
                               array$orig.ident,
                               colnames(array))
colnames(plot_array) <- c("PC1", "PC2", "PC3", "Batch", "Cell_ID")

# visualize PC weights and gene contributions to PCs
ElbowPlot(array)
VizDimLoadings(array, dims = 1, reduction = "pca")
VizDimLoadings(array, dims = 10, reduction = "pca")
VizDimLoadings(array, dims = 20, reduction = "pca")
VizDimLoadings(array, dims = 30, reduction = "pca")
VizDimLoadings(array, dims = 40, reduction = "pca")
VizDimLoadings(array, dims = 50, reduction = "pca")

# visualize 2D and 3D PCA
plot_ly(data = plot_array, x = ~PC1, y = ~PC2, color = ~Batch, text = ~Cell_ID) %>%
  add_markers(marker = list(size = 2, sizemode = "area")) %>%
  layout(title = paste("2D PCA:", id))

plot_ly(data = plot_array, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Batch, text = ~Cell_ID) %>%
  add_markers(marker = list(size = 2, sizemode = "area")) %>%
  layout(title = paste("3D PCA:", id))


## ----5-umap, message = FALSE, warning = FALSE----------------------------------------------------
array <- RunUMAP(array, dims = 1:50, n.epochs = 1000, verbose = FALSE)
array <- RunUMAP(array, dims = 1:50, n.epochs = 1000, verbose = FALSE,
                 n.components = 3, reduction.name = "umap3d", reduction.key = "UMAP3D_",)

# visualize 2D UMAP for batch
plot_array_2d <- cbind.data.frame(as.data.frame(array[["umap"]]@cell.embeddings),
                                  array$orig.ident,
                                  colnames(array))
colnames(plot_array_2d) <- c("UMAP1", "UMAP2", "Batch", "Cell_ID")
plot_ly(data = plot_array_2d, x = ~UMAP1, y = ~UMAP2, color = ~Batch, text = ~Cell_ID) %>%
  add_markers(marker = list(size = 2, sizemode = "area")) %>%
  layout(title = paste("2D UMAP:", id))

# visualize 3D UMAP for batch
plot_array_3d <- cbind.data.frame(as.data.frame(array[["umap3d"]]@cell.embeddings),
                                  array$orig.ident,
                                  colnames(array))
colnames(plot_array_3d) <- c("UMAP1", "UMAP2", "UMAP3", "Batch", "Cell_ID")
plot_ly(data = plot_array_3d, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3, color = ~Batch, text = ~Cell_ID) %>%
  add_markers(marker = list(size = 2, sizemode = "area")) %>%
  layout(title = paste("3D UMAP:", id))


## ----8-save_output, warning = FALSE, message = FALSE---------------------------------------------
saveRDS(array, file = paste0("../../datasets/merged/", id, "_merged.rds"))

