## ----setup---------------------------------------------------------------------------------------
require(Seurat)
require(plotly)
require(DropletUtils)
require(scDblFinder)
require(tidyverse)

source("analysis/misc.R")


## ----1-read_raw_mtx, message = FALSE, warning = FALSE--------------------------------------------

args <- commandArgs(TRUE)
path <- args[1]
id <- args[2]
id_general <- substr(id, 0, 3)
print(path)
print(id)

array <- Read10X(data.dir = path, strip.suffix = TRUE)
colnames(array) <- paste0(id, "_", colnames(array))


## ----2-cell_barcode_call, message = FALSE, warning = FALSE---------------------------------------
# run emptyDrops
barcode_filter <- emptyDrops(array, retain = 1000)
barcode_filter$FDR[is.na(barcode_filter$FDR)] <- 1
knee_ranks <- barcodeRanks(array)

# export metrics
plot_array <- list(cbind.data.frame(knee_ranks$rank, 
                                 knee_ranks$total, 
                                 barcode_filter$FDR),
                     sum(barcode_filter$FDR < 0.001),
                     (sum(array[,barcode_filter$FDR < 0.001]) / sum(array)),
                     knee_ranks@metadata$knee,
                     knee_ranks@metadata$inflection)
names(plot_array) <- c("Metrics", "Total", "UMIs_in_cells", "Knee_threshold", "Inflection_threshold")
colnames(plot_array$Metrics) <- c("Rank", "UMIs", "FDR")
rownames(plot_array$Metrics) <- c(knee_ranks@rownames)

# actually filter the array
array <- array[, plot_array$Metrics$FDR < 0.001]

# give first stage QC metrics
cat(paste("Called barcodes:", length(colnames(array))), 
    paste("Mean UMIs/called barcode:", sum(array)/length(colnames(array))),
    paste("Mean Features/called barcode:", sum(array != 0)/length(colnames(array))),
    paste("Percent reads in called barcodes:", plot_array$UMIs_in_cells), 
    paste("Knee point:", plot_array$Knee_threshold),
    paste("Inflection point:", plot_array$Inflection_threshold), sep = "\n")

# remove duplicates, will crash notebook and computer if retained in plot
plot_array$Metrics <- plot_array$Metrics[!duplicated(plot_array$Metrics$Rank),]
plot_array$Metrics <- plot_array$Metrics[order(plot_array$Metrics$Rank),]

# produce knee plot for output of emptyDrops
plot_ly(data = plot_array$Metrics, x = ~Rank, y = ~UMIs) %>%
  add_markers(color = ~FDR, name = "Barcodes") %>%
  add_trace(y = plot_array$Knee_threshold, mode = "lines", name = "Knee") %>%
  add_trace(y = plot_array$Inflection_threshold, mode = "lines", name = "Inflection") %>%
  layout(title = paste("UMI Elbow Plot:", id), 
         xaxis = list(type = "log"), yaxis = list(type = "log"))


## ----3-format_data-------------------------------------------------------------------------------
array <- CreateSeuratObject(array)
  
# add non-nuclear read proportion metadata
Mito_proportion <- Matrix::colSums(array[grepl("^mt-|-mt-", rownames(array)),]) / array$nCount_RNA
array <- AddMetaData(array, Mito_proportion, col.name = "Mito_proportion")

Ribo_proportion <- Matrix::colSums(array[grepl("rpl|rps", rownames(array)),]) / array$nCount_RNA
array <- AddMetaData(array, Ribo_proportion, col.name = "Ribo_proportion")

# add group information by string splitting ID, removing replicate
array <- AddMetaData(array, strsplit(id, "-")[[1]][1], col.name = "region")

head(array@meta.data, 3)

paste("Cells passing emptyDrops filter:", length(colnames(array)))
count_all("called barcode")


## ----4-complexity_filter, warning = FALSE, message = FALSE---------------------------------------

# Store pre-filter population metrics for plotting
plot_array <- cbind.data.frame(array$nCount_RNA, 
                               array$nFeature_RNA, 
                               ifelse(array$nFeature_RNA < 1000, "<1000 Features", ">1000 Features"),
                               colnames(array))
colnames(plot_array) <- c("UMIs", "Genes", "Filter", "Cell_ID")

# Remove all cells with fewer than 1000 UMIs
array <- subset(array, cells = colnames(array)[array$nFeature_RNA >= 1000])

# Output new QC metrics
paste("Cells passing 1000-feature complexity filter:", length(colnames(array)))
count_all("cell with >1000 features")

plot_ly(data = plot_array, x = ~UMIs, y = ~Genes, color = ~Filter, text = ~Cell_ID) %>%
  add_markers(marker = list(size = 2, sizemode = "area")) %>%
  layout(title = paste("1000 Feature Complexity Filter:", id))


## ----5-mito_filter, message = FALSE, warning = FALSE---------------------------------------------

# Store pre-filter population metrics for plotting
plot_array <- cbind.data.frame(array$nCount_RNA, 
                               array$Mito_proportion, 
                               ifelse(array$Mito_proportion >=
                                        quantile(array$Mito_proportion)[4]+5*IQR(array$Mito_proportion),
                                      ">Q3+5xIQR", "<Q3+5xIQR"),
                               colnames(array))
colnames(plot_array) <- c("UMIs", "Mitochondrial_proportion", "Filter", "Cell_ID")

# Remove all extreme high outliers for mitochondrial reads, defined as Q3+5*IQR, which is an EXTREMELY lenient standard 
array <- subset(array, cells = colnames(array)[array$Mito_proportion <=
                                                 quantile(array$Mito_proportion)[4]+5*IQR(array$Mito_proportion)])

paste("Cells passing mitochondrial outlier filter:", length(colnames(array)))
count_all("cell passing mito filter")

plot_ly(data = plot_array, x = ~UMIs, y = ~Mitochondrial_proportion, color = ~Filter, text = ~Cell_ID) %>%
  add_markers(marker = list(size = 2, sizemode = "area")) %>%
  layout(title = paste("Mitochondrial Outlier Filter:", id))


## ----6-ribo_filter, message = FALSE, warning = FALSE---------------------------------------------

# Store pre-filter population metrics for plotting
plot_array <- cbind.data.frame(array$nCount_RNA, 
                               array$Ribo_proportion, 
                               ifelse(array$Ribo_proportion >=
                                        quantile(array$Ribo_proportion)[4]+5*IQR(array$Ribo_proportion),
                                      ">Q3+5xIQR", "<Q3+5xIQR"),
                               colnames(array))
colnames(plot_array) <- c("UMIs", "Ribosomal_proportion", "Filter",  "Cell_ID")

# Remove all extreme high outliers for mitochondrial reads, defined as Q3+5*IQR, which is an EXTREMELY lenient standard 
array <- subset(array, cells = colnames(array)[array$Ribo_proportion <=
                                                 quantile(array$Ribo_proportion)[4]+5*IQR(array$Ribo_proportion)])

paste("Cells passing ribosomal outlier filter:", length(colnames(array)))
count_all("cell passing ribo filter")

plot_ly(data = plot_array, x = ~UMIs, y = ~Ribosomal_proportion, color = ~Filter, text = ~Cell_ID) %>%
  add_markers(marker = list(size = 2, sizemode = "area")) %>%
  layout(title = paste("Ribosomal Outlier Filter:", id))


## ----7-doublet_filter, message = FALSE, warning = FALSE------------------------------------------
sce <- as.SingleCellExperiment(array) %>%
    scDblFinder

plot_array <- cbind.data.frame(array$nCount_RNA, 
                               array$nFeature_RNA,
                               sce$scDblFinder.class,
                               colnames(array))
colnames(plot_array) <- c("UMIs", "Features", "Doublet_status", "Cell_ID")

array <- subset(array, cells = colnames(array)[sce$scDblFinder.class == "singlet"])

paste("Cells passing doublet filter:", length(colnames(array)))
count_all("cell passing doublet filter")

plot_ly(data = plot_array, x = ~UMIs, y = ~Features, color = ~Doublet_status, text = ~Cell_ID) %>%
  add_markers(marker = list(size = 2, sizemode = "area")) %>%
  layout(title = paste("Doublet Filter:", id))


## ----8-outlier_filter, warning = FALSE, message = FALSE------------------------------------------

# Store pre-filter population metrics for plotting
plot_array <- cbind.data.frame(array$nCount_RNA, 
                               array$nFeature_RNA, 
                               ifelse(array$nCount_RNA > median(array$nCount_RNA)+5*mad(array$nCount_RNA), 
                                      "High Outlier", "Non-outlier"),
                               colnames(array))
colnames(plot_array) <- c("UMIs", "Genes", "Filter", "Cell_ID")

# Remove all cells with fewer than 1000 UMIs
array <- subset(array, cells = colnames(array)[array$nFeature_RNA >= 1000])

# Output new QC metrics
paste("Final cells passing pre-processing:", length(colnames(array)))
count_all("cell passing all filters")

plot_ly(data = plot_array, x = ~UMIs, y = ~Genes, color = ~Filter) %>%
  add_markers(marker = list(size = 2, sizemode = "area")) %>%
  layout(title = paste("High Outlier Filter:", id))


## ----9-save_output, warning = FALSE, message = FALSE---------------------------------------------
saveRDS(array, file = paste0("/datasets/preprocessed/", id_general, '/' , id, "_preprocessed.rds"))

