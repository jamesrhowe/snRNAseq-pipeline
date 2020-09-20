# contains helper functions for analysis notebooks

# this used to collect summary stats for every dataset
count_all <- function(step, object){
  cat(count_stats(paste("counts per", step), object@meta.data$nCount_RNA),
      count_stats(paste("features per", step), object@meta.data$nFeature_RNA),
      count_stats(paste("mitochondrial proportion per", step), object@meta.data$Mito_proportion),
      count_stats(paste("ribosomal proportion per", step), object@meta.data$Ribo_proportion), sep = "\n")
}

count_stats <- function(x, y){
  cat(paste("Min", x, min(y)),
      paste("Max", x, max(y)),
      paste("Mean", x, mean(y)),
      paste("Median", x, median(y)),
      paste("Std. dev.", x, sd(y)),
      paste("Std. error", x, sd(y)/sqrt(length(y))),
      paste(),
      sep = "\n")
}