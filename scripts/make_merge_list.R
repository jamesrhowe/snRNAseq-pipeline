array <- list()
for (i in 1:length(snakemake@input)){
  array[[i]] <- readRDS(snakemake@input[[i]])
}
saveRDS(array, file = snakemake@output[[1]])
