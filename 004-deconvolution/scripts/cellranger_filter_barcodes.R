library(DropletUtils)

cellranger_h5 <- snakemake@input[["h5"]]
umi_cutoff_min <- snakemake@params[["umi_cutoff_min"]]
barcodes <- snakemake@output[["barcodes"]]

sce <- read10xCounts(
  samples = cellranger_h5,
  col.names = TRUE
)

sce <- sce[, colSums(assay(sce, "counts")) > umi_cutoff_min]

barcodes_gz <- gzfile(barcodes, "w")
write(x = colnames(sce), file = barcodes_gz)
close(barcodes_gz)
