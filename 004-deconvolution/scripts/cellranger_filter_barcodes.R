library(DropletUtils)

umi_cutoff_min <- snakemake@params[["umi_cutoff_min"]]

sce <- read10xCounts(
  samples = snakemake@input[["cellranger_filtered"]],
  col.names = TRUE
)

sce <- sce[, colSums(assay(sce, "counts")) > umi_cutoff_min]

write(x = colnames(sce), file = snakemake@output[["cellranger_filtered"]])