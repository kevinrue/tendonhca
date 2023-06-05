library(jsonlite)
library(fgsea)

options(width = 200)

# snakemake@input[[1]]

# prepare C8 gene sets ----

c8_json <- jsonlite::fromJSON("data/c8.all.v2023.1.Hs.json")
# length(c8_json)

extract_symbols <- function(x) {
  unique(x[["geneSymbols"]])
}

format_leadingEdge <- function(x) {
    format_one_value <- function(genes) {
        paste0(genes, collapse = ";")
    }
    unlist(lapply(x, format_one_value))
}

pathways <- lapply(c8_json, extract_symbols)
# summary(lengths(pathways))
# str(head(pathways))

# prepare feature stats ----

rank_genes_groups <- read.table("results/filtered_genes/OMB1556_Ach_MTJ_H/markers/rank_genes_groups.tsv.gz", sep = "\t", header = TRUE)
# head(rank_genes_groups)

file_out <- "results/filtered_genes/OMB1556_Ach_MTJ_H/markers/fgsea.tsv"
header <- paste0(c("group", "pathway", "pval", "padj", "log2err", "ES", "NES", "size",  "leadingEdge"), collapse = "\t")
writeLines(header, file_out)

for (for_group in unique(rank_genes_groups$group)) {
    message("Processing group: ", for_group)
    # fetch feature stats for the group
    rank_genes <- subset(rank_genes_groups, group == for_group)
    feature_stats <- rank_genes$scores
    names(feature_stats) <- rank_genes$gene_name
    feature_stats <- sort(feature_stats, decreasing = TRUE)
    feature_stats <- feature_stats[!duplicated(names(feature_stats))]
    # fgsea
    set.seed(42)
    fgsea_res <- fgsea(pathways = pathways, 
        stats    = feature_stats,
        minSize  = 5,
        maxSize  = 500)
    fgsea_hits <- subset(fgsea_res, padj < 0.05 & NES > 0)
    if (nrow(fgsea_hits) == 0) {
        next
    }
    fgsea_hits <- as.data.frame(fgsea_hits)
    fgsea_hits <- cbind(group = for_group, fgsea_hits)
    fgsea_hits$leadingEdge <- format_leadingEdge(fgsea_hits$leadingEdge)
    write.table(fgsea_hits, file = file_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}
