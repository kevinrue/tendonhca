library(Seurat)
library(ggplot2)
library(patchwork)
library(magrittr)
library(dplyr)

# === Inputs ===
message("=== Process inputs ===")
spaceranger_h5_file <- snakemake@input
spaceranger_out_dir <- gsub("/raw_feature_bc_matrix.h5", "", spaceranger_h5_file)
sample_name <- basename(dirname(spaceranger_out_dir))

reference_rds <- "/ceph/project/tendonhca/albrecht/003-snakemake/data/HAMSTRING_singlets_ambRNA0.2_res0.15.RDS"

# === Outputs ===
message("=== Process outputs ===")
predictions_png <- snakemake@output[["predictions_png"]]
top_prediction_png <- snakemake@output[["top_prediction_png"]]
out_dir <- dirname(predictions_png)

# Fixed for now
use_pcs <- 15
use_resolution <- 0.3
min_umi <- 100

message("=== Load the slide as a Seurat object ===")
seurat_slide <- Load10X_Spatial(spaceranger_out_dir)
seurat_slide

message("=== Remove spots with fewer than 50 UMI ===")
seurat_slide <- subset(seurat_slide, nCount_Spatial >= min_umi)
seurat_slide

message("=== Run SCT on the slide data ===")
seurat_slide <- SCTransform(seurat_slide, assay = "Spatial", verbose = TRUE)
seurat_slide

# Run standard workflow
message("=== Run standard workflow ===")
seurat_slide <- RunPCA(seurat_slide, assay = "SCT", verbose = FALSE)
seurat_slide <- FindNeighbors(seurat_slide, reduction = "pca", dims = 1:use_pcs)
seurat_slide <- FindClusters(seurat_slide, resolution = use_resolution, verbose = FALSE)
seurat_slide <- RunUMAP(seurat_slide, reduction = "pca", dims = 1:use_pcs)

message("=== Plot clusters in spatial layout ===")
p1 <- DimPlot(seurat_slide, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(seurat_slide, label = TRUE, label.size = 3, alpha = c(0.5, 0.5))
ggsave(
  file.path(out_dir, sprintf("dimplot_spatialplot_cluster_%s.png", sample_name)),
  p1 + p2,
  width = 12,
  height = 6
)

message("=== Load reference data ===")
reference <- readRDS(reference_rds)
reference

message("=== Load and preprocess feature mapping table ===")
## NOTE: reference uses ENSEMBL while slide uses SYMBOL
sample_features_tsv <- sprintf(
  "results/spaceranger_count/%s/outs/filtered_feature_bc_matrix/features.tsv.gz",
  sample_name
)
message("Loading features.tsv.gz for sample ", sample_name)
features_tsv <- read.table(sample_features_tsv, sep = "\t")
features_tsv <- features_tsv[, 1:2]
colnames(features_tsv) <- c("ENSEMBL", "SYMBOL")
rownames(features_tsv) <- features_tsv$ENSEMBL

message("=== Make a new Seurat object after renaming ENSEMBL to SYMBOL ===")
reference_counts <- reference@assays$RNA@counts
reference_counts <- reference_counts[intersect(rownames(reference_counts), features_tsv$ENSEMBL), ]
rownames(reference_counts) <- features_tsv[rownames(reference_counts), "SYMBOL"]
reference_symbols <- CreateSeuratObject(reference_counts)

# Merge certain groups of cell type labels
reference_label <- Idents(reference)
levels(reference_label)[c(3, 4, 5, 7, 8, 9)] <- levels(reference_label)[c(2, 1, 1, 1, 6, 6)]
levels(reference_label)[c(1, 2, 3)] <- c("Muscle cells", "Fibroblasts", "Vessel")
reference_symbols$label <- reference_label
Idents(reference_symbols) <- "label"

message("=== Run standard workflow on reference sample ===")
reference_symbols <- SCTransform(reference_symbols, ncells = 3000, verbose = FALSE)
reference_symbols <- RunPCA(reference_symbols, verbose = FALSE)
reference_symbols <- RunUMAP(reference_symbols, dims = 1:30)

message("=== Transfer data (make predictions) ===")
anchors <- FindTransferAnchors(reference = reference_symbols, query = seurat_slide, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = Idents(reference_symbols),
                                  prediction.assay = TRUE,
                                  k.weight = 40, # default 50 is too high for some samples
                                  weight.reduction = seurat_slide[["pca"]], dims = 1:30)
seurat_slide[["predictions"]] <- predictions.assay

message("=== Save predictions to file ===")
DefaultAssay(seurat_slide) <- "predictions"
p <- SpatialFeaturePlot(
  seurat_slide,
  features = levels(reference_label),
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  ncol = 3,
  crop = TRUE,
  min.cutoff = 0, max.cutoff = 0.5
)

ggsave(predictions_png, p, width=12, height=16)

##
# Custom plot #
##

# For each spot,
# display the prediction associated with the highest probability,
# using alpha to indicate its probability (i.e., confidence).

image.use <- seurat_slide@images$slice1
coordinates <- GetTissueCoordinates(object = image.use)
plot_data <- data.frame(
  coordinates,
  max_celltype = rownames(seurat_slide@assays$predictions@data)[apply(seurat_slide@assays$predictions@data, MARGIN = 2, FUN = which.max)],
  max_probability = apply(seurat_slide@assays$predictions@data, MARGIN = 2, FUN = max)
)

cell_types <- c(
  "Muscle cells",
  "Fibroblasts",
  "Vessel",
  "Satellite cells",
  "Adipocytes",
  "Immune cells",
  "Nerve cells"
)

n <- length(cell_types)
hues <- seq(15, 375, length=(n + 1))
fixed_colors <- hcl(h=hues, l=65, c=100)[seq_len(n)]
names(fixed_colors) <- cell_types

p <- SingleSpatialPlot(
  data = plot_data,
  image = seurat_slide@images$slice1,
  image.alpha = 0.3,
  pt.size.factor = 1,
  col.by = "max_celltype",
  alpha.by = "max_probability"
) +
  scale_alpha(limits = c(0, 0.5)) +
  scale_fill_manual(values = fixed_colors) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.key.size = unit(2, "lines")
  )

ggsave(top_prediction_png, p, width=16, height=12)
