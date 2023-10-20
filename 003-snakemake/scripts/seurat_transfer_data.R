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
out_dir <- dirname(predictions_png)

# Fixed for now
use_pcs <- 15
use_resolution <- 0.3

message("=== Load the slide as a Seurat object ===")
seurat_slide <- Load10X_Spatial(spaceranger_out_dir)
seurat_slide

message("=== Remove spots with fewer than 50 UMI ===")
seurat_slide <- subset(seurat_slide, nCount_Spatial > 50)
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

reference_label <- Idents(reference)
levels(reference_label)[c(3, 4, 5)] <- levels(reference_label)[c(2, 1, 1)]
levels(reference_label)[c(1, 2)] <- c("Muscle cells", "Fibroblasts")
reference_symbols$label <- reference_label
Idents(reference_symbols) <- "cell_type"

message("=== Run standard workflow on reference sample ===")
reference_symbols <- SCTransform(reference_symbols, ncells = 3000, verbose = FALSE)
reference_symbols <- RunPCA(reference_symbols, verbose = FALSE)
reference_symbols <- RunUMAP(reference_symbols, dims = 1:30)

message("=== Transfer data (make predictions) ===")
anchors <- FindTransferAnchors(reference = reference_symbols, query = seurat_slide, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = Idents(reference_symbols),
                                  prediction.assay = TRUE,
                                  weight.reduction = seurat_slide[["pca"]], dims = 1:30)
seurat_slide[["predictions"]] <- predictions.assay

message("=== Save predictions to file ===")
DefaultAssay(seurat_slide) <- "predictions"
p <- SpatialFeaturePlot(
  seurat_slide,
  features = levels(reference_label),
  pt.size.factor = 1.6,
  alpha = c(0.25, 0.5),
  ncol = 3,
  crop = TRUE,
  min.cutoff = 0, max.cutoff = 1
)
ggsave(p, predictions_png, width=12, height=16)
