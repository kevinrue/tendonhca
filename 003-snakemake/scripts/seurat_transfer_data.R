library(Seurat)

# === Inputs ===
spaceranger_h5_file <- snakemake@input
spaceranger_out_dir <- gsub("/raw_feature_bc_matrix.h5", "", spaceranger_h5_file)
sample_name <- basename(dirname(spaceranger_out_dir))

reference_rds <- "/ceph/project/tendonhca/albrecht/003-snakemake/data/HAMSTRING_singlets_ambRNA0.2_res0.15.RDS"

# === Outputs ===
predictions_png <- snakemake@output[["predictions_png"]]
out_dir <- dirname(predictions_png)

# Fixed for now
use_pcs <- 15
use_resolution <- 0.3

# Load the slide as a Seurat object
seurat_slide <- Load10X_Spatial(spaceranger_out_dir)
seurat_slide

# Remove spots with fewer than 50 UMI
seurat_slide <- subset(seurat_slide, nCount_Spatial > 50)
seurat_slide

# Run SCT on the slide data
seurat_slide <- SCTransform(seurat_slide, assay = "Spatial", verbose = TRUE)
seurat_slide

# Run standard workflow
seurat_slide <- FindNeighbors(seurat_slide, reduction = "pca", dims = 1:use_pcs)
seurat_slide <- FindClusters(seurat_slide, resolution = use_resolution, verbose = FALSE)
seurat_slide <- RunUMAP(seurat_slide, reduction = "pca", dims = 1:use_pcs)

# Plot clusters in spatial layout
p1 <- DimPlot(seurat_slide, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(seurat_slide, label = TRUE, label.size = 3, alpha = c(0.5, 0.5))
ggsave(
    file.path(out_dir, sprintf("dimplot_spatialplot_cluster_%s.png", sample_name)),
    p1 + p2,
    width = 12,
    height = 6
)

# Load reference data
reference <- readRDS(reference_rds)
reference

# Load and preprocess feature mapping table
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

# Make a new Seurat object after renaming ENSEMBL to SYMBOL
reference_counts <- reference@assays$RNA@counts
reference_label <- Idents(reference)
levels(reference_label)[c(3, 4, 5)] <- levels(reference_label)[c(2, 1, 1)]
levels(reference_label)[c(1, 2)] <- c("Muscle cells", "Fibroblasts")

reference_counts <- reference_counts[intersect(rownames(reference_counts), features_tsv$ENSEMBL), ]
rownames(reference_counts) <- features_tsv[rownames(reference_counts), "SYMBOL"]
reference_symbols <- CreateSeuratObject(reference_counts)