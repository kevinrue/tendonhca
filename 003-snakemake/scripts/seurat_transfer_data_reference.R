library(Seurat)
library(ggplot2)
library(patchwork)
library(magrittr)
library(dplyr)

sessionInfo()

# === Inputs ===
message("=== Process inputs ===")
spaceranger_h5_file <- snakemake@input
spaceranger_out_dir <- gsub("/raw_feature_bc_matrix.h5", "", spaceranger_h5_file)
sample_name <- basename(dirname(spaceranger_out_dir))

samples <- read.csv(snakemake@config[["samples"]], sep = "\t", header = TRUE)
rownames(samples) <- samples[["sample_name"]]

reference_rds <- samples[sample_name, "reference_rds", drop = TRUE]
message("Using reference: ", reference_rds)

annotation_metadata_column <- samples[sample_name, "annotation_metadata_column", drop = TRUE]
message("Using annotation metadata column: ", annotation_metadata_column)

# === Outputs ===
message("=== Process outputs ===")
predictions_png <- snakemake@output[["predictions_png"]]
top_prediction_png <- snakemake@output[["top_prediction_png"]]
probabilities_diff_png <- snakemake@output[["probabilities_diff_png"]]
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

message("=== Move cell type annotations to Idents() ===")
message("* annotation_metadata_column: ", annotation_metadata_column)
reference_label <- reference[[annotation_metadata_column]][[1]]
reference_label <- as.factor(reference_label)
reference$label <- reference_label
Idents(reference) <- "label"
table(Idents(reference))

message("=== Run standard workflow on reference sample ===")
reference <- SCTransform(reference, ncells = 3000, verbose = TRUE)
reference <- RunPCA(reference, verbose = TRUE)
reference <- RunUMAP(reference, dims = 1:30)

message("=== Transfer data (make predictions) ===")
anchors <- FindTransferAnchors(reference = reference, query = seurat_slide, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = Idents(reference),
                                  prediction.assay = TRUE,
                                  k.weight = 40, # arbitrary; default 50 is too high for some samples
                                  weight.reduction = seurat_slide[["pca"]],
                                  dims = 1:30 # arbitrary
                                  )
seurat_slide[["predictions"]] <- predictions.assay

message("=== Save predictions to file ===")
DefaultAssay(seurat_slide) <- "predictions"
p <- SpatialFeaturePlot(
  seurat_slide,
  features = levels(reference_label),
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  ncol = 5,
  crop = FALSE,
  min.cutoff = 0, max.cutoff = 0.5
)

ggsave(predictions_png, p, width=12, height=9)

##
# Custom plot #
##

# For each spot,
# display the prediction associated with the highest probability,
# using alpha to indicate its probability (i.e., confidence).

which.max.n <- function(x, n=1) {
  xo <- order(x, decreasing = TRUE)
  xo[n]
}

max.n <- function(x, n=1) {
  xo <- order(x, decreasing = TRUE)
  x[xo[n]]
}

# remove row "max"
prediction_data <- seurat_slide@assays$predictions@data[rownames(seurat_slide@assays$predictions@data) != "max", ]

image.use <- seurat_slide@images$slice1
coordinates <- GetTissueCoordinates(object = image.use)
plot_data <- data.frame(
  coordinates,
  max_celltype = rownames(prediction_data)[apply(prediction_data, MARGIN = 2, FUN = which.max.n, n = 1)],
  max_probability = apply(prediction_data, MARGIN = 2, FUN = max.n, n = 1),
  max2_celltype = rownames(prediction_data)[apply(prediction_data, MARGIN = 2, FUN = which.max.n, n = 2)],
  max2_probability = apply(prediction_data, MARGIN = 2, FUN = max.n, n = 2)
)
plot_data$diff_max_max2 <- plot_data$max_probability - plot_data$max2_probability

p <- ggplot(plot_data) +
  geom_histogram(aes(x = diff_max_max2), fill = "grey", color = "black", bins = 50) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw()
ggsave(
  probabilities_diff_png,
  p,
  width = 12,
  height = 6
)

plot_data <- subset(plot_data, diff_max_max2 > 0.5 & max_probability > 0.5)

# Abbreviated:
# cell_types <- c(
#   "Muscle cells",
#   "Fibroblasts",
#   "Vessel",
#   "Satellite cells",
#   "Adipocytes",
#   "Immune cells",
#   "Nerve cells"
# )

# fibroblasts = *fibroblasts {3 => 2}
# muscle = *muscle + satellite {4, 5, 7 => 1}
# vessel = mural + vascular + lymphatic + immune {8, 9, 11 => 6, 6, 6}

# Fibroblasts: 

# cell_types <- c(
#   "Fibroblasts",
#   "Macrophages",
#   "Vascular endothelial cells",
#   "Mural cells",
#   "Adipocytes",
#   "T cells",
#   "Nervous system cells",
#   "Lymphatic endothelial cells",
#   "Dividing fibroblasts / mural cells",
#   "Dendritic cells",
#   "Osteoblasts",
#   "Granulocytes",
#   "Osteoclasts",
#   "Dividing macrophages"
# )

# n <- length(cell_types)
# hues <- seq(15, 375, length=(n + 1))
# fixed_colors <- hcl(h=hues, l=65, c=100)[seq_len(n)]
# names(fixed_colors) <- cell_types
# names = LETTERS[1:n]

fixed_colors <- c(
  "Adipocytes" = "#F3C300", # D
  "Endothelial cells" = "#888888", # E
  "Immune cells" = "#E68FAC", # red
  "Fibroblasts" = "#E25822", # C
  "Lymphatic endothelial cells" = "#888888", # F
  "Macrophages" = "#E68FAC", # A
  "Muscle cells" = "#4E79A7", # B
  "Nervous system cells" = "#8DB600", # G
  "Vascular endothelial cells" = "#888888", # E
  "Mural cells" = "#00BC56", # F
  "T cells" = "#E68FAC", # N
  "Dividing fibroblasts / mural cells" = "#00BC56", # F
  "Dendritic cells" = "#E68FAC", # L
  "Osteoblasts" = "grey70",
  "Granulocytes" = "#E68FAC", # K
  "Osteoclasts" = "grey40", 
  "Dividing macrophages" = "#E68FAC" # A
)

p <- SingleSpatialPlot(
  data = plot_data,
  image = seurat_slide@images$slice1,
  image.alpha = 0.3,
  pt.size.factor = 1,
  col.by = "max_celltype",
  alpha.by = "max_probability",
  crop = FALSE
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
