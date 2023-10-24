# sample_names <- list.files("003-snakemake/fastqs")
# sample_name <- "OMB0793_Quad_Enth_T"

# === Inputs ===
message("=== Process inputs ===")
spaceranger_h5_file <- snakemake@input
spaceranger_out_dir <- gsub("/raw_feature_bc_matrix.h5", "", spaceranger_h5_file)
sample_name <- basename(dirname(spaceranger_out_dir))

mito_gene_data <- read.table("annotations/genes_mitochondrial.tsv", header = TRUE)
mito_gene_names <- mito_gene_data$gene_name
# head(mito_gene_names)

ribo_gene_data <- read.table("annotations/genes_ribosomal.tsv", header = TRUE)
ribo_gene_names <- ribo_gene_data$gene_name
# head(ribo_gene_names)

# === Outputs ===
message("=== Process outputs ===")
histograms_png <- snakemake@output[["histograms_png"]]
out_dir <- dirname(histograms_png)

library(Seurat)
library(ggplot2)

message("=== Load the slide as a Seurat object ===")
message("Sample: ", sample_name)
seurat_slide <- Load10X_Spatial(spaceranger_out_dir)
seurat_slide

message("=== Format UMI data table ===")
seurat_plot_data <- data.frame(
  sample_name = sample_name,
  total_umi = seurat_slide$nCount_Spatial
)

message("=== Compute mitochondrial content ===")
seurat_slide$nCount <- seurat_slide$nCount_Spatial
seurat_slide$nFeature <- seurat_slide$nFeature_Spatial
seurat_slide$mito_pct <- colSums(seurat_slide@assays$Spatial@counts[mito_gene_names,]) / colSums(seurat_slide@assays$Spatial@counts)
seurat_slide$ribo_pct <- colSums(seurat_slide@assays$Spatial@counts[ribo_gene_names,]) / colSums(seurat_slide@assays$Spatial@counts)

p1 <- SpatialFeaturePlot(seurat_slide, features = "nCount") +
  scale_fill_viridis_c(limits = c(0, 1000), oob = scales::squish) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 6)
  )
p1
p2 <- SpatialFeaturePlot(seurat_slide, features = "nFeature") +
  scale_fill_viridis_c(limits = c(0, 1000), oob = scales::squish) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 6)
  )
p3 <- SpatialFeaturePlot(seurat_slide, features = "mito_pct") +
  scale_fill_viridis_c(limits = c(0, 0.1), oob = scales::squish) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 6)
  )
p4 <- SpatialFeaturePlot(seurat_slide, features = "ribo_pct") +
  scale_fill_viridis_c(limits = c(0, 0.2), oob = scales::squish) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 6)
  )
p <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 4)
ggsave(histograms_png, p, width = 16, height = 4)
