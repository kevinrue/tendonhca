sample_names <- list.files("003-snakemake/fastqs")
# sample_name <- "OMB0793_Quad_Enth_T"

mito_gene_data <- read.table("003-snakemake/annotations/genes_mitochondrial.tsv", header = TRUE)
mito_gene_names <- mito_gene_data$gene_name
# head(mito_gene_names)

ribo_gene_data <- read.table("003-snakemake/annotations/genes_ribosomal.tsv", header = TRUE)
ribo_gene_names <- ribo_gene_data$gene_name
# head(ribo_gene_names)

library(Seurat)
library(ggplot2)

plot_data <- data.frame(
  sample_name = character(0),
  total_umi = integer(0)
)

for (sample_name in sample_names) {
  spaceranger_out_dir <- sprintf("003-snakemake/results/spaceranger_count/%s/outs", sample_name)

  message("=== Load the slide as a Seurat object ===")
  message("Sample: ", sample_name)
  seurat_slide <- Load10X_Spatial(spaceranger_out_dir)
  seurat_slide

  message("=== Format UMI data table ===")
  seurat_plot_data <- data.frame(
    sample_name = sample_name,
    total_umi = seurat_slide$nCount_Spatial
  )

  plot_data <- rbind(plot_data, seurat_plot_data)

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
  ggsave(sprintf("003-snakemake/tests/seurat_spatial_qc/%s.png", sample_name), p, width = 16, height = 4)
}

p <- ggplot(plot_data, aes(total_umi)) +
  geom_histogram(bins = 100, fill = "lightblue", color = "lightblue") + # binwidth = 100
  geom_vline(xintercept = 100, color = "red", linetype = "dashed", linewidth = 0.2) +
  facet_wrap(~sample_name) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_log10(expand = c(0, 0)) +
  cowplot::theme_cowplot() +
  theme(plot.background = element_rect(fill = "white"), strip.background = element_rect(fill = "white", color = "black"))

ggsave("003-snakemake/tests/seurat_qc_umi_histogram.png", p, width = 12, height = 10)
