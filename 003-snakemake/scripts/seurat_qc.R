sample_names <- list.files("003-snakemake/fastqs")
# sample_name <- "OMB0793_Quad_Enth_T"

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
