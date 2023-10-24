library(ggplot2)
library(cowplot)

#sample_name <- "OMB0793_Quad_Enth_T"
# slide_tissue_lowres_png <- sprintf("003-snakemake/results/spaceranger_count/%s/outs/spatial/tissue_lowres_image.png", sample_name)
# slide_predictions_png <- sprintf("003-snakemake/figures/seurat_transfer/predictions/%s", sample_name)

slide_tissue_lowres_png <- snakemake@input[["lowres_png"]]
slide_predictions_png <- snakemake@input[["predictions_png"]]
slide_qcmetrics_png <- snakemake@input[["qcmetrics_png"]]
slide_top_prediction_png <- snakemake@input[["slide_top_prediction_png"]]
montage_png <- snakemake@output[["montage_png"]]

# 1 2 4 4
# 3 3 4 4

gg_1 <- ggdraw() +
  draw_image(image = slide_tissue_lowres_png)
gg_2 <- ggdraw() +
  draw_image(image = slide_predictions_png)
gg_3 <- ggdraw() +
  draw_image(image = slide_qcmetrics_png)
gg_4 <- ggdraw() +
  draw_image(image = slide_top_prediction_png)
gg_12 <- plot_grid(
  gg_1,
  gg_2,
  ncol = 2, rel_widths = c(0.55, 0.45)
)
gg_123 <- plot_grid(
  gg_12,
  gg_3,
  nrow = 2, rel_heights = c(0.7, 0.3)
)
gg_1234 <- plot_grid(
  gg_123,
  gg_4,
  ncol = 2, rel_widths = c(1, 1)
)
ggsave(montage_png, gg_1234, width = 20, height = 8)
