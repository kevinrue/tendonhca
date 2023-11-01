library(ggplot2)
library(cowplot)

#sample_name <- "OMB0793_Quad_Enth_T"
# slide_tissue_lowres_png <- sprintf("003-snakemake/results/spaceranger_count/%s/outs/spatial/tissue_lowres_image.png", sample_name)
# slide_predictions_png <- sprintf("003-snakemake/figures/seurat_transfer/predictions/%s", sample_name)

slide_tissue_lowres_png <- snakemake@input[["lowres_png"]]
slide_predictions_png <- snakemake@input[["predictions_png"]]
slide_probabilities_diff_png <- snakemake@input[["probabilities_diff_png"]]
slide_qcmetrics_png <- snakemake@input[["qcmetrics_png"]]
slide_top_prediction_png <- snakemake@input[["slide_top_prediction_png"]]
montage_png <- snakemake@output[["montage_png"]]

# TODO: update
# 1 1 2 2 5 5 5
# 1 1 3 3 5 5 5
# 4 4 4 4 5 5 5

gg_1 <- ggdraw() +
  draw_image(image = slide_tissue_lowres_png)
gg_2 <- ggdraw() +
  draw_image(image = slide_predictions_png)
gg_3 <- ggdraw() +
  draw_image(image = slide_probabilities_diff_png)
gg_4 <- ggdraw() +
  draw_image(image = slide_qcmetrics_png)
gg_5 <- ggdraw() +
  draw_image(image = slide_top_prediction_png)
gg_23 <- plot_grid(
  gg_2,
  gg_3,
  nrow = 2, rel_heights = c(0.6, 0.4)
)
gg_123 <- plot_grid(
  gg_1,
  gg_23,
  ncol = 2, rel_widths = c(0.55, 0.45)
)
gg_1234 <- plot_grid(
  gg_123,
  gg_4,
  nrow = 2, rel_heights = c(0.7, 0.3)
)
gg_12345 <- plot_grid(
  gg_1234,
  gg_5,
  ncol = 2, rel_widths = c(1, 1)
)
ggsave(montage_png, gg_12345, width = 20, height = 8)
