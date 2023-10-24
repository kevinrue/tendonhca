library(ggplot2)
library(cowplot)

#sample_name <- "OMB0793_Quad_Enth_T"
# slide_tissue_lowres_png <- sprintf("003-snakemake/results/spaceranger_count/%s/outs/spatial/tissue_lowres_image.png", sample_name)
# slide_predictions_png <- sprintf("003-snakemake/figures/seurat_transfer/predictions/%s", sample_name)

slide_tissue_lowres_png <- snakemake@input[["lowres_png"]]
slide_predictions_png <- snakemake@input[["predictions_png"]]
slide_qcmetrics_png <- snakemake@input[["qcmetrics_png"]]
montage_png <- snakemake@output[["montage_png"]]

gg_topleft <- ggdraw() +
  draw_image(image = slide_tissue_lowres_png)
gg_topright <- ggdraw() +
  draw_image(image = slide_predictions_png)
gg_bottom <- ggdraw() +
  draw_image(image = slide_qcmetrics_png)
gg_top <- plot_grid(
  gg_topleft,
  gg_topright, ncol = 2, rel_widths = c(0.55, 0.45)
)
gg_grid <- plot_grid(
  gg_top,
  gg_bottom, nrow = 2, rel_heights = c(0.7, 0.3)
)
ggsave(montage_png, gg_grid, width = 10, height = 8)
