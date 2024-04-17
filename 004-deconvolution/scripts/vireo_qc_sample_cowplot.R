library(ggplot2)
library(cowplot)

donor_umi <- snakemake@input[["donor_umi"]]
donor_n_vars <- snakemake@input[["donor_n_vars"]]
prob_donor_doublet <- snakemake@input[["prob_donor_doublet"]]
montage <- snakemake@output[["montage"]]

print(donor_umi)
print(donor_n_vars)
print(prob_donor_doublet)
print(montage)

# layout:
# 1
# 2
# 3

gg_1 <- ggdraw() +
  draw_image(image = donor_umi)
gg_2 <- ggdraw() +
  draw_image(image = donor_n_vars)
gg_3 <- ggdraw() +
  draw_image(image = prob_donor_doublet)
gg_123 <- plot_grid(
  gg_1,
  gg_2,
  gg_3,
  ncol = 1
)

ggsave(montage, gg_123, width = 7, height = 15)
