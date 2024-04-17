library(DropletUtils)
library(ggplot2)
library(tidyverse)

sce <- read10xCounts(
  samples = snakemake@input[["cellranger_filtered"]],
  col.names = TRUE
)

donor_ids_data <- read_delim(
  file = snakemake@input[["vireo_donor_ids"]],
  delim = "\t"
)

gg <- donor_ids_data %>% 
  mutate(
    total_umi = colSums(assay(sce[, donor_ids_data$cell], "counts"))
  ) %>% 
  ggplot(aes(donor_id, total_umi)) +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank())

ggsave(
  filename = snakemake@output[["donor_umi"]],
  plot = gg,
  width = 7,
  height = 5
)

gg <- donor_ids_data %>% 
  mutate(
    total_umi = colSums(assay(sce[, donor_ids_data$cell], "counts"))
  ) %>% 
  ggplot(aes(donor_id, n_vars)) +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank())

ggsave(
  filename = snakemake@output[["donor_n_vars"]],
  plot = gg,
  width = 7,
  height = 5
)
