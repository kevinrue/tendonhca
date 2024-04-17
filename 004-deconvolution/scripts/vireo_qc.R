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
  theme(panel.grid.minor.y = element_blank()) +
  labs(
    title = snakemake@params[["sample"]],
    y = "Total UMI",
    x = "Donor ID"
  )

ggsave(
  filename = snakemake@output[["donor_umi"]],
  plot = gg,
  width = 7,
  height = 5
)

gg <- donor_ids_data %>% 
  ggplot(aes(donor_id, n_vars)) +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank()) +
  labs(
    title = snakemake@params[["sample"]],
    y = "Number of variants",
    x = "Donor ID"
  )

ggsave(
  filename = snakemake@output[["donor_n_vars"]],
  plot = gg,
  width = 7,
  height = 5
)

gg <- donor_ids_data %>%
  ggplot(aes(prob_doublet, prob_max)) +
  geom_point() +
  facet_wrap(~donor_id) +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank()) +
  labs(
    title = snakemake@params[["sample"]],
    y = "Probability of assigned donor",
    x = "Probability of doublet"
  )

ggsave(
  filename = snakemake@output[["prob_donor_doublet"]],
  plot = gg,
  width = 7,
  height = 7
)