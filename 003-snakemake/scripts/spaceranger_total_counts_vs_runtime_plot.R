library(ggplot2)
library(ggrepel)

total_counts <- read.table("results/spaceranger_stats/total_counts.tsv", header = TRUE)
runtime <- read.table("results/spaceranger_stats/runtime.tsv", header = TRUE)
full_data <- merge(x = total_counts, y = runtime)
full_data$runtime_seconds <- full_data$hours * 3600 + full_data$minutes * 60 + full_data$seconds
full_data$runtime_hours <- full_data$runtime_seconds / 3600

full_data$total_counts_outlier <- full_data$total_counts > median(full_data$total_counts) + 5 * mad(full_data$total_counts)

ggplot(full_data, aes(runtime_hours, total_counts)) +
  geom_smooth(method = "lm", alpha = 0.3) +
  geom_point() +
  geom_label_repel(aes(label = sample_name), subset(full_data, total_counts_outlier)) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  labs(
    y = "Total counts",
    x = "Runtime (hours)"
  ) +
  theme_minimal()
