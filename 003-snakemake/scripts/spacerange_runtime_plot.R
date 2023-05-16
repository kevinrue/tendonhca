library(ggplot2)

x <- read.table("results/spaceranger_stats/runtime.tsv", header = TRUE)

x$total_seconds <- x$hours*60*60 + x$minutes * 60 + x$seconds
x$total_hours <- x$total_seconds / (60*60)

ggplot(x) +
  geom_col(aes(x = sample_name, y = total_hours)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    y = "Runtime (hours)"
  )
