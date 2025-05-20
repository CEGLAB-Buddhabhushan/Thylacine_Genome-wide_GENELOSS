data <- read.table("S.harrisii.info.plot.tsv", header = T, sep = "\t")

mean_df <- data %>%
  group_by(TOGA_status) %>%
  summarise(mean_gc = mean(GC_Content, na.rm = TRUE))
pairwise.wilcox.test(data$GC_Content, data$TOGA_status, p.adjust.method = "BH")
# Define pairwise combinations of interest
comparisons <- list(
  c("L", "I"),
  c("M", "I"),
  c("PI", "I"),
  c("PM", "I"),
  c("UL", "I"),
  c("M", "L"),
  c("PI", "L"),
  c("PM", "L"),
  c("UL", "L"),
  c("PI", "M"),
  c("PM", "M"),
  c("UL", "M"),
  c("PM", "PI"),
  c("UL", "PI"),
  c("UL", "PM")
)


png("GC_content_raincloud_plot_with_pairwise_wilcox.png", units="in", width=16, height=9, res=900)

# Basic raincloud plot
ggplot(data, aes(x = TOGA_status, y = GC_Content, fill = TOGA_status)) +
  
  # Half-violin plot (distribution)
  geom_half_violin(
    aes(fill = TOGA_status),
    side = "l", 
    alpha = 0.6,
    color = NA,
    trim = FALSE
  ) +
  
  # Boxplot in the middle
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.7
  ) +

  # Jittered individual data points
  geom_half_point(
    side = "r",
    shape = 21,
    size = 1,
    alpha = 0.4,
    position = position_jitter(width = 0.07)
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
               fill = "white", color = "black", stroke = 1.2) +
  geom_text(data = mean_df, aes(x = TOGA_status, y = mean_gc, 
                                label = sprintf("%.2f", mean_gc)),
            vjust = 1,hjust =1.5,  size = 4, color = "black") +
  labs(x = "TOGA Status",
    y = "GC Content (%)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()


