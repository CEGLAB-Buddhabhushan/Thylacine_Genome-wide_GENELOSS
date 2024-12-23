library(devtools)
library(ggside)
library(tidyverse)
library(tidyquant)
library(dplyr)
library(ggplot2)
library(ggrepel)

df <- read.table("S.harrisii.info.plot.tsv", header = TRUE, sep = '\t')
highlight_genes <- c("CUZD1", "VWA7", "SAMD9L", "HSD17B13")
q3_values <- df %>%
  group_by(TOGA_status) %>%
  summarize(
    Q3_GC_Stretch = quantile(GC_Stretch, 0.75, na.rm = TRUE),
    Q3_GC_Content = quantile(GC_Content, 0.75, na.rm = TRUE)
  )
p2 <- df %>%
  ggplot(aes(GC_Stretch, GC_Content, color = TOGA_status)) +
  geom_point(size = 2, alpha = 0.2) +
  geom_xsidedensity(
    aes(
      y    = after_stat(density),
      fill = TOGA_status
    ),
    alpha    = 0.5,
    size     = 1,
    position = "stack"
  ) +
  geom_ysidedensity(
    aes(
      x    = after_stat(density),
      fill = TOGA_status
    ),
    alpha    = 0.5,
    size     = 1,
    position = "stack"
  ) +
  scale_color_tq() +
  scale_fill_tq() +
  theme_tq() +
  labs(x = "Average Length of G/C-stretches", y = "GC Content (%)") +
  theme(
    ggside.panel.scale.x = 0.4,
    ggside.panel.scale.y = 0.4
  ) +
  # Add horizontal and vertical lines for Q3
  geom_vline(data = q3_values, aes(xintercept = Q3_GC_Stretch, color = TOGA_status), 
             linetype = "dashed", size = 1) +
  geom_hline(data = q3_values, aes(yintercept = Q3_GC_Content, color = TOGA_status), 
             linetype = "dashed", size = 1) +
  # Highlight specific genes with red points and adjust label positioning with ggrepel
  geom_point(data = df %>% filter(Gene_name %in% highlight_genes), 
             aes(alpha = 0.9, fill = "red"), size = 3, shape = 17) +  # Adjust point size
  # Use ggrepel to prevent label overlap
  geom_text_repel(data = df %>% filter(Gene_name %in% highlight_genes), 
                  aes(label = Gene_name), 
                  color = "red", 
                  fontface = "bold.italic", 
                  size = 4, 
                  box.padding = 0.3, 
                  point.padding = 0.5, 
                  max.overlaps = 1000,  # Avoid excessive overlap
                  segment.size = 0.5)  # Adjust the arrow size if using segments


png('Global_pattern_lost_genes.png', units="in", width=8, height=5, res=900)
print(p2)
dev.off()
