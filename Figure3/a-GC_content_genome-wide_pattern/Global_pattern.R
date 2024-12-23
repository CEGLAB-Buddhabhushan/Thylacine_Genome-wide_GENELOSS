library(ggplot2)
library(dplyr)
library(ggpubr)
# Load your data
df <- read.table("S.harrisii.info.tsv", header = TRUE, sep = '\t')


# Define colors and shapes
TOGA_Status <- c("I", "L", "M", "PI", "PM", "UL")
colors <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#b15928")
shapes <- c(1, 17, 15, 3, 16, 18)

# Create mapping for colors and shapes
color_mapping <- setNames(colors, TOGA_Status)
shape_mapping <- setNames(shapes, TOGA_Status)

# Plot
p1 <- ggplot(df, aes(x = GC_Stretch, y = GC_Content, color = `TOGA_status`, shape = `TOGA_status`)) +
  geom_point(size = 3) +
  scale_color_manual(values = color_mapping, name = "TOGA_Status") +
  scale_shape_manual(values = shape_mapping, name = "TOGA_Status") +
  labs(x = "Average Length of G/C-stretches",
       y = "GC Content (%)") +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), 
        axis.title = element_text(face = "bold")) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE),
         shape = guide_legend(ncol = 1, byrow = TRUE))+
  theme_classic() +
  scale_fill_brewer(palette = "Set3")+
  facet_wrap(~ TOGA_status)

tiff('GC_Content_vs_GC-stretch.tiff', units="in", width=16, height=9, res=900, compression = 'lzw')
print(p1)
dev.off()

mean_values <- df %>%
  group_by(TOGA_status) %>%
  summarize(mean_GC_Content = mean(GC_Content, na.rm = TRUE))

# Create the first plot for TOGA_status
p2 <- ggplot(df, aes(x = TOGA_status, y = GC_Content, fill = TOGA_status)) +
  geom_violin() +
  geom_boxplot(width = 0.2, fill = "white", color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "black", fill = "red") + # Add mean points
  geom_text(data = mean_values, aes(label = round(mean_GC_Content, 2), y = mean_GC_Content), 
            vjust = -0.5, color = "black", size = 4) + # Add text labels for means
  labs(x = "Gene Status (TOGA Prediction)",
       y = "GC Content (%)") +
  geom_hline(yintercept = 55, linetype = "dashed", color = "red", size = 1) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

# Convert Missing_exon and Deleted_exon to factors
df$Missing_exon <- factor(df$Missing_exon)
df$Deleted_exon <- factor(df$Deleted_exon)

# Create a combined status variable for better visualization
df$Exon_Status <- interaction(df$Missing_exon, df$Deleted_exon)

# Calculate mean GC content for Exon_Status
summary_stats <- df %>%
  group_by(Exon_Status) %>%
  summarise(mean_gc = mean(GC_Content, na.rm = TRUE))

# Create the second plot for Exon_Status
p3 <- ggplot(df, aes(x = Exon_Status, y = GC_Content, fill = Exon_Status)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "black", fill = "red") + # Add mean points
  geom_text(data = summary_stats, aes(label = round(mean_gc, 2), y = mean_gc), 
            vjust = -0.5, color = "black", size = 4) + # Corrected to use mean_gc for labels
  labs(x = "Exon Status (Missing & Deleted)",
       y = "GC Content (%)") +
  geom_hline(yintercept = 55, linetype = "dashed", color = "red", size = 1) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

# Arrange the plots vertically
p23<-ggarrange(p2, p3, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
tiff('Desity_plot_GC.tiff', units="in", width=16, height=9, res=900, compression = 'lzw')

ggarrange(p2, p3, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

dev.off()

