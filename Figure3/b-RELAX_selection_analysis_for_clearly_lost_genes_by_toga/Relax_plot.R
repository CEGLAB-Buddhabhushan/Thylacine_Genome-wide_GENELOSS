setwd("/home/ceglab358/BUDDHA/Tasmanian_wolf/RELAX_selection_L_fig/")

library(ggplot2)
library(ggrepel)

# Manually specify the colors for each group
group_colors <- c(
  "Highlight" = "red", 
  "Lost" = "red", 
  "Low quality protein" = "black", 
  "OR" = "purple",      
  "Pseudogene" = "green", 
  "Unclear" = "black"   
)

# Manually set the shapes for each group
group_shapes <- c(
  "Highlight" = 16,
  "Lost" = 16, 
  "Low quality protein" = 16, 
  "OR" = 15,          
  "Pseudogene" = 17, 
  "Unclear" = 16      
)

# Create the plot

p<-ggplot(data, aes(x = K, y = log_p_value, color = Group, shape = Group)) +
  geom_point(size = 3) + 
  scale_color_manual(values = group_colors) + 
  scale_shape_manual(values = group_shapes) +  
  geom_vline(xintercept = 1, color = "black", linetype = "dashed", size = 1) +  # Vertical line at x = 1
  geom_hline(yintercept = 1.3, color = "black", linetype = "dashed", size = 1) +  # Horizontal line at y = 1.3
  labs(
    x = "Relaxation/intensification parameter (K)",
    y = "-log10(p-value)"
  ) +
  theme(legend.position = "none") +  # Position the legend at the top
  geom_label_repel(
    data = subset(data, Group == "Highlight"),  # Add labels only for "Highlight" group
    aes(label = Gene),  # Label the points with Gene names
    color = "red", 
    fontface = "italic", 
    size = 4, 
    box.padding = 0.9,  # Padding around the box
    point.padding = 0,  # Padding around the points
    segment.color = "grey50",  # Color of the connecting lines
    max.overlaps = 10000,  # Limit the number of overlapping labels
    arrow = arrow(type = "closed", length = unit(0.1, "inches"))  # Add arrows pointing to the text
  ) +
  annotate("text", x = 1, y = 3.2, label = "K = 1", color = "black", size = 3, angle = 90, vjust = -1) +
  annotate("text", x = 4, y = 1.3, label = "p-value = 0.05", color = "black", size = 3, angle = 0, vjust = -1) +
  coord_cartesian(xlim = c(-1, 12), ylim = c(-0.5, 5)) + theme_publish()+
  theme(legend.position = "none")
ggsave("Relax_plot.png", plot = p, dpi = 900, width = 6, height = 6, units = "in")

