library(ggtree)
library(ggplot2)

phy_tree <- read.tree("./Species_name.nwk")
df_tip_data1 <- read.table("tree_SAMD9-9L_gene_status.tsv", sep = "\t", header = TRUE)
df_tip_data1 <- df_tip_data1 %>%
  mutate(
    SAMD9L_gene_status = ifelse(SAMD9L_gene_status == "", "No_chain", SAMD9L_gene_status),
    SAMD9_gene_status = ifelse(SAMD9_gene_status == "", "No_chain", SAMD9_gene_status)
  )
b <- ggtree(phy_tree,layout='circular') %<+% df_tip_data1 + 
  geom_hilight(node=109, fill="#ff8282", alpha=.5) + # Light Orange
  geom_hilight(node=108, fill="#ADD8E6", alpha=.5) + # Light Blue
  geom_hilight(node=151, fill="#98FB98", alpha=.5) # Light Green

    # Adding labels and points with clear separation and appropriate legend settings
b_SAMD9_9L <- b + 
  geom_tiplab(size= 4, offset = 12, hjust = 0.15) + 
  geom_tippoint(aes(color = SAMD9_gene_status, size = 0.3), position = position_nudge(x = 3)) + 
  geom_tippoint(aes(color = SAMD9L_gene_status, size = 0.3), position = position_nudge(x = 6)) + 
  scale_size_continuous(range = c(1, 6)) +
  scale_color_manual(values = c(I = "green", L = "red", M = "blue", PG = "lightblue", UL = "orange", No_chain = "grey")) +
  labs(color = "Gene status (TOGA classification)") +  # Set the legend title for color
  theme(
    legend.position = c(0.9, 0.97),  # Position legend at the top right
    legend.background = element_rect(fill = "white", color = "black"),  # Add background to the legend for clarity
    legend.key = element_rect(fill = "white")  # Clear background for legend keys
  )
b_SAMD9_9L
ggsave("SAMD9_9L.png", units="in", width=18, height=18, dpi=800)
b_SAMD9_9L
dev.off()
  