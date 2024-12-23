library(ggplot2)
library(reshape2)
library(gridExtra)

data <- read.csv("Devil-DUNNART.ortho.counts.csv", header = TRUE, row.names = 1)

plots <- list()

for (gene in rownames(data)) {
  # Subset data for the gene and reshape
  gene_data <- data[gene, , drop = FALSE]
  gene_melted <- data.frame(Sample = colnames(gene_data), TPM = as.numeric(gene_data))
  
# Create the heatmap
p <- ggplot(gene_melted, aes(x = Sample, y = gene, fill = TPM)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "yellow", high = "blue", name = "TPM") +
    labs(x = "Samples", y = "", title = gene) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, face = "bold"))
 
  plots[[gene]] <- p
}

png("Gene_expression.png", units="in", width=16, height=10.5, res=900)
grid.arrange(grobs = plots, ncol = 1) 
dev.off()
