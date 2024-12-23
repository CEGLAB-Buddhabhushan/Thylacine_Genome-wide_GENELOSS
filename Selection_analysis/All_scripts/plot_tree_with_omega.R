#!/usr/bin/env Rscript
library(phytools)
library(ape)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_tree_with_omega.R <tree_file> <data_file>")
}

tree_file <- args[1]  
data_file <- args[2]  

tree <- read.tree(tree_file)

omega_data <- read.delim(data_file, header = TRUE, sep = "\t")

omega_data <- omega_data[omega_data$Species_name %in% tree$tip.label, ]

omega <- setNames(omega_data$omega, omega_data$Species_name)

if (!all(tree$tip.label %in% names(omega))) {
  missing <- setdiff(tree$tip.label, names(omega))
  stop("Missing omega values for species: ", paste(missing, collapse = ", "))
}


cont_map <- phytools::contMap(tree, omega, outline = FALSE, lwd = 7)


output_file <- paste0(tools::file_path_sans_ext(tree_file), "_with_omega.png")

png(output_file, width = 800, height = 800)
par(mar = c(5, 4, 4, 2) + 0.5)
plot(cont_map, legend = 0.7 * max(nodeHeights(tree)), fsize = 0.8, leg.txt = "ω (dN/dS)", leg.cex = 1)
title("Phylogenetic Tree with ω (dN/dS) Values")
dev.off()

