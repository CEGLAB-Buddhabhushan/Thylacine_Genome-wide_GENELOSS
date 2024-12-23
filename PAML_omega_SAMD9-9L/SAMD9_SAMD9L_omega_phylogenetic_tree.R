setwd("./SAMD9-9L/")

library(phytools)

tree <- read.tree("Species_no_internode.nwk")
omega_data <- read.table("omega.tsv", header = TRUE, sep = "\t")

# Extract the omega values for SAMD9 and SAMD9L
omega_SAMD9 <- setNames(omega_data$ωSAMD9, omega_data$Species_name)
omega_SAMD9L <- setNames(omega_data$ωSAMD9L, omega_data$Species_name)

# Match species names in the tree with omega values
tree$tip.label <- gsub(" ", "_", tree$tip.label)  # Ensure compatibility
if (!all(tree$tip.label %in% names(omega_SAMD9))) {
  stop("Mismatch between tree labels and omega data.")
}

# Set consistent color gradient and range
custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)  # Blue to red gradient

# Map omega values onto the tree with fixed color scale
SAMD9_map <- contMap(tree, omega_SAMD9, plot = FALSE)
SAMD9L_map <- contMap(tree, omega_SAMD9L, plot = FALSE)

# Apply consistent color map and scale
SAMD9_map <- setMap(SAMD9_map, colors = custom_colors)
SAMD9L_map <- setMap(SAMD9L_map, colors = custom_colors)

png("SAMD9_SAMD9L_omega_phylogenetic_tree.png", width=16, height=9, units="in",  res = 900)
# Set up layout: 3 panels (left tree, text, right tree)
layout(matrix(1:3, 1, 3), widths = c(0.43, 0.14, 0.43))

# Plot tree with SAMD9 values
plot(SAMD9_map, lwd = 7, ftype = "off", legend = 60, outline = TRUE, fsize = c(0, 1.2), leg.txt = "ω (SAMD9)")
tiplabels(round(omega_SAMD9[tree$tip.label], 2), frame = "none", adj = 0.5, bg = "white", cex = 1.5)

# Create middle panel with species names
plot.new()
plot.window(xlim = c(-0.1, 0.1), ylim = get("last_plot.phylo", envir = .PlotPhyloEnv)$y.lim)
par(cex = 0.8)
text(rep(0, length(tree$tip.label)), 1:Ntip(tree), gsub("_", " ", tree$tip.label), font = 3)

# Plot tree with SAMD9L values facing left
plot(SAMD9L_map, lwd = 7, outline = TRUE, direction = "leftwards", ftype = "off", legend = 60, 
     fsize = c(0, 1.2), leg.txt = "ω (SAMD9L)")
tiplabels(round(omega_SAMD9L[tree$tip.label], 2), frame = "none", adj = +0.5, bg = "white", cex = 1.5)
dev.off()

