library(dplyr)

# Read the data from the TSV file
df <- read.csv("HYPHY_RELAX.Results.p_val.adj.tsv", sep="\t", header=TRUE)

# Calculate -log10 of the adj.p_value column
df <- df %>%
  mutate(log_adj_p_value = -log10(adj.p_value))

# Save the result to a new TSV file
write.table(df, "HYPHY_RELAX.Results.p_val_with_log10.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# Print a message
cat("New file with -log10(adj.p_value) saved as 'HYPHY_RELAX.Results.p_val_with_log10.tsv'")
