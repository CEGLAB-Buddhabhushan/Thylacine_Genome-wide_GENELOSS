hyphy<- read.table("HYPHY_RELAX.Results.sorted.txt", sep = "\t", header  = T)
df <- data.frame(hyphy)
df$adj.p_value <- p.adjust(df$pval, method = "fdr")
write.table(df, file='HYPHY_RELAX.Results.p_val.adj.tsv', quote=T, sep='\t', col.names = NA)
