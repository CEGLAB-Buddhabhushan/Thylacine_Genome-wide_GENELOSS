getwd()
library(devtools)
library(ggside)
library(tidyverse)
library(tidyquant)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gghalves)

df <- read.table("TOGA_info.nonempty.Tasmanian_devil_paralogue_mart_export.tsv", header = TRUE, sep = '\t')
df<-na.omit(df)


facet_GC<-ggplot(df, aes(
  x = Paralogue_.id._target_Tasmanian_devil_gene_identical_to_query_gene,
  y = Gene_._GC_content,
  color = Paralogue_last_common_ancestor_with_Tasmanian_devil
)) +
  geom_point(alpha = 0.6, size = 1.5) +
  facet_wrap(~ TOGA_status) +
  theme_classic2() +
  labs(
    x = "% Identity (Target Tasmanian devil gene vs Paralog)",
    y = "GC Content (%)",
  )+
  theme(legend.position = "bottom")

png("GC_content_Vs_per.identity.png", units="in", width=16, height=9, res=900)
facet_GC
dev.off()

########## rain plot ########
data <- read.table("./TOGA_info.nonempty.Tasmanian_devil_paralogue_mart_export.tsv", header = T, sep = "\t")
data<-na.omit(data)
data[1,]
mean_df <- data %>%
  group_by(TOGA_status) %>%
  summarise(mean_iden = mean(`Paralogue_.id._target_Tasmanian_devil_gene_identical_to_query_gene`, na.rm = TRUE))
pairwise.wilcox.test(data$Paralogue_.id._target_Tasmanian_devil_gene_identical_to_query_gene, data$TOGA_status, p.adjust.method = "BH")
count_df <- data %>%
  group_by(TOGA_status) %>%
  summarise(n = n(),
            mean_iden = mean(`Paralogue_.id._target_Tasmanian_devil_gene_identical_to_query_gene`, na.rm = TRUE))

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



# Basic raincloud plot
paralog_iden<-ggplot(data, aes(x = TOGA_status, y = Paralogue_.id._target_Tasmanian_devil_gene_identical_to_query_gene, fill = TOGA_status)) +
  
  # Half-violin plot (distribution)
  geom_half_violin(
    aes(fill = TOGA_status),
    side = "l", 
    alpha = 0.6,
    color = NA,
    trim = FALSE
  ) +
  geom_text(data = count_df, aes(x = TOGA_status, y = 95, label = paste0("n=", n)),
            size = 3.5, color = "black")+
  
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
  geom_text(data = mean_df, aes(x = TOGA_status, y = mean_iden, 
                                label = sprintf("%.2f", mean_iden)),
            vjust = 1,hjust =1.5,  size = 4, color = "black") +
  labs(x = "TOGA Status",
       y = "% Identity (Target Tasmanian devil gene vs Paralog)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


data <- read.table("S.harrisii.info.plot.tsv", header = T, sep = "\t")
count_df <- data %>%
  group_by(TOGA_status) %>%
  summarise(n = n(),
            mean_iden = mean(GC_Content, na.rm = TRUE))
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



# Basic raincloud plot
GC_content<-ggplot(data, aes(x = TOGA_status, y = GC_Content, fill = TOGA_status)) +
  
  # Half-violin plot (distribution)
  geom_half_violin(
    aes(fill = TOGA_status),
    side = "l", 
    alpha = 0.6,
    color = NA,
    trim = FALSE
  ) +
  geom_text(data = count_df, aes(x = TOGA_status, y = 95, label = paste0("n=", n)),
            size = 3.5, color = "black")+
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


p1 <- GC_content +
  rremove("x.text") +
  rremove("xlab") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal" )+
      guides(fill = guide_legend(nrow = 1))
    
p2 <- paralog_iden +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  guides(fill = guide_legend(nrow = 1))

# Arrange with common one-row legend
png("GC_content_Vs_per.identity.raincloud_plot_with_pairwise_wilcox.png", units="in", width=16, height=9, res=900)
ggarrange(
  p1,
  p2,
  align = "v",
  ncol = 1,
  nrow = 2,
  common.legend = TRUE,
  legend = "bottom"
)
dev.off()


############# Gene rank plot #################
getwd()
library(ggplot2)
library(dplyr)
library(tidyverse)
library(forcats)
library(ggpubr)


data <- read.table("./TOGA_info.nonempty.Tasmanian_devil_paralogue_mart_export.tsv", header = T, sep = "\t")
data<-na.omit(data)
df<- data[,c(1,4,6,8,10,15,16,17,19)]
head(df)
df_unique <- df[!duplicated(df), ] 
head(df_unique)
##############
gene_gc_ranked <- df_unique %>%
  select(Gene_stable_ID, gene_gc = Gene_._GC_content) %>%
  distinct() %>%
  arrange(gene_gc) %>%
  mutate(Gene_stable_ID_rank = row_number())
length(unique(gene_gc_ranked$Gene_stable_ID_rank))
head(gene_gc_ranked)
# Step 2: Rank unique paralogue GC content
paralog_gc_ranked <- df_unique %>%
  select(Tasmanian_devil_paralogue_gene_stable_ID, gene_gc = Gene_._GC_content) %>%
  distinct() %>%
  arrange(gene_gc) %>%
  mutate(Tasmanian_devil_paralogue_gene_stable_ID_rank = row_number())
length(unique(paralog_gc_ranked$Tasmanian_devil_paralogue_gene_stable_ID_rank))
head(paralog_gc_ranked)



df_plot <- df_unique %>%
  left_join(gene_gc_ranked, by = "Gene_stable_ID") %>%
  left_join(paralog_gc_ranked, by = "Tasmanian_devil_paralogue_gene_stable_ID")

head(df_plot)
View(df_plot)
##############
df_plot <- df_plot %>%
  mutate(
    identity_bin_raw = cut(
      Paralogue_.id._target_Tasmanian_devil_gene_identical_to_query_gene,
      breaks = c(0, 25, 50, 75, 100),
      include.lowest = TRUE,
      right = FALSE  # [0–25), [25–50), etc.
    )
  )

# Step 2: Create readable labels
df_plot <- df_plot %>%
  mutate(
    identity_bin_labeled = recode_factor(identity_bin_raw,
                                         "[0,25)" = "0–25%",
                                         "[25,50)" = "25–50%",
                                         "[50,75)" = "50–75%",
                                         "[75,100]" = "75–100%"
    )
  )

# Step 3: Reorder levels (for legend ordering)
df_plot <- df_plot %>%
  mutate(
    identity_bin = fct_relevel(identity_bin_labeled, "0–25%", "25–50%", "50–75%", "75–100%")
  )

View(df_plot)
head(df_plot)
ggplot(df_plot, aes(x = Gene_stable_ID_rank, y = Tasmanian_devil_paralogue_gene_stable_ID_rank)) +
  geom_point(aes(
    color = TOGA_status,
    size = identity_bin
  ), alpha = 0.3) +
  scale_color_manual(values = c("skyblue", "blue", "green", "purple", "orange", "pink", "yellow")) +  # Manual color scale
  labs(
    x = "Gene stable ID (rank by GC content)",
    y = "Tasmanian devil paralogue transcript ID (rank by GC content)",
    color = "TOGA status",
    size = "% Identity (Target Tasmanian devil gene vs Paralog)"
  ) +
  theme_bw()

gene_rank<-ggplot(df_plot, aes(x = Gene_stable_ID_rank, y = Tasmanian_devil_paralogue_gene_stable_ID_rank)) +
  geom_point(aes(
    color = TOGA_status,
    size = identity_bin,
  ), alpha = 0.3) +
  labs(
    x = "Gene stable ID (rank by GC content)",
    y = "Tasmanian devil paralogue transcript ID (rank by GC content)",
    color = "TOGA status",
    size = "% Identity (Target Tasmanian devil gene vs Paralog)"
  ) +
  facet_wrap(~ TOGA_status) +
  theme_classic2()+
  theme(legend.position = "bottom")


##### for all####
all_label_df <- df_plot %>%
  group_by(Gene_name) %>%
  slice_max(order_by = Paralogue_.id._target_Tasmanian_devil_gene_identical_to_query_gene, n = 1) %>%
  ungroup()
library(ggplot2)
library(ggpubr)

All_with_label<-ggplot(df_plot, aes(x = Gene_stable_ID_rank, y = Tasmanian_devil_paralogue_gene_stable_ID_rank)) +
  geom_point(aes(
    color = TOGA_status,
    shape = identity_bin
  ), alpha = 0.3) +
  geom_text(
    data = all_label_df,
    aes(label = Gene_name),
    size = 2.5,hjust=1.5,
    check_overlap = TRUE  # avoids overlapping labels
  ) +
  labs(
    x = "Gene stable ID (rank by GC content)",
    y = "Tasmanian devil paralogue transcript ID (rank by GC content)",
    color = "TOGA status",
    size = "% Identity (Target Tasmanian devil gene vs Paralog)"
  ) +
  theme_classic2() +
  facet_wrap(~ TOGA_status) +
  theme(legend.position = "bottom")

png("All_class.gene_rank.plot.png", units="in", width=16, height=9, res=900)
All_with_label
dev.off()



### with shapes
All_with_label <- ggplot(df_plot, aes(x = Gene_stable_ID_rank, y = Tasmanian_devil_paralogue_gene_stable_ID_rank)) +
  geom_point(aes(
    shape = identity_bin
  ), alpha = 0.6, size = 2.5) +
  geom_text(
    data = all_label_df,
    aes(label = Gene_name),
    size = 2.5,
    hjust = 1.5,
    check_overlap = TRUE
  ) +
  labs(
    x = "Gene stable ID (rank by GC content)",
    y = "Tasmanian devil paralogue transcript ID (rank by GC content)",
    shape = "% Identity Bin"
  ) +
  scale_shape_manual(values = c(4, 2, 8, 20)) +  # adjust as per number of identity_bin levels
  theme_classic2() +
  facet_wrap(~ TOGA_status) +
  theme(legend.position = "bottom")
All_with_label <- ggplot(df_plot, aes(x = Gene_stable_ID_rank, y = Tasmanian_devil_paralogue_gene_stable_ID_rank)) +
  geom_point(aes(
    shape = identity_bin,
    fill = identity_bin  # use fill instead of color
  ), color = "black", alpha = 0.6, size = 2.5) +
  geom_text(
    data = all_label_df,
    aes(label = Gene_name),
    size = 2.5,
    hjust = 1.5,
    check_overlap = TRUE
  ) +
  labs(
    x = "Gene stable ID (rank by GC content)",
    y = "Tasmanian devil paralogue gene stable ID (rank by GC content)",
    shape = "% Identity Bin",
    fill = "% Identity Bin"
  ) +
  scale_shape_manual(values = c(4, 1, 21, 24)) +  # fillable shapes
  theme_classic2() +
  facet_wrap(~ TOGA_status) +
  theme(legend.position = "bottom")
All_with_label

library(ggplot2)
library(dplyr)

# Get unique TOGA_status values
toga_levels <- unique(df_plot$TOGA_status)

# Loop through each TOGA_status and save individual facet plots
for (status in toga_levels) {
  df_subset <- df_plot %>% filter(TOGA_status == status)
  label_subset <- all_label_df %>% filter(TOGA_status == status)
  
  p <- ggplot(df_subset, aes(x = Gene_stable_ID_rank, y = Tasmanian_devil_paralogue_gene_stable_ID_rank)) +
    geom_point(aes(
      shape = identity_bin,
      fill = identity_bin
    ), color = "black", alpha = 0.6, size = 2.5) +
    geom_text(
      data = label_subset,
      aes(label = Gene_name),
      size = 2.5,
      hjust = 1.5,
      check_overlap = TRUE
    ) +
    labs(
      x = "Gene stable ID (rank by GC content)",
      y = "Tasmanian devil paralogue gene stable ID (rank by GC content)",
      shape = "% Identity Bin",
      fill = "% Identity Bin",
      title = paste("TOGA Status:", status)
    ) +
    scale_shape_manual(values = c(4, 1, 21, 24)) +
    theme_classic2() +
    theme(legend.position = "bottom")
  
  # Save the plot
  ggsave(filename = paste0("facet_", gsub("[^A-Za-z0-9]", "_", status), ".png"),
         plot = p, width = 16, height = 9, units = "in", dpi = 300)
}

write.table(df_plot, file = "df_plot.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
