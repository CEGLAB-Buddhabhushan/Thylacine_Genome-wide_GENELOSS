setwd("/home/ceglab358/BUDDHA/upload-GITHUB_Thylacine_Genome-wide_GENELOSS/Phylogenetic_correlation/SAMD9_sorted/")
getwd()

# Read the dataset
df <- read.table("final_output.csv", header = TRUE, sep = ",", check.names = FALSE)
colnames(df)
################## Phylogenetic signal ####################
library(phylolm)
library(ape)
View(df)

df$SAMD9_gene_status <- as.character(df$SAMD9_gene_status)

# Replace "I" with 1 and "L" with 0
df$SAMD9_gene_status[df$SAMD9_gene_status == "I"] <- 1
df$SAMD9_gene_status[df$SAMD9_gene_status == "L"] <- 0

# Convert to numeric if needed
df$SAMD9_gene_status <- as.numeric(df$SAMD9_gene_status)
#View(df)
colnames(df)
tree=read.tree("Species.nwk")
rownames(df)=df[,1]
head(df)
colnames(df)
str(df)
df <- df[!(rownames(df) %in% c("Canis_familiaris", "Spalax_galili", "Spermophilus_tridecemlineatus")), ]

# Identify tips with zero branch length
zero_length_tips <- tree$tip.label[
  tree$edge.length[match(1:length(tree$tip.label), tree$edge[,2])] == 0
]

# Drop those tips from the tree
tree <- drop.tip(tree, zero_length_tips)

# Also remove them from the df
df <- df[!rownames(df) %in% zero_length_tips, ]
keep_species <- rownames(df)

# Drop tips that are NOT in your data
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, keep_species))

# Confirm the tree has only your focal species
num_tips<-length(pruned_tree$tip.label) 
head(df,4)


#num_tips <- length(tree$tip.label)
#num_internal_nodes <- Nnode(tree)
#total_nodes <- num_tips + num_internal_nodes
#total_nodes
# Extract row names and tree tip labels
data_species <- rownames(df)
tree_species <- pruned_tree$tip.label
num_tips<- length(pruned_tree$tip.label)
# Check for mismatches
setdiff(data_species, tree_species)  # Species in df but not in tree
setdiff(tree_species, data_species)  # Species in tree but not in df
df$carnivore <- df$Inv + df$Vend + df$Vect + df$Vfish + df$Vunk + df$Scav
df$non_carnivore <- df$Fruit + df$Nect + df$Seed + df$PlantO
View(df)
condiation<-colnames(df)[4:13]
+Fruit+Nect+Seed+PlantO 
###############MPLE###############
fitwmple=phyloglm(SAMD9_gene_status ~ Vend ,df,pruned_tree, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4,start.beta=NULL, start.alpha=NULL,boot = 2000, full.matrix = TRUE)
cc1=coef(fitwmple)
summary(fitwmple)
mpleinter=round(cc1[1],2)
mplew=round(cc1[2],2)
summary_fit_mple <- summary(fitwmple)
mplepvals <- summary_fit_mple$coefficients[, "p.value"][2]
#####################IG10#########
fitwig10=phyloglm(SAMD9_gene_status ~ Vend ,df,pruned_tree, method = c("logistic_IG10"), btol = 20, log.alpha.bound = 4,start.beta=NULL, start.alpha=NULL,boot = 2000, full.matrix = TRUE)
cc2=coef(fitwig10)
summary(fitwig10)
ig10inter=round(cc2[1],2)
ig10w=round(cc2[2],2)
summary_fit_ig10 <- summary(fitwig10)
ig10pvals <- summary_fit_ig10$coefficients[, "p.value"][2]
####################################
t(table(df$SAMD9_gene_status,df$Vend))->focal
first_col_vector <- as.numeric(rownames(focal))
x_vals <- rep(first_col_vector, 2) 
y_vals <- c(rep(1, length(first_col_vector)), rep(0, length(first_col_vector)))
cex_vals <- c(as.vector(focal[,2]), as.vector(focal[,1])) / 3
##plot
png("Phylogenetic_signal_of_diet.png",height=9,width=10,units = "in",res=600)
plot(x_vals, y_vals, cex=cex_vals, pch=19, xlab="", ylab="", xlim=c(0,100), axes=FALSE)
axis(1)
axis(2, at = seq(0,1,by=1), las = 1)
box()
mtext("Loss",side=2,line=1.5,at=0)
mtext(bquote("The proportion of endothermic vertebrates (%)"),
      side = 1, line = 2, at = 50)
mtext(substitute(paste(italic(SAMD9), " gene status")),side=2,line=2,at=0.5)
mtext("Intact",side=2,line=1.5,at=1)
curve(plogis(cc1[1]+cc1[2]*x),col="red",add=TRUE,lwd=2)
textleg=paste("logistic_MPLE\nIntercept = ",mpleinter,"\nSlope = ",mplew,"\np = ",mplepvals,"\nn = ",num_tips)
legend(85,0.8,textleg,xjust = 0.5,yjust = 0.5,x.intersp = 0.2,y.intersp = 0.8,adj = c(0, 0.5),bty='n', text.col = "red")
curve(plogis(cc2[1]+cc2[2]*x),col="blue",add=TRUE,lwd=2)
text2leg=paste("logistic_IG10\nIntercept = ",ig10inter,"\nSlope = ",ig10w,"\np = ",ig10pvals,"\nn = ",num_tips)
legend(40,0.4,text2leg,xjust = 0.5,yjust = 0.5,x.intersp = 0.2,y.intersp = 0.8,adj = c(0, 0.5),bty='n', text.col = "blue")
dev.off()

##### Data from both datasets ####
setwd("/home/ceglab358/BUDDHA/upload-GITHUB_Thylacine_Genome-wide_GENELOSS/Phylogenetic_correlation/")
getwd()

# Read the dataset
df <- read.table("Diet_data_from_Jiao, Hengwu et al_and_Kapsetaki, Stefania E et al.tsv", header = TRUE, sep = "\t", check.names = FALSE)
colnames(df)
################## Phylogenetic signal ####################
library(phylolm)
library(ape)
View(df)

df$SAMD9_gene_status <- as.character(df$SAMD9_gene_status)

# Replace "I" with 1 and "L" with 0
df$SAMD9_gene_status[df$SAMD9_gene_status == "I"] <- 1
df$SAMD9_gene_status[df$SAMD9_gene_status == "L"] <- 0

# Convert to numeric if needed
df$SAMD9_gene_status <- as.numeric(df$SAMD9_gene_status)
#View(df)
colnames(df)
tree=read.tree("Species.nwk")
rownames(df)=df[,1]
head(df)
colnames(df)
str(df)
df <- df[!(rownames(df) %in% c("Canis_familiaris", "Spalax_galili", "Spermophilus_tridecemlineatus")), ]

# Identify tips with zero branch length
zero_length_tips <- tree$tip.label[
  tree$edge.length[match(1:length(tree$tip.label), tree$edge[,2])] == 0
]

# Drop those tips from the tree
tree <- drop.tip(tree, zero_length_tips)

# Also remove them from the df
df <- df[!rownames(df) %in% zero_length_tips, ]
keep_species <- rownames(df)

# Drop tips that are NOT in your data
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, keep_species))

# Confirm the tree has only your focal species
num_tips<-length(pruned_tree$tip.label) 
head(df,4)


#num_tips <- length(tree$tip.label)
#num_internal_nodes <- Nnode(tree)
#total_nodes <- num_tips + num_internal_nodes
#total_nodes
# Extract row names and tree tip labels
data_species <- rownames(df)
tree_species <- pruned_tree$tip.label
num_tips<- length(pruned_tree$tip.label)
# Check for mismatches
setdiff(data_species, tree_species)  # Species in df but not in tree
setdiff(tree_species, data_species)  # Species in tree but not in df

View(df)
set.seed(35842)

###############MPLE###############
fitwmple=phyloglm(SAMD9_gene_status ~ Vend ,df,pruned_tree, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4,start.beta=NULL, start.alpha=NULL,boot = 2000, full.matrix = TRUE)
cc1=coef(fitwmple)
summary(fitwmple)
mpleinter=round(cc1[1],2)
mplew=round(cc1[2],2)
summary_fit_mple <- summary(fitwmple)
mplepvals <- summary_fit_mple$coefficients[, "p.value"][2]
#####################IG10#########
fitwig10=phyloglm(SAMD9_gene_status ~ Vend ,df,pruned_tree, method = c("logistic_IG10"), btol = 20, log.alpha.bound = 4,start.beta=NULL, start.alpha=NULL,boot = 2000, full.matrix = TRUE)
cc2=coef(fitwig10)
summary(fitwig10)
ig10inter=round(cc2[1],2)
ig10w=round(cc2[2],2)
summary_fit_ig10 <- summary(fitwig10)
ig10pvals <- summary_fit_ig10$coefficients[, "p.value"][2]
####################################
t(table(df$SAMD9_gene_status,df$Vend))->focal
first_col_vector <- as.numeric(rownames(focal))
x_vals <- rep(first_col_vector, 2) 
y_vals <- c(rep(1, length(first_col_vector)), rep(0, length(first_col_vector)))
cex_vals <- c(as.vector(focal[,2]), as.vector(focal[,1])) / 3
legend_species_vals <-  c(1,2,3,4,15)
legend_cex_vals <- legend_species_vals / 3
xpos <- 10  # x location
ypos_start <- 0.11  # y starting point
y_spacing <- 0.06  


##plot
png("Phylogenetic_signal_of_diet.png",height=9,width=10,units = "in",res=600)
plot(x_vals, y_vals, cex=cex_vals, pch=19, xlab="", ylab="", xlim=c(0,100), axes=FALSE)
axis(1)
axis(2, at = seq(0,1,by=1), las = 1)
box()
mtext("Loss",side=2,line=1.5,at=0)
mtext(bquote("Diet of vertebrate endotherms (%)"),
      side = 1, line = 2, at = 50)
mtext(substitute(paste(italic(SAMD9), " gene status")),side=2,line=2,at=0.5)
mtext("Intact",side=2,line=1.5,at=1)
curve(plogis(cc1[1]+cc1[2]*x),col="red",add=TRUE,lwd=3)
textleg=paste("logistic_MPLE\nIntercept = ",mpleinter,"\nSlope = ",mplew,"\np = ",mplepvals,"\nn = ",num_tips)
legend(85,0.7,textleg,xjust = 0.5,yjust = 0.5,x.intersp = 0.2,y.intersp = 0.8,adj = c(0, 0.5),bty='n', text.col = "red")
curve(plogis(cc2[1]+cc2[2]*x),col="blue",add=TRUE,lwd=3)
text2leg=paste("logistic_IG10\nIntercept = ",ig10inter,"\nSlope = ",ig10w,"\np = ",ig10pvals,"\nn = ",num_tips)
legend(75,0.3,text2leg,xjust = 0.5,yjust = 0.5,x.intersp = 0.2,y.intersp = 0.8,adj = c(0, 0.5),bty='n', text.col = "blue")
text(xpos-5, ypos_start + (length(legend_species_vals) + 1)*y_spacing, 
     "Dot size  Number of species", font=1, adj=0)
for (i in seq_along(legend_species_vals)) {
  points(xpos, ypos_start + (i * y_spacing), 
         pch=19, cex=legend_cex_vals[i])
  text(xpos + 10, ypos_start + (i * y_spacing), 
       paste0(legend_species_vals[i], ""), adj=0, cex=0.9)
}
dev.off()

