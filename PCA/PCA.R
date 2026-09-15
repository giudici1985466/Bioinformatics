rm(list)
options(stringsAsFactors = F)
library(ggbiplot)
library(qcc)
library(ggpubr)
library(factoextra)
library(corrplot)
library(FactoMineR)
library(RColorBrewer)
library(survminer)
library(stringr)





setwd("/home/simone/Desktop/BIO_PROJ/")
path_in <- "/home/simone/Desktop/BIO_PROJ/Results/"


dataset <-"UCEC"


dirRes <- "PCA/Results/"

if (!dir.exists(dirRes)){
  dir.create(dirRes)
} else {
  print(paste("The directory", dirRes, "already exists"))
}




dirDataset <- paste0(dirRes,dataset,"/")  #change it to make it have sense

if (!dir.exists(dirDataset)){
  dir.create(dirDataset)
} else {
  print(paste("The directory", dirDataset, "already exists"))
}




################################################

filename_in <- paste0(path_in,dataset,"/matrix_DEG.txt")
filename_list_normal <- paste0(path_in,dataset,"/normal.txt")
filename_list_tumor <- paste0(path_in,dataset,"/tumor.txt")

file_score_plot <- paste0(dirDataset,"score_plot.pdf")
file_pareto_scree_plot <- paste0(dirDataset,"pareto_scree_plots.pdf")
file_loading_plot <- paste0(dirDataset, "loadings_plot.pdf")
file_contribution_plot <- paste0(dirDataset,"PC_contributions_plot.pdf")


##########################################
# 1. Importing data 


data <- read.table(filename_in, header = T, sep = "\t", quote = "", check.names = F, row.names = 1)
name_DEG <- data.frame(str_split_fixed(rownames(data),"\\|",2))

# Save gene codes and names
gene_codes <- name_DEG$X1
gene_names <- name_DEG$X2
gene_labels <- paste(gene_names, gene_codes, sep = " | ")

# Assign to data
rownames(data) <- gene_labels


list_normal <- read.table (filename_list_normal, header = F, quote ="", check.names = F)
list_normal <- list_normal$V1

list_tumor <- read.table(filename_list_tumor, header = F, quote = "", check.names = F)
list_tumor <- list_tumor$V1

data <- t(data[,c(list_tumor,list_normal)])     #transposing for the analysis 

groups <- c(rep("tumor", length(list_tumor)),
            rep("normal", length(list_normal)))


################################################
# 2. Apply PCA
# Rows of data correspond to observations (samples), columns to variables (genes)

pca <- prcomp(data, center = T, scale. = T, retx = T) 
################################################
# 3. Compute score and score plot
# (scores = the coordinates of old data (observations) in the new systems, that are the PCs)

# pca$x = t(data)*pca$rotation 
score <- pca$x 

# alternative for score computation
# scores <- get_pca_ind(pca)$coord
# colnames(scores) <- paste0('PC', seq(1,ncol(scores)))

pdf(file_score_plot,width=5, height=5) 
g <- ggbiplot(pca, 
              obs.scale = 1, 
              var.axes = F, 
              ellipse = T, 
              groups = groups)
print(g)
dev.off() 


################################################
# 4. Compute eigenvalue
# eigenvalues of the covariance matrix ordered in decreasing order (from the largest to the smallest)
eigenvalue = pca$sdev^2 

# variance explained by each PC
varS <- round(eigenvalue/sum(eigenvalue)*100, 2)
names(varS) = paste0("PC", seq(1, length(varS)))

df <- get_eigenvalue(pca)
rownames(df) <- paste0("PC", seq(1,nrow(df)))

# pareto chart
pareto.chart(varS[1:10])

# scree plot
fviz_eig(pca,addlabels = TRUE)

dev.off()

# Pareto chart + scree plot in same PDF
pdf(file_pareto_scree_plot, width = 7, height = 5)

# Pareto chart
pareto.chart(varS[1:10],
             ylab = "Percentage of Variance",
             col = rev(colorRampPalette(brewer.pal(9, "BuPu"))(10)))
             

# Scree plot
fviz_eig(pca, addlabels = TRUE)

dev.off()


loadings <- pca$rotation


contrib_var <- get_pca_var(pca)$contrib 
colnames(contrib_var) <- paste0('PC', seq(1,ncol(contrib_var)))

contrib_var <- contrib_var[order(contrib_var[,"PC1"], decreasing = T),]
pdf(file_contribution_plot, width= 10, height = 10)


corrplot(contrib_var[1:100, 1:10], is.corr=FALSE, 
         tl.col = "black", 
         method = "circle",
         col = brewer.pal(n = 9, name = "BuPu"),
         addCoef.col = NULL)

# to PC1
fviz_contrib(pca, choice = "var", axes = 1, top = 10)
# to PC2
fviz_contrib(pca, choice = "var", axes = 2, top = 10)
#SUM
fviz_contrib(pca, choice = "var", axes = 1:2, top = 10)


dev.off()




