# https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html

library(enrichR)
library(ggplot2)
library(forcats)
library(stringr)
################################################
setwd("/home/simone/Desktop/BIO_PROJ")
source("enrichR/getEnrichment.R")
source("enrichR/getEnrichmentPlot.R")
################################################
dataset <- "UCEC"

dirRes <- "Results/"

if (!dir.exists(dirRes)){
  dir.create(dirRes)
}else{
  print(paste("The directory",dirRes,"already exists"))
}

dirDataset <- paste0(dirRes,dataset,"/")

if (!dir.exists(dirDataset)){
  dir.create(dirDataset)
}else{
  print(paste("The directory",dirDataset,"already exists"))
}

dirEnrich <- paste0(dirDataset,"Functional_Enrichment/")

if (!dir.exists(dirEnrich)){
  dir.create(dirEnrich)
}else{
  print(paste("The directory",dirEnrich,"already exists"))
}


################################################
top_term <- 20
thr_pval <- 0.05
################################################
file_input_list <- paste0(dirDataset,"DEG.txt")

dbs <- listEnrichrDbs() #list of all db available in enrichR   (As in the website)

dbs <- c("DisGeNET","GO_Molecular_Function_2021", "GO_Biological_Process_2021", "KEGG_2021_Human", "TRANSFAC_and_JASPAR_PWMs")   #selecting some of them that are relevant

input_list <- read.table(file_input_list, sep = "\t", header = T, check.names = F, quote = "")
# input_list <- input_list$genes
list <- split(input_list$genes, input_list$direction)


df <- lapply(list, function(x){
  enrichr(x, dbs)
})

getEnrichment(df$UP,"UP")
getEnrichment(df$DOWN,"DOWN")


