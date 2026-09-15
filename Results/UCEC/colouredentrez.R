#install.packages("BiocManager")
#BiocManager::install("biomaRt")


# Load biomaRt
library(biomaRt)

setwd("/home/simone/Desktop/BIO_PROJ/Results/UCEC")

# Read up and downregulated gene symbol files
up_genes <- read.table("UP.txt", stringsAsFactors = FALSE)[,1]
down_genes <- read.table("DOWN.txt", stringsAsFactors = FALSE)[,1]

# Connect to Ensembl and choose human dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Function to convert symbols to Entrez IDs using biomaRt
convert_symbols_to_entrez <- function(gene_symbols) {
  getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
        filters = "hgnc_symbol",
        values = gene_symbols,
        mart = ensembl)
}

# Convert each list
up_entrez <- convert_symbols_to_entrez(up_genes)
down_entrez <- convert_symbols_to_entrez(down_genes)

# Add color labels
up_entrez$color <- "orange"
down_entrez$color <- "skyblue"

# Combine and select only ENTREZID and color
kegg_data <- rbind(
  up_entrez[, c("entrezgene_id", "color")],
  down_entrez[, c("entrezgene_id", "color")]
)

# Remove rows with missing IDs
kegg_data <- na.omit(kegg_data)

# Remove duplicates
kegg_data <- unique(kegg_data)

# Write KEGG color file (no header, no quotes)
write.table(kegg_data,
            file = "kegg_colors.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
