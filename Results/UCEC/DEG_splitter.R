# Set your working directory and input file
setwd("/home/simone/Desktop/BIO_PROJ/Results/UCEC/")
input_file <- "DEG.txt" 

# Read the DEG table
deg_data <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Split into UP and DOWN regulated genes
up_genes <- deg_data[deg_data$direction == "UP", ]
down_genes <- deg_data[deg_data$direction == "DOWN", ]

# Write gene names only to individual files (one column)
write.table(up_genes$genes, file = "UP.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(down_genes$genes, file = "DOWN.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Optional: Write full rows if you prefer full metadata
# write.table(up_genes, file = "UP_full.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(down_genes, file = "DOWN_full.txt", sep = "\t", quote = FALSE, row.names = FALSE)
