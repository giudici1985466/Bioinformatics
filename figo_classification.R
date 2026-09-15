# Load required libraries
library(tidyverse)

setwd("/home/simone/Desktop/BIO_PROJ")

# Read the metadata file
metadata <- read.table("clinical_ucec.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(metadata) <- c("submitter_id", "figo_stage")  # Adjust if needed

# Read the matrix file
matrix <- read.table("Results/UCEC/matrix_DEG.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
matrix_ids <- colnames(matrix)
matrix_truncated_ids <- sapply(strsplit(matrix_ids, "\\."), function(x) paste(x[1:3], collapse = "-"))
matrix_truncated_ids <- gsub("\\.", "-", matrix_truncated_ids)

metadata <- metadata[!is.na(metadata$figo_stage), ]
metadata <- metadata[!is.na(metadata$submitter_id), ]
metadata <- metadata[, 1:2]
metadata <- metadata[!duplicated(metadata$submitter_id), ]

# Normalize case, remove whitespace, and strip "Stage" prefix
metadata$figo_stage <- toupper(trimws(metadata$figo_stage))
metadata$figo_stage <- gsub("^STAGE\\s+", "", metadata$figo_stage)

# Collapse sub-stages (e.g., IIIC1, IIIA → III)
metadata$figo_stage <- gsub("^([IVX]+)[A-Z0-9]*$", "\\1", metadata$figo_stage)


# Mark whether each patient is present in the matrix
metadata <- metadata %>%
  mutate(InMatrix = ifelse(submitter_id %in% matrix_truncated_ids, "Present", "Absent"))

metadata_present <- metadata %>% filter(InMatrix == "Present")

figo_counts <- metadata_present %>%
  group_by(figo_stage) %>%
  summarise(Count = n(), .groups = "drop")

ggplot(figo_counts, aes(x = "", y = Count, fill = figo_stage)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "FIGO Stage Distribution", fill = "FIGO Stage") +
  theme_void() +
  theme(legend.position = "right")+
  scale_fill_brewer(palette = "Set3")

