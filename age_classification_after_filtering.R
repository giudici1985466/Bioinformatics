# Load required libraries
library(tidyverse)

# Set working directory (change this if needed)
setwd("/home/simone/Desktop/BIO_PROJ")

# Read the metadata file
metadata <- read.table("clinical_ucec.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Read the DEG matrix
matrix <- read.table("Results/UCEC/matrix_DEG.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

# Extract submitter IDs from matrix column names
matrix_ids <- colnames(matrix)
matrix_truncated_ids <- sapply(strsplit(matrix_ids, "\\."), function(x) paste(x[1:3], collapse = "-"))
matrix_truncated_ids <- gsub("\\.", "-", matrix_truncated_ids)

# Keep only necessary columns and clean metadata
# Select columns: submitter_id, figo_stage, age_at_index
metadata <- metadata[, c("submitter_id", "figo_stage", "age_at_index")]

# Clean NA and duplicates
metadata <- metadata[!is.na(metadata$figo_stage), ]
metadata <- metadata[!is.na(metadata$submitter_id), ]
metadata <- metadata[!duplicated(metadata$submitter_id), ]

# Normalize and clean FIGO stage
metadata$figo_stage <- toupper(trimws(metadata$figo_stage))
metadata$figo_stage <- gsub("^STAGE\\s+", "", metadata$figo_stage)
metadata$figo_stage <- gsub("^([IVX]+)[A-Z0-9]*$", "\\1", metadata$figo_stage)

# Filter only patients present in DEG matrix
metadata <- metadata %>%
  mutate(InMatrix = ifelse(submitter_id %in% matrix_truncated_ids, "Present", "Absent"))

metadata_present <- metadata %>%
  filter(InMatrix == "Present")

# Clean and convert age
metadata_present$age_at_index <- as.numeric(metadata_present$age_at_index)
metadata_present <- metadata_present %>% filter(!is.na(age_at_index))

# Create age bins (10-year intervals)
metadata_present <- metadata_present %>%
  mutate(age_group = cut(age_at_index,
                         breaks = seq(20, 100, by = 10),
                         labels = paste(seq(20, 90, by = 10), seq(29, 99, by = 10), sep = "-"),
                         right = FALSE))

# Plot the age distribution
ggplot(metadata_present, aes(x = age_group)) +
  geom_bar(fill = "pink3") +
  labs(title = "Age Distribution of Patients (filtered)",
       x = "Age Group", y = "Count") +
  theme_minimal()
