library(tidyverse)

setwd("/home/simone/Desktop/BIO_PROJ")


# Read the data
metadata <- read.table("clinical_ucec.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Check the structure
str(metadata$age_at_index)


# Convert to numeric if needed
metadata$age_at_index <- as.numeric(metadata$age_at_index)

# Remove NAs
metadata <- metadata %>% filter(!is.na(age_at_index))


metadata <- metadata %>%
  mutate(age_group = cut(age_at_index, 
                         breaks = seq(20, 100, by = 10),
                         labels = paste(seq(20, 90, by = 10), seq(29, 99, by = 10), sep = "-"),
                         right = FALSE))

ggplot(metadata, aes(x = age_group)) +
  geom_bar(fill = "pink3") +
  labs(title = "Age Group Distribution", x = "Age Group", y = "Count") +
  theme_minimal()