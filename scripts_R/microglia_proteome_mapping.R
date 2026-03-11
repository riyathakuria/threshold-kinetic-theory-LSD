# -------------------------------------------------------------------
# Microglia Proteome Mapping Workflow
# Part of the Threshold Kinetic Theory Project
# -------------------------------------------------------------------

library(dplyr)
library(stringr)
library(tidyr)

# 1. Load the human microglia proteome data (proteinGroups.txt)
# Source: Human Microglia Proteome [cite: 18-19]
human_proteome_data <- read.delim("proteinGroups.txt", sep = "\t")

# 2. Load and clean your specific LSD gene list [cite: 1929]
# Removes commas and whitespace to ensure matching [cite: 1662-1663]
unique_proteins <- read.csv("Protein_list.csv")$x
your_genes <- toupper(str_trim(str_remove_all(unique_proteins, ",.*")))

# 3. Perform the mapping [cite: 1665-1672]
# Splits multi-protein strings and finds overlaps with your list
mapped_data <- human_proteome_data %>%
  mutate(Gene_List = str_split(toupper(Gene.names), ";")) %>%
  unnest(Gene_List) %>%
  filter(Gene_List %in% your_genes) %>%
  distinct(Gene_List, .keep_all = TRUE)

# 4. Save the final mapped set (Targeting the 746-gene subset) [cite: 1775-1780]
write.csv(mapped_data, "Final_Clean_Microglia_Proteome_Mapping.csv", row.names = FALSE)
