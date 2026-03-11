# 1. Load the human microglia proteome data (proteinGroups.txt)
human_proteome_data <- read.delim("proteinGroups.txt", sep = "\t")

# --- THE MISSING CRITICAL STEP: Canonical Filtering [cite: 1681-1682] ---
# This removes isoforms (marked with -\d) to ensure data high-fidelity.
non_canonical_filter <- str_detect(human_proteome_data$Protein.IDs, pattern = "-\\d")
canonical_proteome_data <- human_proteome_data[!non_canonical_filter, ] 
# ------------------------------------------------------------------------

# 2. Load and clean your specific LSD gene list [cite: 1929]
unique_proteins <- read.csv("Protein_list.csv")$x
your_genes <- toupper(str_trim(str_remove_all(unique_proteins, ",.*")))

# 3. Perform the mapping on the CANONICAL data [cite: 1688-1690]
mapped_data <- canonical_proteome_data %>%
  mutate(Gene_List = str_split(toupper(Gene.names), ";")) %>%
  unnest(Gene_List) %>%
  filter(Gene_List %in% your_genes) %>%
  distinct(Gene_List, .keep_all = TRUE)
