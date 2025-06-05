# Load required library
library(readr)
library(dplyr)

# Load your fully annotated BDP catalogue
bdp <- read_csv("BDP_Fully_Annotated.csv")

# Filter lncRNAs on each strand
watson_lncRNA <- bdp %>%
  filter(watson_biotype == "lncRNA") %>%
  dplyr::select(watson_ensembl_id) %>%
  distinct()

crick_lncRNA <- bdp %>%
  filter(crick_biotype == "lncRNA") %>%
  dplyr::select(crick_ensembl_id) %>%
  distinct()

# Combine and get unique lncRNA Ensembl IDs
all_lncRNAs <- bind_rows(
  watson_lncRNA %>% rename(ensembl_id = watson_ensembl_id),
  crick_lncRNA %>% rename(ensembl_id = crick_ensembl_id)
) %>%
  filter(!is.na(ensembl_id)) %>%
  distinct()

# Count total unique lncRNAs across both strands
total_unique_lncRNAs <- nrow(all_lncRNAs)

# Count number of BDPs with at least one lncRNA on either strand
bdp_with_lncRNA <- bdp %>%
  filter(watson_biotype == "lncRNA" | crick_biotype == "lncRNA")

bdp_with_lncRNA_count <- nrow(bdp_with_lncRNA)

# Output the results
cat("Total unique lncRNAs in BDPs:", total_unique_lncRNAs, "\n")
cat("Number of BDPs with at least one lncRNA:", bdp_with_lncRNA_count, "\n")
