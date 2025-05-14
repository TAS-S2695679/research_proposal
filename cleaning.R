library(dplyr)
library(readr)
library(stringr)
library(biomaRt)
library(tidyverse)

# Step 1: Load the Cataloguew
catalogue <- read_csv("Catalogue.csv") %>% as_tibble()

# Step 2: Filter to transcript-level entries only
transcripts <- catalogue %>%
  filter(Feature == "transcript")

# Step 3: Extract base BDP ID by removing '_2' or '_3' suffix from Annotation
transcripts <- transcripts %>%
  mutate(BDP_ID = str_replace(Annotation, "_[23]$", ""))

# Step 4: Separate Watson (+ strand) and Crick (− strand) entries
watson <- transcripts %>%
  filter(Strand_Pos == "+") %>%
  group_by(BDP_ID) %>%
  slice(1) %>%
  ungroup()

crick <- transcripts %>%
  filter(Strand_Neg == "-") %>%
  group_by(BDP_ID) %>%
  slice(1) %>%
  ungroup()

# Step 5: Join Watson and Crick into BDP pairs
bdp_pairs <- inner_join(watson, crick, by = "BDP_ID", suffix = c("_watson", "_crick"))

# Step 6: Select and rename useful fields
bdp_clean <- bdp_pairs %>%
  transmute(
    BDP_ID,
    Chromosome = Chromosome_watson, #Keeping it as singular as its shared between both strands
    watson_ensembl_id = ensembl_id_watson,
    watson_gene_name = external_gene_name_watson,
    watson_start = Start_watson,
    watson_end = End_watson,
    crick_ensembl_id = ensembl_id_crick,
    crick_gene_name = external_gene_name_crick,
    crick_start = Start_crick,
    crick_end = End_crick
  )

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# Collect all unique gene IDs
all_gene_ids <- unique(c(bdp_clean$watson_ensembl_id, bdp_clean$crick_ensembl_id))
all_gene_ids <- na.omit(all_gene_ids)

# Query start and end positions
gene_coords <- getBM(
  attributes = c("ensembl_gene_id", "start_position", "end_position"),
  filters = "ensembl_gene_id",
  values = all_gene_ids,
  mart = ensembl
) %>%
  mutate(gene_length = end_position - start_position + 1) %>%
  dplyr::select(ensembl_gene_id, gene_length)

# Add lengths for both strands
bdp_clean <- bdp_clean %>%
  left_join(gene_coords, by = c("watson_ensembl_id" = "ensembl_gene_id")) %>%
  rename(watson_length = gene_length) %>%
  left_join(gene_coords, by = c("crick_ensembl_id" = "ensembl_gene_id")) %>%
  rename(crick_length = gene_length)

# Step 7: Output to CSV
write_csv(bdp_clean, "BDP_Cleaned_Pairs.csv")

#print(head(bdp_clean))
