library(dplyr)
library(readr)

# Load both datasets
cleaned <- read_csv("BDP_Cleaned_Pairs.csv")
annotated <- read_csv("BDP_Annotated_Pairs.csv")

# --- Compare Watson gene names ---
watson_filled <- annotated %>%
  filter(is.na(watson_gene_name) & !is.na(watson_gene_name_updated)) %>%
  select(BDP_ID, watson_ensembl_id, watson_gene_name_updated, watson_biotype, watson_description)

# --- Compare Crick gene names ---
crick_filled <- annotated %>%
  filter(is.na(crick_gene_name) & !is.na(crick_gene_name_updated)) %>%
  select(BDP_ID, crick_ensembl_id, crick_gene_name_updated, crick_biotype, crick_description)

print("Watson genes updated:")
print(watson_filled)

print("Crick genes updated:")
print(crick_filled)
