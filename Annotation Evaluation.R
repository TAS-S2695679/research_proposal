# ---------------------------------------------------------------
# Annotation Evaluation Script (with % Coverage) for BDP Dataset
# ---------------------------------------------------------------

library(dplyr)
library(readr)

# Load the fully annotated BDP dataset
bdp <- read_csv("BDP_Fully_Annotated.csv")

# Get total number of unique Ensembl gene IDs per strand
n_watson <- bdp %>%
  filter(!is.na(watson_gene_name_updated)) %>%
  distinct(watson_ensembl_id) %>%
  nrow()

n_crick <- bdp %>%
  filter(!is.na(crick_gene_name_updated)) %>%
  distinct(crick_ensembl_id) %>%
  nrow()

# Define a helper to count annotations and calculate coverage
count_and_coverage <- function(data, prefix, total_genes) {
  gene_data <- data %>%
    dplyr::select(starts_with(prefix)) %>%
    rename_with(~ gsub(paste0(prefix, "_"), "", .x)) %>%
    mutate(gene_id = data[[paste0(prefix, "_ensembl_id")]]) %>%
    arrange(!is.na(gene_name_updated)) %>%  # Prefer rows with gene name
    distinct(gene_id, .keep_all = TRUE)
  
  tibble(
    metric = c("gene_name", "biotype", "description", "uniprot_id", "keyword", 
               "rnacentral", "gencode_biotype", "canonical_transcript"),
    count = c(
      sum(!is.na(gene_data$gene_name_updated)),
      sum(!is.na(gene_data$biotype)),
      sum(!is.na(gene_data$description)),
      sum(!is.na(gene_data$uniprot)),
      sum(!is.na(gene_data$keyword)),
      sum(!is.na(gene_data$rnacentral_id)),
      sum(!is.na(gene_data$biotype_gencode)),
      sum(!is.na(gene_data$canonical))
    )
  ) %>%
    mutate(
      coverage_percent = round(100 * count / total_genes, 1),
      strand = prefix
    )
}

# Run for Watson and Crick
watson_eval <- count_and_coverage(bdp, "watson", n_watson)
crick_eval  <- count_and_coverage(bdp, "crick", n_crick)

# Combine and arrange
annotation_eval <- bind_rows(watson_eval, crick_eval) %>%
  dplyr::select(strand, metric, count, coverage_percent)

# View results
print(annotation_eval)

bdp_pair_coverage <- bdp %>%
  mutate(watson_valid = !is.na(watson_ensembl_id),
         crick_valid = !is.na(crick_ensembl_id),
         coverage_class = case_when(
           watson_valid & crick_valid ~ "both",
           watson_valid | crick_valid ~ "one",
           TRUE ~ "none"
         )) %>%
  count(coverage_class)

print(bdp_pair_coverage)

