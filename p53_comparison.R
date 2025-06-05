# ------------------ Libraries ------------------
library(dplyr)
library(readr)
library(tibble)

# ------------------ Load Data ------------------

# Full set of expressed genes
expression_data <- read_csv("expression_dataset.csv")
ensembl_ids <- unique(expression_data$ensembl_gene_id)

# Set up mouse BioMart
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# Retrieve gene symbols for those Ensembl IDs
mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                 filters = "ensembl_gene_id",
                 values = ensembl_ids,
                 mart = ensembl)

# Clean up: remove empty mappings and duplicates
expressed_genes <- mapping$external_gene_name %>%
  na.omit() %>%
  unique() %>%
  toupper()

# Your BDP-associated genes
bdp_genes <- read_csv("BDP_Fully_Annotated.csv") %>%
  dplyr::select(watson_gene_name, crick_gene_name) %>%
  pivot_longer(everything(), values_to = "gene") %>%
  pull(gene) %>%
  na.omit() %>%
  toupper() %>%
  unique()

# Load full p53 target dataset from Table S2 (from PMID:22387025)
p53_table <- read_csv("data/NIHMS353782-supplement-03.csv")

# ------------------ Split p53 Target List ------------------

# Column names may vary slightly — adjust if needed
p53_targets <- p53_table %>%
  dplyr::select(Gene = `Gene Symbol`, Ratio = 4) %>%  # column 4 is fold change
  filter(!is.na(Gene)) %>%
  mutate(
    Gene = toupper(Gene),
    Regulation = case_when(
      Ratio > 1 ~ "Activated",
      Ratio < 1 ~ "Repressed",
      TRUE ~ NA_character_
    )
  )

# Split into gene symbol vectors
p53_activated <- p53_targets %>%
  filter(Regulation == "Activated") %>%
  pull(Gene) %>%
  unique()

p53_repressed <- p53_targets %>%
  filter(Regulation == "Repressed") %>%
  pull(Gene) %>%
  unique()

# ------------------ Run Fisher Tests ------------------

run_fisher_test <- function(bdp, target, universe) {
  a <- sum(bdp %in% target)
  b <- sum(!(bdp %in% target))
  c <- sum((universe %in% target) & !(universe %in% bdp))
  d <- sum(!(universe %in% union(bdp, target)))
  
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                dimnames = list(BDP = c("p53 target", "non-target"),
                                Background = c("BDP", "non-BDP")))
  
  result <- fisher.test(mat)
  list(
    p_value = result$p.value,
    odds_ratio = result$estimate,
    count_BDP_overlap = a,
    total_BDP = length(bdp),
    total_target = length(target),
    matrix = mat
  )
}

# Run enrichment tests
activated_result <- run_fisher_test(bdp_genes, p53_activated, expressed_genes)
repressed_result <- run_fisher_test(bdp_genes, p53_repressed, expressed_genes)

# ------------------ Summarise Results ------------------

summary_df <- tibble(
  Set = c("p53 Activated", "p53 Repressed"),
  Overlap_with_BDPs = c(activated_result$count_BDP_overlap, repressed_result$count_BDP_overlap),
  Total_Expressed_BDP_genes = length(bdp_genes), # Only BDP genes found in the expressed gene universe (after BioMart mapping)
  Total_p53_targets = c(length(p53_activated), length(p53_repressed)),
  Odds_Ratio = c(activated_result$odds_ratio, repressed_result$odds_ratio),
  P_Value = c(activated_result$p_value, repressed_result$p_value)
)

print(summary_df)

write_csv(summary_df, "outputs/p53_comparison/summary.csv")
