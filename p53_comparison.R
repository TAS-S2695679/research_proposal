# ------------------ Libraries ------------------
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(stringr)
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



# ------------------ Break Down p53 Hits in Reactome Output ------------------
# Reactome enrichment output
reactome <- read_csv("reactome_enrich.csv")

# Filter for TP53-related terms
p53_terms <- reactome %>%
  filter(str_detect(Description, "TP53"))

# Fully annotated BDP catalogue
bdp <- read_csv("BDP_Fully_Annotated.csv")

# ------------------ Prepare Long-Format BDP Metadata ------------------

# Watson strand genes
watson <- bdp %>%
  transmute(
    gene = toupper(watson_gene_name),
    ensembl_id = watson_ensembl_id,
    biotype = watson_biotype,
    strand = "Watson"
  )

# Crick strand genes
crick <- bdp %>%
  transmute(
    gene = toupper(crick_gene_name),
    ensembl_id = crick_ensembl_id,
    biotype = crick_biotype,
    strand = "Crick"
  )

# Combine both strands
bdp_long <- bind_rows(watson, crick) %>%
  filter(!is.na(gene) & gene != "")

# ------------------ Match Reactome Genes to BDPs ------------------

# Parse and annotate p53 pathway genes
p53_hits <- list()

for (i in seq_len(nrow(p53_terms))) {
  pathway_name <- p53_terms$Description[i]
  pathway_id <- p53_terms$ID[i]
  gene_list <- str_split(p53_terms$geneID[i], pattern = "/", simplify = TRUE)[1, ]
  gene_list <- toupper(gene_list)
  
  matched <- bdp_long %>%
    filter(gene %in% gene_list) %>%
    mutate(
      Reactome_Term = pathway_name,
      Reactome_ID = pathway_id
    )
  
  p53_hits[[i]] <- matched
}

# Combine all results
p53_annotated_hits <- bind_rows(p53_hits)

# ------------------ Export or View ------------------
print(p53_annotated_hits)

write_csv(p53_annotated_hits, "p53_BDP_reactome_hits.csv")


p53_hits <- read_csv("p53_BDP_reactome_hits.csv")

# ------------------ Biotype Summary ------------------
biotype_summary <- p53_hits %>%
  count(biotype, sort = TRUE)

print("Summary of Biotypes in p53-Associated BDP Genes:")
print(biotype_summary)

# ------------------ Extract Coordination Annotations ------------------
coord_data <- read_csv("BDP_coordination_results.csv")
# Watson strand mapping
watson_status <- coord_data %>%
  dplyr::select(watson_ensembl_id, brg1_coordination, oct4_coordination) %>%
  rename(
    ensembl_id = watson_ensembl_id,
    brg1_sync = brg1_coordination,
    oct4_sync = oct4_coordination
  )

# Crick strand mapping
crick_status <- coord_data %>%
  dplyr::select(crick_ensembl_id, brg1_coordination, oct4_coordination) %>%
  rename(
    ensembl_id = crick_ensembl_id,
    brg1_sync = brg1_coordination,
    oct4_sync = oct4_coordination
  )

# Combine and continue with joining to p53 hits as before
sync_status_long <- bind_rows(watson_status, crick_status) %>%
  filter(!is.na(ensembl_id)) %>%
  distinct()

# ------------------ Join with p53 Hit Table ------------------
p53_hits_annotated <- p53_hits %>%
  left_join(sync_status_long, by = "ensembl_id")

# ------------------ Summarise Results ------------------

# Brg1 summary
brg1_summary <- p53_hits_annotated %>%
  count(brg1_sync, sort = TRUE) %>%
  rename(BRG1_Coordination = brg1_sync, Count = n)

# Oct4 summary
oct4_summary <- p53_hits_annotated %>%
  count(oct4_sync, sort = TRUE) %>%
  rename(Oct4_Coordination = oct4_sync, Count = n)

# ------------------ Output ------------------
print("Coordination of p53-BDP Genes under BRG1 depletion:")
print(brg1_summary)

print("Coordination of p53-BDP Genes under OCT4 depletion:")
print(oct4_summary)

write_csv(p53_hits_annotated, "outputs/p53_hits_with_coordination.csv")
write_csv(brg1_summary, "outputs/p53_brg1_coordination_summary.csv")
write_csv(oct4_summary, "outputs/p53_oct4_coordination_summary.csv")
