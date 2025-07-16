# ------------------ Libraries ------------------
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(stringr)
library(biomaRt)
library(tidyverse)
library(clusterProfiler)
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

# ------------------ Make Unique Gene List ------------------
p53_annotated_hits_unique <- p53_annotated_hits %>%
  distinct(gene, .keep_all = TRUE)

# ------------------ Export or View ------------------
print(p53_annotated_hits)
print(p53_annotated_hits_unique)
write_csv(p53_annotated_hits, "outputs/p53_BDP_reactome_hits.csv")
write_csv(p53_annotated_hits_unique, "outputs/p53_BDP_reactome_hits_unique.csv")



# ------------------ Biotype Summary ------------------
biotype_summary <- p53_annotated_hits_unique %>%
  dplyr::count(biotype, sort = TRUE)

print("Summary of Biotypes in p53-Associated BDP Genes:")
print(biotype_summary)

# ------------------ Extract Coordination Annotations ------------------
coord_data <- read_csv("BDP_coordination_results.csv")
# Watson strand mapping
watson_status <- coord_data %>%
  dplyr::select(watson_ensembl_id, brg1_coordination, oct4_coordination) %>%
  dplyr::rename(
    ensembl_id = watson_ensembl_id,
    brg1_sync = brg1_coordination,
    oct4_sync = oct4_coordination
  )

# Crick strand mapping
crick_status <- coord_data %>%
  dplyr::select(crick_ensembl_id, brg1_coordination, oct4_coordination) %>%
  dplyr::rename(
    ensembl_id = crick_ensembl_id,
    brg1_sync = brg1_coordination,
    oct4_sync = oct4_coordination
  )

# Combine and continue with joining to p53 hits as before
sync_status_long <- bind_rows(watson_status, crick_status) %>%
  filter(!is.na(ensembl_id)) %>%
  distinct()

# ------------------ Join with p53 Hit Table ------------------
p53_hits_annotated <- p53_annotated_hits_unique %>%
  left_join(sync_status_long, by = "ensembl_id")

# ------------------ Summarise Results ------------------

# Brg1 summary
brg1_summary <- p53_hits_annotated %>%
  dplyr::count(brg1_sync, sort = TRUE) %>%
  dplyr::rename(BRG1_Coordination = brg1_sync, Count = n)

# Oct4 summary
oct4_summary <- p53_hits_annotated %>%
  dplyr::count(oct4_sync, sort = TRUE) %>%
  dplyr::rename(Oct4_Coordination = oct4_sync, Count = n)

# ------------------ Output ------------------
print("Coordination of p53-BDP Genes under BRG1 depletion:")
print(brg1_summary)

print("Coordination of p53-BDP Genes under OCT4 depletion:")
print(oct4_summary)

write_csv(p53_hits_annotated, "outputs/p53_hits_with_coordination.csv")
write_csv(brg1_summary, "outputs/p53_brg1_coordination_summary.csv")
write_csv(oct4_summary, "outputs/p53_oct4_coordination_summary.csv")



asymmetric_hits <- p53_hits_annotated %>%
  filter(oct4_sync == "Asymmetric")
print(asymmetric_hits)

bdp <- read_csv("BDP_Fully_Annotated.csv")
asymmetric_matches <- asymmetric_hits %>%
  left_join(bdp, by = c("ensembl_id" = "watson_ensembl_id")) %>%
  mutate(partner_ensembl = crick_ensembl_id,
         partner_symbol = crick_gene_name) %>%
  bind_rows(
    asymmetric_hits %>%
      left_join(bdp, by = c("ensembl_id" = "crick_ensembl_id")) %>%
      mutate(partner_ensembl = watson_ensembl_id,
             partner_symbol = watson_gene_name)
  ) %>%
  filter(!is.na(partner_ensembl)) %>%
  dplyr::select(gene, ensembl_id, strand, partner_symbol, partner_ensembl)

print(asymmetric_matches)


genes_to_check <- c("ENSMUSG00000025358",  # Cdk2
                    "ENSMUSG00000025359",  # Pmel
                    "ENSMUSG00000064128",  # Cenpj
                    "ENSMUSG00000054507")  # Parp4

oct4_clean %>%
  filter(ensembl_gene_id %in% genes_to_check) %>%
  dplyr::select(ensembl_gene_id, oct4_log2fc, oct4_padj)

# ------------------ Sanity Check: Verify Classification ------------------

# Helper: assess if gene is significantly changed
is_significant <- function(log2fc, padj, threshold = 0.58, pval_cutoff = 0.05) {
  !is.na(log2fc) && !is.na(padj) && abs(log2fc) > threshold && padj < pval_cutoff
}

# Add flags for significance
bdp_checked <- bdp_expr %>%
  mutate(
    watson_sig = mapply(is_significant, oct4_watson_log2fc, oct4_watson_padj),
    crick_sig  = mapply(is_significant, oct4_crick_log2fc, oct4_crick_padj)
  )

# Filter cases where classification ≠ what the data actually shows
bdp_mismatch <- bdp_checked %>%
  mutate(
    expected = case_when(
      watson_sig & crick_sig &
        sign(oct4_watson_log2fc) == sign(oct4_crick_log2fc) ~ "Synchronised",
      watson_sig & crick_sig &
        sign(oct4_watson_log2fc) != sign(oct4_crick_log2fc) ~ "Opposite",
      xor(watson_sig, crick_sig) ~ "Asymmetric",
      TRUE ~ "Unchanged"
    )
  ) %>%
  filter(oct4_coordination != expected) %>%
  dplyr::select(BDP_ID,
                watson_ensembl_id, crick_ensembl_id,
                oct4_watson_log2fc, oct4_crick_log2fc,
                oct4_watson_padj, oct4_crick_padj,
                oct4_coordination, expected)



# CDK2 and PMEL
oct4_clean %>%
  filter(ensembl_gene_id %in% c("ENSMUSG00000025358", "ENSMUSG00000025359")) %>%
  select(ensembl_gene_id, oct4_log2fc, oct4_padj)

# CENPJ and PARP4
oct4_clean %>%
  filter(ensembl_gene_id %in% c("ENSMUSG00000064128", "ENSMUSG00000054509")) %>%
  select(ensembl_gene_id, oct4_log2fc, oct4_padj)


# ------------------ Cross-Referencing GO BP Terms with BDP and p53 Hits for validation ------------------

go_bp <- read_csv("GO_BP_combined.csv")  # BP enrichment results (default background)
bdp <- read_csv("BDP_Fully_Annotated.csv")   # Annotated BDP catalogue
p53_hits <- read_csv("p53_BDP_reactome_hits.csv")  # Genes from TP53-related Reactome pathways

# ------------------ 1. Filter for p53-related GO terms ------------------
keywords <- c("apoptosis", "cell cycle arrest", "DNA damage", "TP53", "p53", "stress", "senescence", "death")

filtered_go <- go_bp %>%
  filter(str_detect(tolower(Description), str_c(keywords, collapse = "|")))

# ------------------ 2. Parse and flatten gene lists ------------------
go_gene_list <- filtered_go %>%
  separate_rows(geneID, sep = "/") %>%
  mutate(gene = toupper(geneID)) %>%
  distinct(gene)

# ------------------ 3. Prepare long-format BDP gene list ------------------
bdp_long <- bind_rows(
  bdp %>%
    transmute(
      gene = toupper(watson_gene_name),
      ensembl_id = watson_ensembl_id,
      strand = "Watson"
    ),
  bdp %>%
    transmute(
      gene = toupper(crick_gene_name),
      ensembl_id = crick_ensembl_id,
      strand = "Crick"
    )
) %>%
  filter(!is.na(gene) & gene != "")

# ------------------ 4. Match GO genes to BDPs ------------------
go_bdp_hits <- go_gene_list %>%
  inner_join(bdp_long, by = "gene")  # brings in Ensembl ID and strand

# ------------------ 5. Compare with TP53 Reactome hits ------------------
p53_genes <- unique(toupper(p53_hits$gene))

go_bdp_hits <- go_bdp_hits %>%
  mutate(is_in_p53_reactome = gene %in% p53_genes)

# ------------------ 6. Summary statistics ------------------
summary_counts <- go_bdp_hits %>%
  summarise(
    total_GO_genes_in_BDPs = n(),
    overlap_with_p53_reactome = sum(is_in_p53_reactome)
  )

# ------------------ 7. Save results ------------------
write_csv(go_bdp_hits, "outputs/GO_BP_BDP_TP53_overlap_genes.csv")
write_csv(summary_counts, "outputs/GO_BP_BDP_TP53_overlap_summary.csv")

# Print summary
print(summary_counts)



# ------------------ ENRICHMENT OF P53 REACTOME RESULTS ------------------

p53_bdp <- read_csv("outputs/p53_BDP_reactome_hits_unique.csv")

# Assume column is called "ensembl_id" (change if needed)
p53_genes <- unique(p53_bdp$ensembl_id)

# Load your full BDP background (Watson + Crick Ensembl IDs)
bdp <- read_csv("BDP_Fully_Annotated.csv")

bdp_genes <- unique(c(bdp$watson_ensembl_id, bdp$crick_ensembl_id)) %>%
  na.omit()

# Run GO BP enrichment using BDPs as background
p53_enrich <- enrichGO(
  gene          = p53_genes,
  universe      = bdp_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  readable      = TRUE
)

# Save results
write.csv(as.data.frame(p53_enrich), "outputs/p53_bdp_go_enrichment.csv", row.names = FALSE)

# Quick view
head(as.data.frame(p53_enrich))
