if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

library(dplyr)
library(readr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tibble)

# Load data
bdp <- read_csv("BDP_coordination_results.csv")

# Flatten into a long format: one row per gene
bdp_long <- bdp %>%
  transmute(
    ensembl_id = watson_ensembl_id,
    coordination = oct4_coordination
  ) %>%
  bind_rows(
    bdp %>%
      transmute(
        ensembl_id = crick_ensembl_id,
        coordination = oct4_coordination
      )
  ) %>%
  filter(!is.na(ensembl_id)) %>%
  distinct()

# Prepare gene sets
asymmetric_genes <- bdp_long %>%
  filter(coordination == "Asymmetric") %>%
  pull(ensembl_id)

synchronised_genes <- bdp_long %>%
  filter(coordination == "Synchronised") %>%
  pull(ensembl_id)

bdp_background_genes <- bdp_long %>%
  pull(ensembl_id) %>%
  unique()

unchanged_genes <- bdp_long %>%
  filter(coordination == "Unchanged") %>%
  pull(ensembl_id)

# Run GO enrichment: ASYMMETRIC
asym_enrich <- enrichGO(
  gene          = asymmetric_genes,
  universe      = bdp_background_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

# Run GO enrichment: SYNCHRONISED
sync_enrich <- enrichGO(
  gene          = synchronised_genes,
  universe      = bdp_background_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

# Save results
write.csv(as.data.frame(asym_enrich), "outputs/asymmetric_bdp_go_enrichment.csv", row.names = FALSE)
write.csv(as.data.frame(sync_enrich), "outputs/synchronised_bdp_go_enrichment.csv", row.names = FALSE)

# Quick check
print(head(asym_enrich))
print(head(sync_enrich))

#From custom background creation in original enrichment
custom_valid_universe <- union(
  foreground_metadata$ensembl_id, matched_background$ensembl_id)

asym_enrich <- enrichGO(
  gene          = asymmetric_genes,
  universe      = custom_valid_universe,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

# Run GO enrichment: SYNCHRONISED
sync_enrich <- enrichGO(
  gene          = synchronised_genes,
  universe      = custom_valid_universe,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)


unchanged_enrich <- enrichGO(
  gene          = unchanged_genes,
  universe      = custom_valid_universe,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)


write.csv(as.data.frame(asym_enrich), "outputs/asymmetric_bdp_go_enrichment_custom.csv", row.names = FALSE)
write.csv(as.data.frame(sync_enrich), "outputs/synchronised_bdp_go_enrichment_custom.csv", row.names = FALSE)
write.csv(as.data.frame(unchanged_enrich), "outputs/unchanged_bdp_go_enrichment_custom.csv", row.names = FALSE)
