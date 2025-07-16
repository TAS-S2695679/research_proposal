if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

library(dplyr)
library(readr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tibble)
library(biomaRt)
library(tidyr)


# Load data
bdp <- read_csv("BDP_coordination_results.csv")
bdp_genes <- unique(c(bdp$watson_ensembl_id, bdp$crick_ensembl_id))

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

custom_universe <- unique(c(brg1_clean$ensembl_gene_id, oct4_clean$ensembl_gene_id))

all_bdp <- read_csv("BDP_Fully_Annotated.csv")
all_bdp_genes <- unique(c(all_bdp$watson_ensembl_id, all_bdp$crick_ensembl_id))
all_bdp_genes <- all_bdp_genes[!is.na(all_bdp_genes)]


asymmetric_genes_symbol <- bitr(
  asymmetric_genes,
  fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db
)

BDP_genes_symbol <- bitr(
  all_bdp_genes,
  fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db
)
as# Run GO enrichment: ASYMMETRIC
asym_enrich <- enrichGO(
  gene          = asymmetric_genes,
  universe      = all_bdp_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "none",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  readable      = TRUE
)

# Run GO enrichment: SYNCHRONISED
sync_enrich <- enrichGO(
  gene          = synchronised_genes,
  universe      = all_bdp_genes,
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

mapped_entrez <- bitr(
  asymmetric_genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

mapped_bg_entrez <- bitr(
  all_bdp_genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

bg_entrez_ids <- mapped_bg_entrez$ENTREZID
entrez_ids <- mapped_entrez$ENTREZID

asym_reactome_enrich <- enrichPathway(
  gene         = entrez_ids,
  universe     = bg_entrez_ids,
  organism     = "mouse",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

write.csv(as.data.frame(asym_reactome_enrich), "outputs/unchanged_bdp_go_enrichment_custom_reactome.csv", row.names = FALSE)



#-----------Descriptive Functional Labelling-----------
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")


mapping <- getBM(
  attributes = c("ensembl_gene_id",
                 "external_gene_name",
                 "description",
                 "gene_biotype"
  ),
  filters = "ensembl_gene_id",
  values = asymmetric_genes,
  mart = ensembl
)
mapping



# Step 1: Load your CSV
df <- read_csv("BDP_coordination_results.csv")

# Step 2: Connect to Ensembl BioMart
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# Step 3: Get unique gene IDs
all_genes <- unique(c(df$watson_ensembl_id, df$crick_ensembl_id))

# Step 4: Fetch gene biotypes
biotypes <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = all_genes,
  mart = ensembl
)

# Step 5: Join biotypes back to main dataframe
df <- df %>%
  left_join(biotypes, by = c("watson_ensembl_id" = "ensembl_gene_id")) %>%
  rename(watson_biotype = gene_biotype) %>%
  left_join(biotypes, by = c("crick_ensembl_id" = "ensembl_gene_id")) %>%
  rename(crick_biotype = gene_biotype)

# Step 6: Classify pair types
df <- df %>%
  mutate(pair_type = case_when(
    watson_biotype == "protein_coding" & crick_biotype == "protein_coding" ~ "PCG–PCG",
    watson_biotype == "protein_coding" | crick_biotype == "protein_coding" ~ "PCG–lncRNA",
    TRUE ~ "lncRNA–lncRNA"
  ))

# Step 7: Count combinations per coordination class
summary_table <- df %>%
  count(oct4_coordination, pair_type) %>%
  pivot_wider(names_from = pair_type, values_from = n, values_fill = 0)

# Step 8: View the result
print(summary_table)
