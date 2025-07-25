if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("org.Mm.eg.db")
}

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
if (!require("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
BiocManager::install("ReactomePA")

library(ReactomePA)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)
library(readr)
library(gprofiler2)
library(tidyr)
library(ggplot2)


#--------- Foreground Setups for Both Strands and Meta Data ----------
# Extracting the individual Ensembl ids from master dataset FOR EACH STRAND and removing any NA values.
bdp <- read_csv("BDP_Fully_Annotated.csv")

watson_meta <- bdp %>%
  dplyr::select(watson_ensembl_id, watson_biotype, watson_length) %>%
  rename(ensembl_id = watson_ensembl_id,
         biotype = watson_biotype,
         length = watson_length)

crick_meta <- bdp %>%
  dplyr::select(crick_ensembl_id, crick_biotype, crick_length) %>%
  rename(ensembl_id = crick_ensembl_id,
         biotype = crick_biotype,
         length = crick_length)

foreground_metadata <- bind_rows(watson_meta, crick_meta) %>%
  filter(!is.na(ensembl_id)) %>%
  distinct() %>%
  rename(gene_length = length)

bdp_genes <- foreground_metadata$ensembl_id

#--------- Evaluation of Coverage for Enrichment ---------
unique_bdp_genes <- bdp %>%
  dplyr::select(watson_ensembl_id, crick_ensembl_id) %>%
  pivot_longer(everything(), values_to = "ensembl_id") %>%
  filter(!is.na(ensembl_id)) %>%
  distinct(ensembl_id)

total_unique_bdp_genes <- nrow(unique_bdp_genes)
inclusion_rate <- length(bdp_genes) / total_unique_bdp_genes * 100

cat("Total unique BDP-associated genes:", total_unique_bdp_genes, "\n")
cat("Genes included in enrichment:", length(bdp_genes), "\n")
cat("Inclusion rate:", round(inclusion_rate, 1), "%\n")


#--------- GO Enrichment using clusterProfiler ----------
# --------- Creating Custom Background (Universal Binning) ---------

expression_data <- read_csv("expression_dataset.csv") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Split into foreground and background gene sets
control_pool <- expression_data %>%
  filter(!ensembl_gene_id %in% bdp_genes)

foreground_metadata <- foreground_metadata %>%
  left_join(expression_data, by = c("ensembl_id" = "ensembl_gene_id")) %>%
  mutate(
    # Combine length columns (prefer expression_data's if available)
    gene_length = coalesce(gene_length.y, gene_length.x),
    # Same for biotype
    gene_biotype = coalesce(gene_biotype, biotype)
  ) %>%
  select(-gene_length.x, -gene_length.y, -biotype) %>%  # Remove duplicate columns
  filter(!is.na(gene_length) & !is.na(gene_biotype) & !is.na(baseMean))

# Add set labels and combine
combined <- bind_rows(
  foreground_metadata %>% mutate(set = "foreground"),
  control_pool %>% mutate(set = "background", ensembl_id = ensembl_gene_id)
)

# Universal binning across both sets
combined <- combined %>%
  mutate(
    expression_bin = ntile(baseMean, 4),
    length_bin = ntile(gene_length, 4),
    biotype_bin = gene_biotype
  )

# Split back into foreground and background
foreground_metadata <- combined %>%
  filter(set == "foreground") %>%
  dplyr::select(-set)

control_pool <- combined %>%
  filter(set == "background") %>%
  dplyr::select(-set)

# Count foreground genes per bin
bin_counts <- foreground_metadata %>%
  count(biotype_bin, length_bin, expression_bin, name = "n")

# Match background genes 1:1 by bin
set.seed(42)
matched_background <- bin_counts %>%
  group_by(biotype_bin, length_bin, expression_bin) %>%
  group_modify(~ {
    candidates <- control_pool %>%
      filter(
        gene_biotype == .y$biotype_bin,
        length_bin == .y$length_bin,
        expression_bin == .y$expression_bin
      ) %>%
      dplyr::select(ensembl_id, gene_biotype, gene_length, baseMean)
    
    if (nrow(candidates) >= .x$n[1]) {
      sample_n(candidates, .x$n[1])
    } else {
      candidates
    }
  }) %>%
  ungroup()

# --------- Validate bin matching ---------

foreground_bins <- foreground_metadata %>%
  count(biotype_bin, length_bin, expression_bin, name = "foreground_count")

background_bins <- matched_background %>%
  count(gene_biotype, length_bin, expression_bin, name = "background_count") %>%
  rename(biotype_bin = gene_biotype)

bin_comparison <- full_join(foreground_bins, background_bins,
                            by = c("biotype_bin", "length_bin", "expression_bin")) %>%
  mutate(match = foreground_count == background_count)

print(bin_comparison)

ggplot(bin_comparison, aes(x = factor(expression_bin), y = factor(length_bin), fill = background_count)) +
  geom_tile(color = "white") +
  facet_wrap(~ biotype_bin) +
  scale_fill_gradient(low = "lightblue", high = "steelblue") +
  labs(
    title = "Distribution of Matched Background Genes",
    x = "Expression Bin",
    y = "Length Bin",
    fill = "Background Count"
  ) +
  theme_minimal()

ggplot(bin_comparison, aes(x = factor(expression_bin), y = factor(length_bin), fill = match)) +
  geom_tile(color = "white") +
  facet_wrap(~ biotype_bin) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick")) +
  labs(
    title = "Bin Matching Status by Biotype",
    x = "Expression Bin",
    y = "Length Bin",
    fill = "Match"
  ) +
  theme_minimal

bin_comparison %>%
  group_by(biotype_bin) %>%
  summarise(total_background = sum(background_count, na.rm = TRUE)) %>%
  ggplot(aes(x = biotype_bin, y = total_background, fill = biotype_bin)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Background Gene Count by Biotype Bin", x = "Biotype Bin", y = "Total Count") +
  theme_minimal() +
  theme(legend.position = "none")
# A matched background set was generated by stratifying BDP-associated genes and non-BDP genes into bins based on gene biotype, expression quartiles, and length quartiles. One-to-one matching was performed for each bin, resulting in a background distribution identical to the foreground. No unmatched bins were detected.

write.csv(matched_background, "matched_background.csv", row.names = FALSE)

#--------- Enrichment with Custom Background ----------
#Need to map the new foreground IDs to symbols using org.Mm.eg.db to ensure they are up-to-date.

mapped_foreground <- bitr(
  foreground_metadata$ensembl_id,
  fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db
)
#length(mapped_foreground$SYMBOL)

# Map background too
mapped_background <- bitr(
  matched_background$ensembl_id,
  fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db
)
#length(mapped_background$SYMBOL)
valid_bdp_symbols <- mapped_foreground$SYMBOL

# Create combined valid universe
custom_valid_universe <- union(
  mapped_foreground$SYMBOL,
  mapped_background$SYMBOL
)
    
#Default background
GO_combined <- enrichGO(
  gene          = valid_bdp_symbols,
  universe      = keys(org.Mm.eg.db, keytype = "SYMBOL"),
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

write.csv(as.data.frame(GO_combined), "GO_BP_combined.csv", row.names = FALSE)

#Custom Background
GO_combined_custom <- enrichGO(
  gene          = valid_bdp_symbols,
  universe      = custom_valid_universe,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

write.csv(as.data.frame(GO_combined_custom), "GO_combined_custom.csv", row.names = FALSE)


GO_MF_custom <- enrichGO(
  gene          = valid_bdp_symbols,
  universe      = custom_valid_universe,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

write.csv(as.data.frame(GO_MF_custom), "GO_MF_custom.csv", row.names = FALSE)


GO_CC_custom <- enrichGO(
  gene          = valid_bdp_symbols,
  universe      = custom_valid_universe,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

write.csv(as.data.frame(GO_CC_custom), "GO_CC_custom.csv", row.names = FALSE)

#Map foreground to Entrez IDs
mapped_entrez <- bitr(
  foreground_metadata$ensembl_id,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

mapped_bg_entrez <- bitr(
  custom_valid_universe,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

bg_entrez_ids <- mapped_bg_entrez$ENTREZID
entrez_ids <- mapped_entrez$ENTREZI

kegg_enrich <- enrichKEGG(
  gene         = entrez_ids,
  universe     = bg_entrez_ids,
  organism     = "mmu",
  keyType      = "kegg",
  pvalueCutoff = 0.05
)

write.csv(as.data.frame(kegg_enrich), "kegg_enrich.csv", row.names = FALSE)


reactome_enrich <- enrichPathway(
  gene         = entrez_ids,
  universe     = bg_entrez_ids,
  organism     = "mouse",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

write.csv(as.data.frame(reactome_enrich), "reactome_enrich.csv", row.names = FALSE)

#--------- Complementary Analysis Using gProfiler to Highlight LncRNA ----------
#As 26 lncRNA genes were dropped to form the custom_valid universe, leaving us only with...

all_genes <- unique(foreground_metadata$ensembl_id)

# Non-protein-coding genes only
non_protein_coding_genes <- foreground_metadata %>%
  filter(gene_biotype != "protein_coding") %>%
  pull(ensembl_id) %>%
  unique()

cat("All BDP genes:", length(all_genes), "\n")
cat("Non-protein-coding BDP genes:", length(non_protein_coding_genes), "\n")


#------------------ Run gProfiler: All BDP Genes ------------------
gost_all <- gost(
  query = all_genes,
  organism = "mmusculus",
  sources = c("GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG", "TRANSFAC", "TF", "MIRNA"),
  correction_method = "fdr",
  evcodes = FALSE
)

flat_result <- gost_all$result %>%
  dplyr::select(-parents)  # or use purrr::map_chr if you want to keep them

# Write to CSV
write.csv(flat_result, "gprofiler_coding.csv", row.names = FALSE)

# ------------------ Run gProfiler: Non-Protein-Coding Only ------------------

gost_non_pc <- gost(
  query = non_protein_coding_genes,
  organism = "mmusculus",
  sources = c("GO:BP", "TF", "MIRNA"),
  correction_method = "fdr",
  evcodes = FALSE,
  significant=FALSE
)

flat_result <- gost_non_pc$result %>%
  dplyr::select(-parents)  # or use purrr::map_chr if you want to keep them
# Write to CSV
write.csv(flat_result, "gprofiler_non_protein_coding.csv", row.names = FALSE)
