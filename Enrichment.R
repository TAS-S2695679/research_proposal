if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("org.Mm.eg.db")
}

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}


library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)
library(readr)
library(gprofiler2)

#--------- Foreground Setups for Both Strands and Meta Data ----------
bdp <- read_csv("BDP_Fully_Annotated.csv")

# Extracting the individual Ensembl ids from master dataset FOR EACH STRAND and removing any NA values.
watson_genes <- bdp %>%
  dplyr::select(watson_ensembl_id) %>%
  distinct() %>%
  filter(!is.na(watson_ensembl_id)) %>%
  pull(watson_ensembl_id)

crick_genes <- bdp %>%
  dplyr::select(crick_ensembl_id) %>%
  distinct() %>%
  filter(!is.na(crick_ensembl_id)) %>%
  pull(crick_ensembl_id)

watson_meta <- bdp %>%
  dplyr::select(watson_ensembl_id, watson_biotype, watson_length) %>%
  distinct() %>%
  dplyr::rename(ensembl_id = watson_ensembl_id,
         biotype = watson_biotype,
         length = watson_length)

crick_meta <- bdp %>%
  dplyr::select(crick_ensembl_id, crick_biotype, crick_length) %>%
  distinct() %>%
  dplyr::rename(ensembl_id = crick_ensembl_id,
         biotype = crick_biotype,
         length = crick_length)

# Combine both
foreground_metadata <- bind_rows(watson_meta, crick_meta) %>%
  filter(!is.na(ensembl_id)) %>%
  distinct()



#--------- Default Background Setup ----------
background_genes <- AnnotationDbi::keys(org.Mm.eg.db, keytype = "ENSEMBL")

bdp_genes <- unique(c(watson_genes, crick_genes))

background_genes_filtered <- setdiff(background_genes, bdp_genes)

#Map ensembl ids to symbols
mapped_bdp <- bitr(bdp_genes, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
valid_bdp_symbols <- mapped_bdp$SYMBOL
#head(watson_genes)
#head(background_genes_filtered)

default_universe_symbols <- keys(org.Mm.eg.db, keytype = "SYMBOL")

#--------- Evaluation of Coverage for Enrichment ---------
total_bdp_genes <- nrow(bdp) * 2
inclusion_rate <- length(bdp_genes) / total_bdp_genes * 100
cat("Total BDP-associated genes:", total_bdp_genes, "\n")
cat("Genes included in enrichment:", length(bdp_genes), "\n")
cat("Inclusion rate:", round(inclusion_rate, 1), "%\n")

#--------- GO Enrichment using clusterProfiler ----------

GO_combined <- enrichGO(
  gene          = valid_bdp_symbols,
  universe      = default_universe_symbols,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

write.csv(as.data.frame(GO_combined), "GO_BP_combined.csv", row.names = FALSE)

GO_watson <- enrichGO(
  gene          = watson_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",                  # You can repeat for MF and CC separately
  pAdjustMethod = "BH",                 # Benjamini-Hochberg FDR control
  pvalueCutoff  = 0.05,                 # Strict threshold
  qvalueCutoff  = 0.05,                 # Adds redundancy control
  minGSSize     = 10,                   # Ignore overly small GO categories
  maxGSSize     = 500,                  # Exclude overly general ones
  readable      = TRUE                  # Convert Ensembl to gene symbols
)

GO_watson

#head(GO_watson)

as.data.frame(GO_watson)[100:200, c("ID", "Description", "GeneRatio", "p.adjust")]

# Dotplot
dotplot(GO_watson, showCategory = 20, title = "GO:BP Enrichment – Watson Strand")

write.csv(as.data.frame(GO_watson), "GO_BP_Watson.csv", row.names = FALSE)


GO_crick <- enrichGO(
  gene          = crick_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",                  # You can repeat for MF and CC separately
  pAdjustMethod = "BH",                 # Benjamini-Hochberg FDR control
  pvalueCutoff  = 0.05,                 # Strict threshold
  qvalueCutoff  = 0.05,                 # Adds redundancy control
  minGSSize     = 10,                   # Ignore overly small GO categories
  maxGSSize     = 500,                  # Exclude overly general ones
  readable      = TRUE                  # Convert Ensembl to gene symbols
)

GO_crick

#head(GO_crick)

as.data.frame(GO_crick)[100:200, c("ID", "Description", "GeneRatio", "p.adjust")]

# Dotplot
dotplot(GO_crick, showCategory = 20, title = "GO:BP Enrichment – Crick Strand")


write.csv(as.data.frame(GO_crick), "GO_BP_Crick.csv", row.names = FALSE)


#--------- GO:MF Enrichment ----------

GO_watson_MF <- enrichGO(
  gene          = watson_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = TRUE
)

write.csv(as.data.frame(GO_watson_MF), "GO_MF_Watson.csv", row.names = FALSE)
dotplot(GO_watson_MF, showCategory = 20, title = "GO:MF Enrichment – Watson Strand")


GO_crick_MF <- enrichGO(
  gene          = crick_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = TRUE
)

write.csv(as.data.frame(GO_crick_MF), "GO_MF_Crick.csv", row.names = FALSE)
dotplot(GO_crick_MF, showCategory = 20, title = "GO:MF Enrichment – Crick Strand")

#--------- GO:CC Enrichment ----------

GO_watson_CC <- enrichGO(
  gene          = watson_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = TRUE
)

write.csv(as.data.frame(GO_watson_CC), "GO_CC_Watson.csv", row.names = FALSE)
dotplot(GO_watson_CC, showCategory = 20, title = "GO:CC Enrichment – Watson Strand")


GO_crick_CC <- enrichGO(
  gene          = crick_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = TRUE
)

write.csv(as.data.frame(GO_crick_CC), "GO_CC_Crick.csv", row.names = FALSE)
dotplot(GO_crick_CC, showCategory = 20, title = "GO:CC Enrichment – Crick Strand")


#--------- Creating Custom Background ----------
expression_data <- read_csv("expression_dataset.csv")

expression_data <- expression_data %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

control_pool <- expression_data %>%
  filter(!ensembl_gene_id %in% bdp_genes)

#count(control_pool)
#length(bdp_genes)
#count(expression_data)

foreground_metadata <- foreground_metadata %>%
  left_join(expression_data, by = c("ensembl_id" = "ensembl_gene_id")) %>%
  filter(!is.na(gene_length) & !is.na(gene_biotype) & !is.na(baseMean)) %>%
  mutate(
    expression_bin = ntile(baseMean, 4),
    length_bin = ntile(gene_length, 4),
    biotype_bin = gene_biotype
  )

control_pool <- control_pool %>%
  mutate(
    expression_bin = ifelse(is.na(baseMean), "missing", ntile(baseMean, 4)),
    length_bin = ntile(gene_length, 4),
    biotype_bin = gene_biotype
  )

bin_counts <- foreground_metadata %>%
  group_by(biotype_bin, length_bin, expression_bin) %>%
  summarise(n = n(), .groups = "drop")

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
      dplyr::select(-biotype_bin, -length_bin, -expression_bin)  #strip grouping cols
    
    if (nrow(candidates) >= .x$n[1]) {
      sample_n(candidates, .x$n[1])
    } else {
      candidates
    }
  }) %>%
  ungroup()

#Quick check to see if there are any missmatches in the bins

foreground_bins <- foreground_metadata %>%
  count(biotype_bin, length_bin, expression_bin, name = "foreground_count")

# Count genes per bin in the matched background
background_bins <- matched_background %>%
  count(gene_biotype, length_bin, expression_bin, name = "background_count") %>%
  rename(biotype_bin = gene_biotype)

# Join and compare the two sets
bin_comparison <- full_join(foreground_bins, background_bins,
                            by = c("biotype_bin", "length_bin", "expression_bin")) %>%
  mutate(match = foreground_count == background_count)

# View the results
print(bin_comparison)



# A matched background set was generated by stratifying BDP-associated genes and non-BDP genes into bins based on gene biotype, expression quartiles, and length quartiles. One-to-one matching was performed for each bin, resulting in a background distribution identical to the foreground. No unmatched bins were detected.

write.csv(matched_background, "matched_background.csv", row.names = FALSE)

#Of the 5,356 genes associated with bidirectional promoters (BDPs), 5,335 (99.6%) were successfully mapped from RefSeq transcript IDs to Ensembl gene identifiers using biomaRt. A small number (21 genes) could not be mapped and were excluded from downstream analysis.
#Following mapping, biotype annotations were retrieved from Ensembl. A total of 1,438 BDP-mapped genes lacked a defined biotype and were filtered out to ensure compatibility with background matching. This resulted in 3,897 BDP genes retained for custom enrichment analysis using a matched control background. These represent the subset with both valid expression measurements and functional annotations.
#The remaining 1,459 BDP genes (21 unmapped + 1,438 unannotated) were included in the initial default-background enrichment analysis (using org.Mm.eg.db) but excluded from custom background enrichment to avoid bias introduced by unmatched gene properties.


#--------- Enrichment with Custom Background ----------
#Need to map the new foreground IDs to symbols
mapped_bdp <- bitr(foreground_metadata$ensembl_id, 
                   fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)

valid_bdp_symbols <- mapped_bdp$SYMBOL

custom_valid_universe <- bitr(
  union(foreground_metadata$ensembl_id, matched_background$ensembl_gene_id),
  fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db
)$SYMBOL

GO_combined_custom <- enrichGO(
  gene          = valid_bdp_symbols,
  universe      = custom_valid_universe,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

write.csv(as.data.frame(GO_combined_custom), "GO_combined_custom.csv", row.names = FALSE)


#--------- Complementary Analysis Using gProfiler to Highlight LncRNA ----------
#As 26 lncRNA genes were dropped to form the custom_valid universe, leaving us only with 

