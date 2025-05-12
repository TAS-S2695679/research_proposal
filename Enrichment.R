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

#--------- Foreground Setups for Both Strands and Meta Data ----------
bdp <- read_csv("BDP_Fully_Annotated.csv")

# Extracting the individual Ensembl ids from master dataset FOR EACH STRAND and removing any NA values.
watson_genes <- bdp %>%
  select(watson_ensembl_id) %>%
  distinct() %>%
  filter(!is.na(watson_ensembl_id)) %>%
  pull(watson_ensembl_id)

crick_genes <- bdp %>%
  select(crick_ensembl_id) %>%
  distinct() %>%
  filter(!is.na(crick_ensembl_id)) %>%
  pull(crick_ensembl_id)

watson_meta <- bdp %>%
  select(watson_ensembl_id, watson_biotype, watson_length) %>%
  distinct() %>%
  rename(ensembl_id = watson_ensembl_id,
         biotype = watson_biotype,
         length = watson_length)

crick_meta <- bdp %>%
  select(crick_ensembl_id, crick_biotype, crick_length) %>%
  distinct() %>%
  rename(ensembl_id = crick_ensembl_id,
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


#head(watson_genes)
#head(background_genes_filtered)

#--------- GO Enrichment using clusterProfiler ----------

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