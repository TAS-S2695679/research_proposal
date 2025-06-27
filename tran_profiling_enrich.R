# ---- Load packages ----
library(dplyr)
library(readr)
library(tidyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

# ---- Load data ----
bdp <- read_csv("BDP_Fully_Annotated.csv")
coord <- read_csv("BDP_coordination_results.csv")

# ---- Join and reshape gene lists per category ----
bdp_long <- bdp %>%
  dplyr::select(BDP_ID, watson_ensembl_id, crick_ensembl_id)

# Merge coordination class from Oct4 results
bdp_classified <- coord %>%
  pivot_longer(cols = c(watson_ensembl_id, crick_ensembl_id),
               names_to = "strand", values_to = "ensembl_id") %>%
  filter(!is.na(ensembl_id))


# ---- Define categories ----
categories <- unique(bdp_classified$oct4_coordination)

# ---- Map Ensembl IDs to SYMBOLs ----
all_ids <- unique(bdp_classified$ensembl_id)
mapped_genes <- bitr(all_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)

# ---- Run enrichment per category ----
for (cat in categories) {
  cat("\n🔍 Enriching category:", cat, "\n")
  
  # Get genes in this category
  gene_set <- bdp_classified %>%
    filter(oct4_coordination == cat) %>%
    pull(ensembl_id) %>%
    unique()
  
  # Map to gene symbols
  gene_symbols <- mapped_genes %>%
    filter(ENSEMBL %in% gene_set) %>%
    pull(SYMBOL)
  
  if (length(gene_symbols) >= 10) {
    enrich_res <- enrichGO(
      gene          = gene_symbols,
      OrgDb         = org.Mm.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      readable      = TRUE
    )
    
    out_file <- paste0("GO_BP_enrichment_", cat, ".csv")
    write.csv(as.data.frame(enrich_res), out_file, row.names = FALSE)
    cat("✅ Saved to:", out_file, "\n")
  } else {
    cat("⚠️  Skipping", cat, "- not enough genes for enrichment.\n")
  }
}
