library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)

genes_by_class <- bdp_expr %>%
  dplyr::select(oct4_coordination, watson_ensembl_id, crick_ensembl_id) %>%
  pivot_longer(cols = c(watson_ensembl_id, crick_ensembl_id), values_to = "ensembl_id") %>%
  filter(!is.na(ensembl_id)) %>%
  distinct() %>%
  group_by(oct4_coordination) %>%
  summarise(gene_ids = list(unique(ensembl_id)))

run_enrichment <- function(genes) {
  enrichGO(
    gene         = genes,
    OrgDb        = org.Mm.eg.db,
    keyType      = "ENSEMBL",
    ont          = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
}

enrichment_results <- lapply(setNames(genes_by_class$gene_ids, genes_by_class$oct4_coordination), run_enrichment)
dotplot(enrichment_results[["Synchronised"]], showCategory = 20, title = "GO BP Enrichment in Oct4-Synchronised BDPs")
s
write.csv(as.data.frame(enrichment_results[["Synchronised"]]), "oct4_synchronised_GO.csv")

oct4_sync_genes <- bdp_expr %>%
  filter(oct4_coordination == "Synchronised") %>%
  dplyr::select(watson_ensembl_id, crick_ensembl_id) %>%
  pivot_longer(everything(), values_to = "ensembl_id") %>%
  distinct(ensembl_id) %>%
  filter(!is.na(ensembl_id)) %>%
  pull(ensembl_id)

entrez_ids <- bitr(oct4_sync_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
kegg_results <- enrichKEGG(
  gene = entrez_ids$ENTREZID,
  organism = 'mmu',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

head(kegg_results)
