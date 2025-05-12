library(dplyr)
library(biomaRt)

expr_raw<- read.delim("GSE87821_nucRNAseq_ReadCount_Quantitation.txt.gz")

expr_raw$baseMean <- rowMeans(expr_raw[, grep("ReadCount", colnames(expr_raw))])

head(expr_raw$baseMean)

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

mapping <- getBM(
  attributes = c("refseq_mrna", "ensembl_gene_id"),
  filters = "refseq_mrna",
  values = expr_raw$ID,
  mart = ensembl
)

expr_mapped <- expr_raw %>%
  left_join(mapping, by = c("ID" = "refseq_mrna")) %>%
  filter(!is.na(ensembl_gene_id))

head(expr_mapped)

expression_data <- expr_mapped %>%
  select(ensembl_gene_id, baseMean, gene_length = Size) %>%
  distinct()

biotypes <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = expression_data$ensembl_gene_id,
  mart = ensembl
)

expression_data <- expression_data %>%
  left_join(biotypes, by = "ensembl_gene_id") %>%
  filter(!is.na(gene_biotype))


write.csv(as.data.frame(expression_data), "expression_dataset.csv", row.names = FALSE)