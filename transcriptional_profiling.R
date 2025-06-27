library(dplyr)
library(biomaRt)
library(readr)

brg1_results <- read.delim("data/GSE87821_BRG1fl.nucRNAseq.DESeq2_Results.txt.gz", sep = ",", header = TRUE)
oct4_results <- read.delim("data/GSE87821_ZHBTC4.nucRNAseq.DESeq2_Results.txt.gz", sep = ",", header = TRUE)

# BRG1
brg1_results <- brg1_results %>%
  rename(
    refseq_id = refGeneID,
    brg1_log2fc = BRG1fl.DESeq.DOX.UNT_log2FoldChange,
    brg1_padj = BRG1fl.DESeq.PvalueAdjusted
  )

oct4_results <- oct4_results %>%
  rename(
    refseq_id = refGeneID,
    oct4_log2fc = ZHBTC4.DESeq.DOX.UNT_log2FoldChange,
    oct4_padj = ZHBTC4.DESeq.PvalueAdjusted
  )

all_refseq_ids <- unique(c(brg1_results$refseq_id, oct4_results$refseq_id))


ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

mapping <- getBM(
  attributes = c("refseq_mrna", "ensembl_gene_id"),
  filters = "refseq_mrna",
  values = all_refseq_ids,
  mart = ensembl
)

# Merge mapping
brg1_results <- brg1_results %>%
  left_join(mapping, by = c("refseq_id" = "refseq_mrna")) %>%
  filter(!is.na(ensembl_gene_id))

oct4_results <- oct4_results %>%
  left_join(mapping, by = c("refseq_id" = "refseq_mrna")) %>%
  filter(!is.na(ensembl_gene_id))


head(brg1_results)
head(oct4_results)

# Filter BRG1 results to 1 row per Ensembl ID
brg1_clean <- brg1_results %>%
  group_by(ensembl_gene_id) %>%
  arrange(brg1_padj) %>%
  slice(1) %>%
  ungroup()

# Filter oct4 results to 1 row per Ensembl ID
oct4_clean <- oct4_results %>%
  group_by(ensembl_gene_id) %>%
  arrange(oct4_padj) %>%
  slice(1) %>%
  ungroup()

#--------- Now to merge and classigy BDP by Coordination ----------
bdp <- read_csv("BDP_Fully_Annotated.csv")

bdp_expr <- bdp %>%
  left_join(
    brg1_clean %>%
      dplyr::select(ensembl_gene_id, brg1_log2fc, brg1_padj) %>%
      rename(watson_ensembl_id = ensembl_gene_id),
    by = "watson_ensembl_id"
  ) %>%
  rename(watson_log2fc = brg1_log2fc, watson_padj = brg1_padj) %>%
  
  left_join(
    brg1_clean %>%
      dplyr::select(ensembl_gene_id, brg1_log2fc, brg1_padj) %>%
      rename(crick_ensembl_id = ensembl_gene_id),
    by = "crick_ensembl_id"
  ) %>%
  rename(crick_log2fc = brg1_log2fc, crick_padj = brg1_padj) %>%
  
  # Step 2: Classify BRG1 coordination
  mutate(
    brg1_coordination = case_when(
      watson_padj < 0.05 & crick_padj < 0.05 &
        sign(watson_log2fc) == sign(crick_log2fc) ~ "Synchronised",
      watson_padj < 0.05 & crick_padj < 0.05 &
        sign(watson_log2fc) != sign(crick_log2fc) ~ "Opposite",
      xor(watson_padj < 0.05, crick_padj < 0.05) ~ "Asymmetric",
      TRUE ~ "Unchanged"
    )
  ) %>%
  
  # Step 3: Join Oct4 results for both strands
  left_join(
    oct4_clean %>%
      dplyr::select(ensembl_gene_id, oct4_log2fc, oct4_padj) %>%
      rename(watson_ensembl_id = ensembl_gene_id),
    by = "watson_ensembl_id"
  ) %>%
  rename(oct4_watson_log2fc = oct4_log2fc, oct4_watson_padj = oct4_padj) %>%
  
  left_join(
    oct4_clean %>%
      dplyr::select(ensembl_gene_id, oct4_log2fc, oct4_padj) %>%
      rename(crick_ensembl_id = ensembl_gene_id),
    by = "crick_ensembl_id"
  ) %>%
  rename(oct4_crick_log2fc = oct4_log2fc, oct4_crick_padj = oct4_padj) %>%
  
  # Step 4: Classify Oct4 coordination
  mutate(
    oct4_coordination = case_when(
      oct4_watson_padj < 0.05 & oct4_crick_padj < 0.05 &
        sign(oct4_watson_log2fc) == sign(oct4_crick_log2fc) ~ "Synchronised",
      oct4_watson_padj < 0.05 & oct4_crick_padj < 0.05 &
        sign(oct4_watson_log2fc) != sign(oct4_crick_log2fc) ~ "Opposite",
      xor(oct4_watson_padj < 0.05, oct4_crick_padj < 0.05) ~ "Asymmetric",
      TRUE ~ "Unchanged"
    )
  )


table(bdp_expr$brg1_coordination)
table(bdp_expr$oct4_coordination)


bdp_expr_export <- bdp_expr %>%
  dplyr::select(BDP_ID,
         watson_ensembl_id, crick_ensembl_id,
         brg1_coordination, oct4_coordination)

write_csv(bdp_expr_export, "BDP_coordination_results.csv")


# ------------------ Coordination Summary ------------------

# BRG1 coordination counts
brg1_summary <- bdp_expr %>%
  count(brg1_coordination, sort = TRUE) %>%
  rename(Category = brg1_coordination, Count = n)

# OCT4 coordination counts
oct4_summary <- bdp_expr %>%
  count(oct4_coordination, sort = TRUE) %>%
  rename(Category = oct4_coordination, Count = n)

# Print summaries
cat("BRG1 Coordination Summary:\n")
print(brg1_summary)

cat("\nOCT4 Coordination Summary:\n")
print(oct4_summary)

oct4_results %>% filter(ensembl_gene_id == "ENSMUSG00000059552") %>% dplyr::select(oct4_log2fc, oct4_padj)
