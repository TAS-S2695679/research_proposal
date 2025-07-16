library(dplyr)


valid_oct4_ids <- unique(oct4_results$ensembl_gene_id)

watson_check <- bdp_expr %>%
  dplyr::select(BDP_ID, watson_ensembl_id, oct4_watson_log2fc, oct4_watson_padj) %>%
  mutate(
    watson_na_reason = case_when(
      is.na(watson_ensembl_id) ~ "No Ensembl ID (unmapped)",
      !watson_ensembl_id %in% valid_oct4_ids ~ "Not in DE results (annotation issue)",
      is.na(oct4_watson_log2fc) | is.na(oct4_watson_padj) ~ "Low expression (not tested)",
      TRUE ~ "OK"
    )
  )

crick_check <- bdp_expr %>%
  dplyr::select(BDP_ID, crick_ensembl_id, oct4_crick_log2fc, oct4_crick_padj) %>%
  mutate(
    crick_na_reason = case_when(
      is.na(crick_ensembl_id) ~ "No Ensembl ID (unmapped)",
      !crick_ensembl_id %in% valid_oct4_ids ~ "Not in DE results (annotation issue)",
      is.na(oct4_crick_log2fc) | is.na(oct4_crick_padj) ~ "Low expression (not tested)",
      TRUE ~ "OK"
    )
  )

coordination_diagnostics <- bdp_expr %>%
  left_join(watson_check %>% dplyr::select(BDP_ID, watson_na_reason), by = "BDP_ID") %>%
  left_join(crick_check %>% dplyr::select(BDP_ID, crick_na_reason), by = "BDP_ID")
watson_summary <- table(coordination_diagnostics$watson_na_reason)
crick_summary <- table(coordination_diagnostics$crick_na_reason)


cat("Watson NA reasons:\n")
print(watson_summary)

cat("\nCrick NA reasons:\n")
print(crick_summary)



total_oct4 <- nrow(oct4_results)

# Reload full Oct4 before filtering
oct4_raw <- read.delim("data/GSE87821_ZHBTC4.nucRNAseq.DESeq2_Results.txt.gz", sep = ",", header = TRUE)

# Rename to match previous
oct4_raw <- oct4_raw %>%
  rename(
    refseq_id = refGeneID,
    oct4_log2fc = ZHBTC4.DESeq.DOX.UNT_log2FoldChange,
    oct4_padj = ZHBTC4.DESeq.PvalueAdjusted
  )

# Add mapping
oct4_mapped <- oct4_raw %>%
  left_join(mapping, by = c("refseq_id" = "refseq_mrna"))

# Count how many did NOT map
unmapped_count <- sum(is.na(oct4_mapped$ensembl_gene_id))

# Of those that mapped, count how many have NA log2FC or padj (i.e. low expression)
low_expr_count <- oct4_mapped %>%
  filter(!is.na(ensembl_gene_id)) %>%
  filter(is.na(oct4_log2fc) | is.na(oct4_padj)) %>%
  nrow()

# Total mapped
mapped_total <- sum(!is.na(oct4_mapped$ensembl_gene_id))

cat("Total DESeq2 entries:", nrow(oct4_mapped), "\n")
cat("Unmapped (no Ensembl ID):", unmapped_count, "\n")
cat("Mapped but low expression (NA values):", low_expr_count, "\n")
cat("Mapped and usable:", mapped_total - low_expr_count, "\n")






bdp_pair_coverage <- bdp_expr %>%
  mutate(
    watson_status = case_when(
      is.na(watson_ensembl_id) ~ "Unmapped",
      is.na(oct4_watson_log2fc) | is.na(oct4_watson_padj) ~ "Low expression",
      TRUE ~ "OK"
    ),
    crick_status = case_when(
      is.na(crick_ensembl_id) ~ "Unmapped",
      is.na(oct4_crick_log2fc) | is.na(oct4_crick_padj) ~ "Low expression",
      TRUE ~ "OK"
    ),
    bdp_coverage = case_when(
      watson_status == "OK" & crick_status == "OK" ~ "Both OK",
      xor(watson_status == "OK", crick_status == "OK") ~ "One OK",
      TRUE ~ "Neither OK"
    )
  )

# Count BDPs by coverage category
bdp_pair_coverage_summary <- bdp_pair_coverage %>%
  count(bdp_coverage) %>%
  arrange(desc(n))

print(bdp_pair_coverage_summary)




bdp_mapped_only <- bdp_expr %>%
  filter(!is.na(watson_ensembl_id) & !is.na(crick_ensembl_id))

joined_df <- oct4_results %>%
  left_join(mapping, by = c("refseq_id" = "refseq_mrna"))

# Create gene status table from DESeq2 + mapping
gene_status <- joined_df %>%
  mutate(status = case_when(
    is.na(ensembl_gene_id.y) ~ "Unmapped",
    TRUE ~ "OK"
  )) %>%
  dplyr::select(refseq_id, ensembl_gene_id = ensembl_gene_id.y, status)

# Merge statuses into BDP expression data
bdp_qc <- bdp_expr %>%
  left_join(gene_status %>% dplyr::select(ensembl_gene_id, status) %>%
              rename(watson_ensembl_id = ensembl_gene_id,
                     watson_status = status),
            by = "watson_ensembl_id") %>%
  left_join(gene_status %>% dplyr::select(ensembl_gene_id, status) %>%
              rename(crick_ensembl_id = ensembl_gene_id,
                     crick_status = status),
            by = "crick_ensembl_id")

# Apply filter:
bdp_filtered <- bdp_qc %>%
  filter(
    watson_status != "Unmapped" & crick_status != "Unmapped",  # Drop any with annotation issues
    watson_status == "OK" | crick_status == "OK"               # Keep only if at least one tested
  )


table(bdp_filtered$watson_status, bdp_filtered$crick_status)

cat("Final filtered BDPs:", nrow(bdp_filtered), "\n")
write.csv(bdp_filtered, "outputs/Filtered_BDP_Coordination_Results.csv", row.names = FALSE)

bdp_filtered_summary <- bdp_filtered %>%
  dplyr::select(
    BDP_ID,
    watson_gene_name, watson_ensembl_id,
    watson_log2fc, watson_padj,
    crick_gene_name, crick_ensembl_id,
    crick_log2fc, crick_padj,
    brg1_coordination, oct4_coordination
  )

write.csv(bdp_filtered_summary, "outputs/Filtered_BDP_Coordination_Results_summary.csv", row.names = FALSE)








# 1. All BDPs with both Watson and Crick Ensembl IDs present
bdps_with_both_ids <- bdp %>%
  filter(!is.na(watson_ensembl_id) & !is.na(crick_ensembl_id))

# 2. Define detected Ensembl IDs (i.e., those that appear in the DESeq2 data)
detected_genes <- unique(c(brg1_clean$ensembl_gene_id, oct4_clean$ensembl_gene_id))

# 3. Identify BDPs that are *dropped* because one or both genes aren't in detected_genes
bdps_dropped <- bdps_with_both_ids %>%
  filter(!(watson_ensembl_id %in% detected_genes & crick_ensembl_id %in% detected_genes))

# 4. How many are dropped?
n_dropped <- nrow(bdps_dropped)
cat("BDPs with both IDs present but dropped due to RNA-seq filtering: ", n_dropped, "\n")



dropped_ids <- unique(c(bdps_dropped$watson_ensembl_id, bdps_dropped$crick_ensembl_id))
dropped_ids <- dropped_ids[!is.na(dropped_ids)]


n_matched_genes <- sum(dropped_ids %in% detected_genes)

n_total_genes <- length(dropped_ids)

cat("Out of", n_total_genes, "unique gene IDs in dropped BDPs,", n_matched_genes, "are present in RNA-seq detected_genes.\n")


