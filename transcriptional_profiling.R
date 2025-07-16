library(dplyr)
library(biomaRt)
library(readr)
library(ggplot2)
library(org.Mm.eg.db)
library(AnnotationDbi)

brg1_results <- read_csv("data/BRG1_DESeq2_results.csv")
oct4_results <- read_csv("data/Oct4_DESeq2_results.csv")


all_refseq_ids <- unique(c(brg1_results$refseq_id, oct4_results$refseq_id))

all_refseq_ids_no_version <- sub("\\..*", "", all_refseq_ids)


refseq_to_ensembl <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = all_refseq_ids_no_version,
  keytype = "REFSEQ",
  columns = c("ENSEMBL", "SYMBOL")
)

mapping <- unique(refseq_to_ensembl[, c("REFSEQ", "ENSEMBL")])
colnames(mapping)[colnames(mapping) == "REFSEQ"] <- "refseq_id"
colnames(mapping)[colnames(mapping) == "ENSEMBL"] <- "ensembl_gene_id"




#ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = 101)

#mapping <- getBM(
 # attributes = c("refseq_mrna", "ensembl_gene_id"),
  #filters = "refseq_mrna",
  #values = all_refseq_ids,
  #mart = ensembl
#)


all_bdp_ids <- unique(c(bdp$watson_ensembl_id, bdp$crick_ensembl_id))

matched <- all_bdp_ids %in% mapping$ensembl_gene_id

table(matched)

valid_pairs <- bdp %>%
  filter(watson_ensembl_id %in% detected_genes,
         crick_ensembl_id %in% detected_genes)

nrow(valid_pairs)


initial_valid_pairs <- bdp %>%
  filter(
    watson_ensembl_id %in% detected_genes,
    crick_ensembl_id %in% detected_genes
  )

nrow(initial_valid_pairs)




# Merge mapping
brg1_results <- brg1_results %>%
  dplyr::left_join(mapping, by = c("refseq_id" = "refseq_id")) %>%
  dplyr::filter(!is.na(ensembl_gene_id))

brg1_results$ensembl_gene_id <- as.character(brg1_results$ensembl_gene_id)


oct4_results <- oct4_results %>%
  left_join(mapping, by = c("refseq_id" = "refseq_id")) %>%
  mutate(status = case_when(
    is.na(ensembl_gene_id) ~ "Unmapped",
    is.na(oct4_log2fc) ~ "Low expression",
    TRUE ~ "OK"
  ))

oct4_results$ensembl_gene_id <- as.character(oct4_results$ensembl_gene_id)


head(brg1_results)
head(oct4_results)

# Filter BRG1 results to 1 row per Ensembl ID
brg1_clean <- brg1_results %>%
  group_by(ensembl_gene_id) %>%
  slice_min(order_by = brg1_padj, n = 1, with_ties = FALSE) %>%
  ungroup()

# Filter oct4 results to 1 row per Ensembl ID

oct4_clean <- oct4_results %>%
  group_by(ensembl_gene_id) %>%
  slice_min(order_by = oct4_padj, n = 1, with_ties = FALSE) %>%
  ungroup()

#--------- Now to merge and classigy BDP by Coordination ----------
bdp <- read_csv("BDP_Fully_Annotated.csv")

bdp_expr <- bdp %>%
  left_join(
    brg1_clean %>%
      dplyr::select(ensembl_gene_id, brg1_log2fc, brg1_padj) %>%
      dplyr::rename(
        watson_log2fc = brg1_log2fc,
        watson_padj = brg1_padj
      ),
    by = c("watson_ensembl_id" = "ensembl_gene_id")
  ) %>%
  left_join(
    brg1_clean %>%
      dplyr::select(ensembl_gene_id, brg1_log2fc, brg1_padj) %>%
      dplyr::rename(
        crick_log2fc = brg1_log2fc,
        crick_padj = brg1_padj
      ),
    by = c("crick_ensembl_id" = "ensembl_gene_id")
  )
  
  # Step 2: Classify BRG1 coordination
  
  bdp_expr <- bdp_expr %>%
  mutate(
    watson_brg1_sig = !is.na(watson_log2fc) & !is.na(watson_padj) &
      abs(watson_log2fc) > 0.30 & watson_padj < 0.05,
    
    crick_brg1_sig = !is.na(crick_log2fc) & !is.na(crick_padj) &
      abs(crick_log2fc) > 0.30 & crick_padj < 0.05
  )

# Apply coordination classification using flags
  bdp_expr <- bdp_expr %>%
    mutate(
      brg1_coordination = case_when(
        watson_brg1_sig & crick_brg1_sig &
          sign(watson_log2fc) == sign(crick_log2fc) ~ "Synchronised",
        
        watson_brg1_sig & crick_brg1_sig &
          sign(watson_log2fc) != sign(crick_log2fc) ~ "Opposite",
        
        xor(watson_brg1_sig, crick_brg1_sig) ~ "Asymmetric",
        
        TRUE ~ "Unchanged"
      )
    ) %>%
    left_join(
      oct4_clean %>%
        dplyr::select(ensembl_gene_id, oct4_log2fc, oct4_padj) %>%
        dplyr::rename(
          oct4_watson_log2fc = oct4_log2fc,
          oct4_watson_padj = oct4_padj
        ),
      by = c("watson_ensembl_id" = "ensembl_gene_id")
    ) %>%
    left_join(
      oct4_clean %>%
        dplyr::select(ensembl_gene_id, oct4_log2fc, oct4_padj) %>%
        dplyr::rename(
          oct4_crick_log2fc = oct4_log2fc,
          oct4_crick_padj = oct4_padj
        ),
      by = c("crick_ensembl_id" = "ensembl_gene_id")
    )
  
  # Step 4: Classify Oct4 coordination
  # 1. Add significance flags for Watson and Crick strands
  bdp_expr <- bdp_expr %>%
    mutate(
      watson_sig = !is.na(oct4_watson_log2fc) & !is.na(oct4_watson_padj) &
        abs(oct4_watson_log2fc) > 0.30 & oct4_watson_padj < 0.05,
      
      crick_sig = !is.na(oct4_crick_log2fc) & !is.na(oct4_crick_padj) &
        abs(oct4_crick_log2fc) > 0.30 & oct4_crick_padj < 0.05
    )

# 2. Classify coordination
bdp_expr <- bdp_expr %>%
  mutate(
    oct4_coordination = case_when(
      watson_sig & crick_sig &
        sign(oct4_watson_log2fc) == sign(oct4_crick_log2fc) ~ "Synchronised",
      
      watson_sig & crick_sig &
        sign(oct4_watson_log2fc) != sign(oct4_crick_log2fc) ~ "Opposite",
      
      xor(watson_sig, crick_sig) ~ "Asymmetric",
      
      TRUE ~ "Unchanged"
    )
  )


table(bdp_expr$brg1_coordination)
table(bdp_expr$oct4_coordination)


bdp_expr <- bdp_expr %>%
  dplyr::select(BDP_ID,
         watson_ensembl_id, crick_ensembl_id,
         brg1_coordination, oct4_coordination, oct4_crick_log2fc, oct4_watson_log2fc)

detected_genes <- unique(c(oct4_clean$ensembl_gene_id, brg1_clean$ensembl_gene_id))

bdp_expr <- bdp_expr %>%
  filter(!is.na(watson_ensembl_id),
         !is.na(crick_ensembl_id),
         watson_ensembl_id %in% detected_genes,
         crick_ensembl_id %in% detected_genes)

write_csv(bdp_expr, "BDP_coordination_results.csv")


bdp_expr <- bdp_expr %>%
  filter(!is.na(watson_ensembl_id) & !is.na(crick_ensembl_id))
bdp_expr <- as_tibble(bdp_expr)

# ------------------ Coordination Summary ------------------

# BRG1 coordination counts
brg1_summary <- bdp_expr %>%
  dplyr::count(brg1_coordination, sort = TRUE) %>%
  dplyr::rename(Category = brg1_coordination, Count = n)

# OCT4 coordination counts
oct4_summary <- bdp_expr %>%
  dplyr::count(oct4_coordination, sort = TRUE) %>%
  dplyr::rename(Category = oct4_coordination, Count = n)

# Print summaries
cat("BRG1 Coordination Summary:\n")
print(brg1_summary)

cat("\nOCT4 Coordination Summary:\n")
print(oct4_summary)

oct4_results %>% filter(ensembl_gene_id == "ENSMUSG00000059552") %>% dplyr::select(oct4_log2fc, oct4_padj)

# ---- BDP gene coverage verification ----

# Get unique BDP gene IDs
bdp_ids <- unique(c(bdp$watson_ensembl_id, bdp$crick_ensembl_id)) %>%
  na.omit()

# Get unique expressed gene IDs from the Oct4 dataset
expr_ids <- unique(oct4_clean$ensembl_gene_id)

# Find how many BDP-linked genes are detected in Oct4 RNA-seq
expressed_bdp_genes <- intersect(bdp_ids, expr_ids)

# Print summary
cat("Total BDP-linked genes:", length(bdp_ids), "\n")
cat("BDP-linked genes detected in Oct4 dataset:", length(expressed_bdp_genes), "\n")

# Calculate and print percentage
coverage_rate <- length(expressed_bdp_genes) / length(bdp_ids) * 100
cat("Coverage rate:", round(coverage_rate, 1), "%\n")

#Of the 6,432 total BDP-linked genes, 3,897 (~61%) were detected in the Oct4 nuclear RNA-seq dataset (GSE87821). Coordination analysis was therefore limited to these expressed genes


example_refseqs <- c(
  "NM_001001333", "NM_001001445", "NM_001001565",
  "NM_001001809", "NM_001001883"
)

mapping %>%
  filter(refseq_mrna %in% example_refseqs)


oct4_synch_bdp <- bdp_expr %>%
  filter(oct4_coordination == "Synchronised") %>%
  dplyr::select(watson_gene_name, crick_gene_name,
                oct4_watson_log2fc, oct4_crick_log2fc,
                oct4_watson_padj, oct4_crick_padj)
oct4_synch_bdp


oct4_opposite_bdp <- bdp_expr %>%
       filter(oct4_coordination == "Opposite") %>%
       dplyr::select(
             BDP_ID,
             watson_gene_name, crick_gene_name,
             oct4_watson_log2fc, oct4_crick_log2fc,
             oct4_watson_padj, oct4_crick_padj)
 print(oct4_opposite_bdp)
 
 
 # Plot for BRG1 coordination
 bdp_expr %>%
   filter(!is.na(watson_log2fc) & !is.na(crick_log2fc)) %>%
   ggplot(aes(x = watson_log2fc, y = crick_log2fc, color = brg1_coordination)) +
   geom_point(alpha = 0.7, size = 2) +
   scale_color_manual(values = c(
     "Synchronised" = "forestgreen",
     "Opposite" = "firebrick",
     "Asymmetric" = "goldenrod",
     "Unchanged" = "grey70"
   )) +
   labs(
     title = "BRG1: Coordination of BDP Gene Pairs",
     x = "Watson log2FC",
     y = "Crick log2FC",
     color = "Coordination"
   ) +
   theme_minimal(base_size = 13)

 
 # Plot for OCT4 coordination
bdp_expr %>%
   filter(!is.na(oct4_watson_log2fc) & !is.na(oct4_crick_log2fc)) %>%
   ggplot(aes(x = oct4_watson_log2fc, y = oct4_crick_log2fc, color = oct4_coordination)) +
   geom_point(alpha = 0.7, size = 2) +
   scale_color_manual(values = c(
     "Synchronised" = "darkblue",
     "Opposite" = "orangered",
     "Asymmetric" = "darkorange",
     "Unchanged" = "grey70"
   )) +
   labs(
     title = "OCT4: Coordination of BDP Gene Pairs",
     x = "Watson log2FC",
     y = "Crick log2FC",
     color = "Coordination"
   ) +
   theme_minimal(base_size = 13)
 



go_overlap <- read_csv("outputs/GO_BP_BDP_TP53_overlap_genes.csv") %>%
  filter(is_in_p53_reactome == TRUE)

coordination <- read_csv("BDP_coordination_results.csv")

# Combine Watson + Crick Ensembl IDs with coordination
coord_long <- coordination %>%
  transmute(BDP_ID, gene = watson_ensembl_id, coordination = oct4_coordination) %>%
  bind_rows(
    coordination %>%
      transmute(BDP_ID, gene = crick_ensembl_id, coordination = oct4_coordination)
  ) %>%
  filter(!is.na(gene))

# Match GO-p53 overlap genes with coordination info
overlap_with_coordination <- go_overlap %>%
  left_join(coord_long, by = c("ensembl_id" = "gene")) %>%
  dplyr::distinct(ensembl_id, .keep_all = TRUE) %>%
  dplyr::count(coordination, sort = TRUE)



write_csv(overlap_with_coordination, "outputs/TP53_overlap_coordination_genes.csv")

print(overlap_with_coordination)









# Combine Watson and Crick into long format with labels
bdp_long <- bdp_expr %>%
  transmute(
    BDP_ID,
    gene_strand = "Watson",
    ensembl_id = watson_ensembl_id,
    log2fc = oct4_watson_log2fc,
    padj = oct4_watson_padj
  ) %>%
  bind_rows(
    bdp_expr %>%
      transmute(
        BDP_ID,
        gene_strand = "Crick",
        ensembl_id = crick_ensembl_id,
        log2fc = oct4_crick_log2fc,
        padj = oct4_crick_padj
      )
  )

# Classify reason for NA
bdp_long <- bdp_long %>%
  mutate(
    na_reason = case_when(
      is.na(ensembl_id) ~ "Unmapped (no Ensembl ID)",
      !is.na(ensembl_id) & (is.na(log2fc) | is.na(padj)) ~ "Low expression (not tested)",
      TRUE ~ "OK"
    )
  )

# Summarise
bdp_na_summary <- bdp_long %>%
  count(na_reason) %>%
  arrange(desc(n))

print(bdp_na_summary)

