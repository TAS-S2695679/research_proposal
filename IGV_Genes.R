library(dplyr)
library(readr)

# Load your corrected coordination file
coord <- read_csv("BDP_coordination_results.csv")  # adjust path as needed

# Step 1: Filter interesting coordination categories
interesting_classes <- c("Asymmetric", "Opposite", "Synchronised")

bdp_filtered <- coord %>%
  filter(oct4_coordination %in% interesting_classes)

# Step 2: Create long format with both strands
watson <- bdp_filtered %>%
  transmute(
    BDP_ID,
    gene = watson_ensembl_id,
    strand = "Watson",
    log2fc = oct4_watson_log2fc,
    partner_gene = crick_ensembl_id,
    partner_log2fc = oct4_crick_log2fc,
    coordination = oct4_coordination
  )

crick <- bdp_filtered %>%
  transmute(
    BDP_ID,
    gene = crick_ensembl_id,
    strand = "Crick",
    log2fc = oct4_crick_log2fc,
    partner_gene = watson_ensembl_id,
    partner_log2fc = oct4_watson_log2fc,
    coordination = oct4_coordination
  )

bdp_long <- bind_rows(watson, crick)

# Step 3: Flag genes with strong change
bdp_long <- bdp_long %>%
  mutate(
    abs_log2fc = abs(log2fc),
    is_strong_change = abs_log2fc > 1
  )

# Step 4: Filter to keep only BDPs where at least one gene is strongly changed
high_interest <- bdp_long %>%
  group_by(BDP_ID) %>%
  filter(any(is_strong_change)) %>%
  ungroup()

# Step 5: Save output
write_csv(high_interest, "outputs/high_interest_bdp_genes_for_IGV.csv")



ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

all_ids <- unique(c(high_interest$gene, high_interest$partner_gene))

mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = all_ids,
  mart = ensembl
)

bdp_named <- high_interest %>%
  dplyr::left_join(mapping, by = c("gene" = "ensembl_gene_id")) %>%
  dplyr::rename(gene_symbol = external_gene_name) %>%
  dplyr::left_join(mapping, by = c("partner_gene" = "ensembl_gene_id")) %>%
  dplyr::rename(partner_symbol = external_gene_name)

# Step 6: Save final annotated table
write_csv(bdp_named, "outputs/high_interest_bdp_genes_with_names.csv")

