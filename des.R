library(DESeq2)

raw_counts <- read.delim("data/GSE87821_nucRNAseq_ReadCount_Quantitation.txt.gz", header = TRUE)

# Read counts (already done)
counts <- raw_counts

# Clean up
rownames(counts) <- counts$ID
counts_matrix <- counts[, grepl("ReadCount", colnames(counts))]

# Create sample metadata
samples <- colnames(counts_matrix)

# Extract condition and cell line from sample names
condition <- ifelse(grepl("UNT", samples), "UNT",
                    ifelse(grepl("DOX", samples), "DOX", "TAM"))

cell_line <- ifelse(grepl("ZHBTC4", samples), "ZHBTC4", "BRG1fl")

# Now build the data frame
sample_info <- data.frame(
  sample = samples,
  condition = factor(condition, levels = c("UNT", "DOX", "TAM")),
  cell_line = factor(cell_line)
)
rownames(sample_info) <- sample_info$sample


dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = sample_info,
                              design = ~ cell_line + condition)

dds <- DESeq(dds)


# Subset to BRG1fl samples
dds_brg1 <- dds[, dds$cell_line == "BRG1fl"]
dds_brg1$condition <- droplevels(dds_brg1$condition)
design(dds_brg1) <- ~ condition

# Run DESeq
dds_brg1 <- DESeq(dds_brg1)

# Get TAM vs UNT results (BRG1 depletion)
res_brg1 <- results(dds_brg1, contrast = c("condition", "TAM", "UNT"))

# Convert to data frame and add RefSeq ID column
brg1_results <- as.data.frame(res_brg1) %>%
  rename(brg1_log2fc = log2FoldChange, brg1_padj = padj)


# Subset to ZHBTC4 samples
dds_oct4 <- dds[, dds$cell_line == "ZHBTC4"]
dds_oct4$condition <- droplevels(dds_oct4$condition)
design(dds_oct4) <- ~ condition

# Run DESeq
dds_oct4 <- DESeq(dds_oct4)

# Get DOX vs UNT results (Oct4 depletion)
res_oct4 <- results(dds_oct4, contrast = c("condition", "DOX", "UNT"))

# Convert to data frame and add RefSeq ID column
oct4_results <- as.data.frame(res_oct4) %>%
  tibble::rownames_to_column("refseq_id") %>%
  rename(oct4_log2fc = log2FoldChange, oct4_padj = padj)


write.csv(brg1_results, "data/BRG1_DESeq2_results.csv", row.names = FALSE)
write.csv(oct4_results, "data/Oct4_DESeq2_results.csv", row.names = FALSE)







raw_counts$baseMean <- rowMeans(raw_counts[, grepl("ReadCount", colnames(raw_counts))])

# Filter genes that are actually detected (baseMean > 0, or use >5 if stricter)
expr_detected <- raw_counts %>%
  filter(baseMean > 0)

# Get RefSeq IDs
detected_refseq_ids <- expr_detected$ID

# Map those to Ensembl using biomaRt
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

refseq_to_ensembl <- getBM(
  attributes = c("refseq_mrna", "ensembl_gene_id"),
  filters = "refseq_mrna",
  values = detected_refseq_ids,
  mart = ensembl
) %>%
  filter(!is.na(ensembl_gene_id)) %>%
  distinct()

# Final list of expressed Ensembl IDs
detected_genes <- unique(refseq_to_ensembl$ensembl_gene_id)

bdp_expr <- bdp_expr %>%
  filter(watson_ensembl_id %in% detected_genes,
         crick_ensembl_id %in% detected_genes)
