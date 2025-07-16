# ---- Load libraries ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("rtracklayer", "GenomicRanges"))

library(rtracklayer)
library(GenomicRanges)
library(biomaRt)
library(dplyr)
library(tidyr)
library(stringr)

# ---- Set input/output ----
input_dir <- "ttseq"  # folder with all your .bw files (plus/minus per rep)
output_dir <- "ttseq_normalised_bw_stranded"
dir.create(output_dir, showWarnings = FALSE)

# ---- Define timepoints ----
timepoints <- c("0h", "3h", "6h", "9h", "12h", "15h")
strands <- c("plus", "minus")

# ---- Step 1: Calculate total signal per strand and timepoint ----
total_signal <- matrix(0, nrow = length(timepoints), ncol = 2,
                       dimnames = list(timepoints, strands))

for (tp in timepoints) {
  for (strand in strands) {
    cat("\n Processing:", tp, "| Strand:", strand, "\n")
    
    pattern <- paste0(tp, ".*", strand, ".*\\.bw$")
    files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
    
    strand_sum <- 0
    for (f in files) {
      cat(" Reading:", basename(f), "... ")
      bw <- import(f, as = "RleList")
      s <- sum(as.numeric(unlist(bw)))
      cat("Done. Signal =", round(s), "\n")
      strand_sum <- strand_sum + s
    }
    
    total_signal[tp, strand] <- strand_sum
    cat("  Total signal for", tp, strand, "=", round(strand_sum), "\n")
  }
}

# ---- Step 2: Compute scale factors per strand using 0h as reference ----
ref_plus <- total_signal["0h", "plus"]
ref_minus <- total_signal["0h", "minus"]

scale_factors <- total_signal
scale_factors[, "plus"] <- ref_plus / total_signal[, "plus"]
scale_factors[, "minus"] <- ref_minus / total_signal[, "minus"]

cat("\n Scale factors:\n")
print(round(scale_factors, 3))

# ---- Step 3: Merge, normalise, and export per strand and timepoint ----
for (tp in timepoints) {
  for (strand in strands) {
    cat("\n  Processing:", tp, "| Strand:", strand, "\n")
    
    pattern <- paste0(tp, ".*", strand, ".*\\.bw$")
    files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
    sf <- scale_factors[tp, strand]
    
    merged <- Reduce(`+`, lapply(files, function(f) import(f, as = "RleList")))
    norm <- merged * sf
    
    out_file <- file.path(output_dir, paste0("Oct4_", tp, "_", strand, "_scaled.bw"))
    export(norm, out_file)
    
    cat(" Saved:", basename(out_file), "\n")
  }
}

cat("\n All strand-specific normalised .bw files saved to:", output_dir, "\n")


#------Numerically quantify genes of interest ------
genes_of_interest = c("Trp53", "Wrap53", 'Fzd10', 'Fzd10os', 'Dnmt3l', 'Aire', 'Pnma2', 'Dpysl2')

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version =102)

gene_coords <- getBM(
  attributes = c("external_gene_name", "refseq_mrna", "ensembl_transcript_id",
                 "chromosome_name", "transcript_start", "transcript_end", "strand"),
  filters = "external_gene_name",
  values = genes_of_interest,
  mart = ensembl
)

gene_coords <- gene_coords %>%
  filter(!is.na(refseq_mrna) & refseq_mrna != "") %>%
  mutate(length = transcript_end - transcript_start) %>%
  group_by(external_gene_name) %>%
  slice_max(order_by = length, n = 1, with_ties = FALSE) %>%
  ungroup()


bw_dir <- "ttseq_normalised_bw_stranded"
timepoints <- c("0h", "3h", "6h", "9h", "12h", "15h")

# Convert gene coordinates into GRanges object
gene_gr <- GRanges(
  seqnames = paste0("chr", gene_coords$chromosome_name),
  ranges = IRanges(start = gene_coords$transcript_start, end = gene_coords$transcript_end),
  strand = ifelse(gene_coords$strand == 1, "+", "-"),
  gene = gene_coords$external_gene_name
)

# Empty list to store results
expression_results <- list()

# Loop through timepoints and strands
for (tp in timepoints) {
  for (strand in c("plus", "minus")) {
    message("Processing: ", tp, " | Strand: ", strand)
    
    # Load BigWig file
    bw_file <- file.path(bw_dir, paste0("Oct4_", tp, "_", strand, "_scaled.bw"))
    if (!file.exists(bw_file)) {
      warning("Missing file: ", bw_file)
      next
    }
    
    # Import BigWig
    bw_data <- import(bw_file, as = "RleList")
    
    # Filter genes that match the strand
    target_genes <- gene_gr[strand(gene_gr) == ifelse(strand == "plus", "+", "-")]
    
    # Extract signal
    for (gene in target_genes$gene) {
      region <- target_genes[target_genes$gene == gene]
      region_data <- import(bw_file, which = region)
      signal <- sum(score(region_data))
      expression_results[[length(expression_results) + 1]] <- data.frame(
        gene = gene,
        strand = as.character(strand(region)),
        timepoint = tp,
        signal = signal
      )
    }
  }
}

# Combine and tidy up
expression_df <- bind_rows(expression_results) %>%
  arrange(gene, timepoint)
expression_df$timepoint <- factor(expression_df$timepoint,
                                  levels = c("0h", "3h", "6h", "9h", "12h", "15h"),
                                  ordered = TRUE
)
expression_df <- expression_df %>%
  arrange(gene, timepoint)
print(expression_df)
