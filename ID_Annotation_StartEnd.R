library(biomaRt)
library(dplyr)
library(readr)

# Load and clean dataset
bdp <- read_csv("BDP_Cleaned_Pairs.csv") %>%
  mutate(Chromosome = gsub("^chr", "", Chromosome))  # standardise chromosome format

# Set up Ensembl BioMart
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# Extract rows with missing Watson and Crick Ensembl IDs but valid coordinates
watson_missing <- bdp %>%
  filter(is.na(watson_ensembl_id), !is.na(watson_start), !is.na(watson_end)) %>%
  mutate(strand = 1, start = watson_start, end = watson_end) %>%
  select(Chromosome, start, end, strand)

crick_missing <- bdp %>%
  filter(is.na(crick_ensembl_id), !is.na(crick_start), !is.na(crick_end)) %>%
  mutate(strand = -1, start = crick_start, end = crick_end) %>%
  select(Chromosome, start, end, strand)

# Combine into one search table
rescue_set <- bind_rows(watson_missing, crick_missing) %>%
  rename(chromosome = Chromosome)

# Prepare output dataframe
rescue_results <- data.frame(
  chromosome = character(),
  start = integer(),
  end = integer(),
  strand = integer(),
  ensembl_gene_id = character(),
  external_gene_name = character(),
  gene_biotype = character(),
  stringsAsFactors = FALSE
)

# Run coordinate-based queries (simple loop)
for (i in 1:nrow(rescue_set)) {
  cat("Query", i, "of", nrow(rescue_set), ":", 
      rescue_set$chromosome[i], rescue_set$start[i], rescue_set$end[i], "strand", rescue_set$strand[i], "\n")
  
  result <- tryCatch({
    getBM(
      attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
      filters = c("chromosome_name", "start", "end", "strand"),
      values = list(
        rescue_set$chromosome[i],
        rescue_set$start[i],
        rescue_set$end[i],
        rescue_set$strand[i]
      ),
      mart = ensembl
    )
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(result) && nrow(result) > 0) {
    first_match <- result[1, ]
    rescue_results <- bind_rows(rescue_results, tibble(
      chromosome = rescue_set$chromosome[i],
      start = rescue_set$start[i],
      end = rescue_set$end[i],
      strand = rescue_set$strand[i],
      ensembl_gene_id = first_match$ensembl_gene_id,
      external_gene_name = first_match$external_gene_name,
      gene_biotype = first_match$gene_biotype
    ))
    print(first_match[1, ])
  } else {
    rescue_results <- bind_rows(rescue_results, tibble(
      chromosome = rescue_set$chromosome[i],
      start = rescue_set$start[i],
      end = rescue_set$end[i],
      strand = rescue_set$strand[i],
      ensembl_gene_id = NA,
      external_gene_name = NA,
      gene_biotype = NA
    ))
    cat(" No gene found\n")
  }
}

print(sum(!is.na(rescue_results$ensembl_gene_id)))

# Save results
write_csv(rescue_results, "Rescued_Gene_Annotations.csv")
cat("\n Finished. Rescued annotations saved to 'Rescued_Gene_Annotations.csv'\n")
