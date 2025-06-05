bdp <- read_csv("BDP_Fully_Annotated.csv")

expression_data <- read_csv("expression_dataset.csv") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Total number of active promoters (proxy = expressed genes)
total_active_promoters <- expression_data %>%
  filter(baseMean > 0) %>%
  pull(ensembl_gene_id) %>%
  unique()

# Number of promoters associated with BDPs = number of BDP rows
bdp_promoters_count <- nrow(bdp)

# Estimate number of non-BDP promoters
non_bdp_promoters_count <- length(total_active_promoters) - bdp_promoters_count

cat("Total active promoters:", length(total_active_promoters), "\n")
cat("BDP-associated promoters:", bdp_promoters_count, "\n")
cat("Non-BDP promoters:", non_bdp_promoters_count, "\n")
