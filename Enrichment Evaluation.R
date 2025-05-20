library(tidyverse)
library(readr)

# ----------------------- Load Data ----------------------------
enrich<- read_csv("GO_combined_custom.csv")
bdp <- read_csv("BDP_Fully_Annotated.csv")

# ----------------------- Extract Gene Symbols -----------------
# Combine all geneID entries and split on "/"
gene_symbols <- enrich$geneID %>%
  na.omit() %>%
  paste(collapse = "/") %>%
  str_split("/") %>%
  unlist() %>%
  str_trim() %>%
  toupper() %>%
  unique()

# ----------------------- Prepare Annotated Dataset ------------
# Extract gene name + biotype from both Watson and Crick
watson <- bdp %>%
  dplyr::select(gene_symbol = watson_gene_name, biotype = watson_biotype) %>%
  filter(!is.na(gene_symbol)) %>%
  mutate(gene_symbol = toupper(gene_symbol))

crick <- bdp %>%
  dplyr::select(gene_symbol = crick_gene_name, biotype = crick_biotype) %>%
  filter(!is.na(gene_symbol)) %>%
  mutate(gene_symbol = toupper(gene_symbol))

# Combine and de-duplicate
annotated_genes <- bind_rows(watson, crick) %>%
  distinct()

# ----------------------- Match & Count ------------------------
enriched_genes <- annotated_genes %>%
  filter(gene_symbol %in% gene_symbols)

biotype_summary <- enriched_genes %>%
  count(biotype, sort = TRUE) %>%
  mutate(percent = round(100 * n / sum(n), 2))

# ----------------------- Output Summary -----------------------
print(biotype_summary)
write_csv(biotype_summary, "outputs/enrichment/enrichment_biotype_summary.csv")

# ----------------------- Optional: Bar Plot -------------------
ggplot(biotype_summary, aes(x = reorder(biotype, -n), y = percent)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Biotype Composition of Enriched Genes",
    x = "Gene Biotype",
    y = "Percentage"
  ) +
  theme_minimal()
