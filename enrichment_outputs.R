# -------------------- Setup --------------------
library(clusterProfiler)
library(ggplot2)
library(DOSE)
library(enrichplot)
library(dplyr)
library(readr)
library(stringr)
library(forcats)

# -------------------- Load Results --------------------
go_bp <- read_csv("GO_combined_custom.csv")        # Output from enrichGO (BP terms, custom background)
kegg <- read_csv("kegg_enrich.csv")                # Output from enrichKEGG
reactome <- read_csv("reactome_enrich.csv")        # Output from enrichPathway

# -------------------- Output 1: Top GO Terms Table --------------------
top_go_terms <- go_bp %>%
  arrange(p.adjust) %>%
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  slice(1:10)

write_csv(top_go_terms, "outputs/enrichment/top_GO_terms_table.csv")
print(top_go_terms)

# -------------------- Output 2: Bar Plot of Top GO:BP Terms --------------------
# Take top 10 GO terms by adjusted p-value
top_go_plot <- go_bp %>%
  arrange(p.adjust) %>%
  slice(1:10) %>%
  mutate(Description = fct_reorder(Description, -log10(p.adjust)))

# Plot using ggplot2
pdf("outputs/enrichment/go_barplot.pdf", width = 8, height = 5)
print(
  ggplot(top_go_plot, aes(x = Description, y = -log10(p.adjust))) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Top Enriched GO:BP Terms", x = "", y = "-log10(FDR)") +
    theme_minimal(base_size = 12)
)
dev.off()

# -------------------- Output 3: Dot Plot --------------------
# Prepare top GO terms
top_go_dot <- go_bp %>%
  arrange(p.adjust) %>%
  slice(1:10) %>%
  mutate(
    Description = fct_reorder(Description, -p.adjust),
    logFDR = -log10(p.adjust),
    GeneRatio = Count / as.numeric(str_extract(BgRatio, "\\d+$"))  # extract denominator
  )

# Plot manually with ggplot2
pdf("outputs/enrichment/go_dotplot.pdf", width = 8, height = 6)
ggplot(top_go_dot, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = Count, color = logFDR)) +
  scale_color_viridis_c(option = "D") +
  labs(
    title = "GO:BP Dot Plot (Top 10 Enriched Terms)",
    x = "Gene Ratio",
    y = "",
    color = "-log10(FDR)",
    size = "Gene Count"
  ) +
  theme_minimal(base_size = 12)
dev.off()

# -------------------- Output 4: KEGG & Reactome Tables --------------------
top_kegg <- kegg %>%
  arrange(p.adjust) %>%
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  slice(1:10)

write_csv(top_kegg, "outputs/enrichment/top_KEGG_terms_table.csv")

top_reactome <- reactome %>%
  arrange(p.adjust) %>%
  dplyr::select(ID, Description, p.adjust, Count, geneID) %>%
  slice(1:10)

write_csv(top_reactome, "outputs/enrichment/top_Reactome_terms_table.csv")

# -------------------- Output 5: Biotype Composition (from separate analysis) --------------------
biotype_counts <- read_csv("outputs/enrichment/enrichment_biotype_summary.csv")

pdf("outputs/enrichment/biotype_barplot.pdf", width = 6, height = 4)
ggplot(biotype_counts, aes(x = reorder(biotype, -n), y = percent)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Biotype Composition of Enriched Genes", x = "Gene Biotype", y = "Percentage") +
  theme_minimal()
dev.off()

# -------------------- Output 6: Summary Stats --------------------
summary_stats <- tibble(
  GO_terms_tested = nrow(go_bp),
  GO_terms_significant = sum(go_bp$p.adjust < 0.05),
  KEGG_terms_significant = sum(kegg$p.adjust < 0.05),
  Reactome_terms_significant = sum(reactome$p.adjust < 0.05)
)

write_csv(summary_stats, "outputs/enrichment/enrichment_summary_stats.csv")
print(summary_stats)
