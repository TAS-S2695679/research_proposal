library(tidyverse)
library(readr)
library(ggplot2)
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

# ----------------------- Bar Plot -------------------
ggplot(biotype_summary, aes(x = reorder(biotype, -n), y = percent)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Biotype Composition of Enriched Genes",
    x = "Gene Biotype",
    y = "Percentage"
  ) +
  theme_minimal()

# ----------------------- Evaluate whether gene lengths should be included? -------------------
ggsave("outputs/gene_length_distribution.png", 
       plot = ggplot(combined, aes(x = gene_length, fill = set)) +
         geom_density(alpha = 0.5) +
         scale_x_log10() +
         labs(title = "Gene Length Distribution: Foreground vs Background",
              x = "Gene Length (log10)",
              y = "Density") +
         theme_minimal(),
       width = 8, height = 5, dpi = 300)


background_short <- control_pool %>%
  filter(gene_length < 2000 & !(ensembl_id %in% foreground_metadata$ensembl_id))

nrow(background_short)

library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
annotations <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = unique(background_short$ensembl_id),
  mart = mart
)
print(annotations)
#Key questions:
#Are they mostly non-coding (e.g., pseudogenes, lncRNAs)?
#Summary of gene biotypes
biotype_summary <- annotations %>% 
  count(gene_biotype, sort = TRUE) %>% 
  mutate(percentage = n/nrow(annotations)*100)

print(biotype_summary)

# Visualize
library(ggplot2)
ggplot(biotype_summary, aes(x = reorder(gene_biotype, -n), y = n)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Biotype Distribution in Short Background Genes",
       x = "Gene Biotype",
       y = "Count")

#Do they belong to specific families (e.g., histones, ribosome-related)?
# Extract common patterns from descriptions
# Extract common patterns from descriptions
family_keywords <- c("olfactory receptor", "taste receptor", "histone", "ribosomal", "keratin", "immunoglobulin")

annotations <- annotations %>% 
  mutate(
    gene_family = sapply(description, function(desc) {
      ifelse(any(sapply(family_keywords, grepl, x = desc, ignore.case = TRUE)),
             paste(family_keywords[sapply(family_keywords, grepl, x = desc, ignore.case = TRUE)], collapse = "; "),
             "Other")
    }
    ))
    
    # Count gene families
    family_counts <- annotations %>% 
      count(gene_family, sort = TRUE)
    print(family_counts)
#Are they poorly annotated or lowly expressed?
poor_annotations <- annotations %>% 
      filter(gene_biotype == "unknown" | 
               is.na(external_gene_name) | 
               grepl("^RIK.*|^Gm.*", external_gene_name))
    
    print(paste("Number of poorly annotated genes:", nrow(poor_annotations)))
    print(table(poor_annotations$gene_biotype))
        
    
    

foreground_short <- foreground_metadata %>% filter(gene_length < 2000)
nrow(foreground_short)

library(clusterProfiler)
enrich_result <- enrichGO(
  gene = background_short$ensembl_id,
  OrgDb = "org.Mm.eg.db",
  keyType = "ENSEMBL",
  ont = "BP"
)
dotplot(enrich_result)

foreground_enrich <- enrichGO(
  gene = foreground_short$ensembl_id,
  OrgDb = "org.Mm.eg.db",
  keyType = "ENSEMBL",
  ont = "BP"
)
dotplot(foreground_enrich)
