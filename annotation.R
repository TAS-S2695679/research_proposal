# -----------------------------------------------------------------------------
# Annotate Bidirectional Promoters (BDPs) in Mouse ES Cells
# Integrates Ensembl gene & transcript info, UniProt keywords, and GENCODE metadata
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(biomaRt)
library(httr)
library(jsonlite)
library(tidyr)
library(progress)

# --------------------------- Step 1: Load BDP Pair Metadata ------------------
bdp_pairs <- read_csv("BDP_Cleaned_Pairs.csv")

# --------------------------- Step 2: Setup BioMart ---------------------------
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# --------------------------- Step 3: Gene-Level Annotation -------------------
gene_attributes <- c("ensembl_gene_id", "external_gene_name", "gene_biotype", 
                     "description", "refseq_ncrna", "uniprotswissprot")

all_gene_ids <- na.omit(unique(c(bdp_pairs$watson_ensembl_id, bdp_pairs$crick_ensembl_id)))

gene_data <- getBM(attributes = gene_attributes, filters = "ensembl_gene_id", 
                   values = all_gene_ids, mart = ensembl)

gene_data <- gene_data %>%
  group_by(ensembl_gene_id) %>%
  slice(1) %>%
  ungroup()

# --------------------------- Step 4: RNAcentral Transcript IDs ---------------
transcript_data <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "rnacentral"),
                         filters = "ensembl_gene_id",
                         values = all_gene_ids,
                         mart = ensembl) %>%
  filter(!is.na(rnacentral)) %>%
  group_by(ensembl_gene_id) %>%
  slice(1) %>%
  ungroup()

# --------------------------- Step 5: Merge Gene + Transcript -----------------
full_annot <- gene_data %>%
  left_join(transcript_data, by = "ensembl_gene_id")

# --------------------------- Step 6: Format Strand-Specific Data -------------
watson_info <- full_annot %>%
  rename(watson_ensembl_id = ensembl_gene_id,
         watson_gene_name_updated = external_gene_name,
         watson_biotype = gene_biotype,
         watson_description = description,
         watson_refseq = refseq_ncrna,
         watson_uniprot = uniprotswissprot,
         watson_transcript_id = ensembl_transcript_id,
         watson_rnacentral_id = rnacentral)

crick_info <- full_annot %>%
  rename(crick_ensembl_id = ensembl_gene_id,
         crick_gene_name_updated = external_gene_name,
         crick_biotype = gene_biotype,
         crick_description = description,
         crick_refseq = refseq_ncrna,
         crick_uniprot = uniprotswissprot,
         crick_transcript_id = ensembl_transcript_id,
         crick_rnacentral_id = rnacentral)

bdp_annotated <- bdp_pairs %>%
  left_join(watson_info, by = "watson_ensembl_id") %>%
  left_join(crick_info, by = "crick_ensembl_id")

# --------------------------- Step 7: UniProt Keyword Annotation --------------
# Helper function
get_uniprot_keywords <- function(id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", id, ".json")
  res <- tryCatch(GET(url), error = function(e) NULL)
  if (!is.null(res) && status_code(res) == 200) {
    data <- httr::content(res, as = "parsed", encoding = "UTF-8")
    if (!is.null(data$keywords)) {
      kw <- vapply(data$keywords, \(k) k$name %||% "", character(1))
      kw <- kw[nzchar(kw)]
      return(data.frame(uniprot_id = id, keyword = paste(kw, collapse = "; ")))
    }
  }
  data.frame(uniprot_id = id, keyword = NA_character_)
}

# Query keywords
all_uniprot_ids <- na.omit(unique(c(bdp_annotated$watson_uniprot, bdp_annotated$crick_uniprot)))
pb <- progress_bar$new(format = "Querying UniProt [:bar] :percent", total = length(all_uniprot_ids))

keyword_list <- lapply(seq_along(all_uniprot_ids), function(i) {
  pb$tick()
  Sys.sleep(0.2)  # avoid rate limiting
  get_uniprot_keywords(all_uniprot_ids[i])
})

uniprot_keywords <- bind_rows(keyword_list)

# Merge UniProt keywords
bdp_annotated <- bdp_annotated %>%
  left_join(uniprot_keywords %>% rename(watson_uniprot = uniprot_id, watson_keyword = keyword),
            by = "watson_uniprot") %>%
  left_join(uniprot_keywords %>% rename(crick_uniprot = uniprot_id, crick_keyword = keyword),
            by = "crick_uniprot")

# --------------------------- Step 8: GENCODE Annotation ----------------------
transcript_ids <- na.omit(unique(c(bdp_annotated$watson_transcript_id, bdp_annotated$crick_transcript_id)))

gencode_attrs <- c("ensembl_transcript_id", "transcript_biotype", "transcript_gencode_basic",
                   "transcript_tsl", "transcript_appris", "transcript_is_canonical")

gencode_data <- getBM(attributes = gencode_attrs,
                      filters = "ensembl_transcript_id",
                      values = transcript_ids,
                      mart = ensembl)

# Prioritise canonical or basic
gencode_data_clean <- gencode_data %>%
  arrange(desc(transcript_is_canonical), desc(transcript_gencode_basic == "GENCODE basic"), transcript_tsl) %>%
  group_by(ensembl_transcript_id) %>%
  slice(1) %>%
  ungroup()

# Format strand-specific
watson_gencode <- gencode_data_clean %>%
  filter(ensembl_transcript_id %in% bdp_annotated$watson_transcript_id) %>%
  rename(watson_transcript_id = ensembl_transcript_id,
         watson_biotype_gencode = transcript_biotype,
         watson_gencode_basic = transcript_gencode_basic,
         watson_tsl = transcript_tsl,
         watson_appris = transcript_appris,
         watson_canonical = transcript_is_canonical)

crick_gencode <- gencode_data_clean %>%
  filter(ensembl_transcript_id %in% bdp_annotated$crick_transcript_id) %>%
  rename(crick_transcript_id = ensembl_transcript_id,
         crick_biotype_gencode = transcript_biotype,
         crick_gencode_basic = transcript_gencode_basic,
         crick_tsl = transcript_tsl,
         crick_appris = transcript_appris,
         crick_canonical = transcript_is_canonical)

# Merge
bdp_annotated <- bdp_annotated %>%
  left_join(watson_gencode, by = "watson_transcript_id") %>%
  left_join(crick_gencode, by = "crick_transcript_id")

# --------------------------- Step 9: Final Column Ordering -------------------
core_cols <- c("BDP_ID", "Chromosome", "watson_start", "watson_end", "crick_start", "crick_end")
watson_cols <- names(bdp_annotated)[startsWith(names(bdp_annotated), "watson_") & !(names(bdp_annotated) %in% core_cols)]
crick_cols <- names(bdp_annotated)[startsWith(names(bdp_annotated), "crick_") & !(names(bdp_annotated) %in% core_cols)]
ordered_cols <- c(core_cols, sort(watson_cols), sort(crick_cols))
bdp_annotated <- bdp_annotated[, ordered_cols]

# --------------------------- Step 10: Export Final Dataset -------------------
write_csv(bdp_annotated, "BDP_Fully_Annotated.csv")
