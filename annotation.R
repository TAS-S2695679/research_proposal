library(dplyr)
library(readr)
library(biomaRt)
library(httr)
library(jsonlite)
library(tidyr)
library(progress)

# Step 1: Load the cleaned BDP pair table
bdp_pairs <- read_csv("BDP_Cleaned_Pairs.csv") 

# Step 2: Set up Ensembl BioMart (mouse genes)
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# Step 3: Define attributes to pull from gene-level annotations
gene_attributes <- c("ensembl_gene_id",
                     "external_gene_name",
                     "gene_biotype",
                     "description",
                     "refseq_ncrna",
                     "uniprotswissprot")

# Step 4: Get unique gene IDs from both strands
all_gene_ids <- unique(c(bdp_pairs$watson_ensembl_id, bdp_pairs$crick_ensembl_id))
all_gene_ids <- na.omit(all_gene_ids)

# Step 5: Query gene-level annotations
# Core attributes that always exist
essential_attrs <- c("ensembl_gene_id", "external_gene_name", "gene_biotype", "description")
essential_annotations <- getBM(
  attributes = essential_attrs,
  filters = "ensembl_gene_id",
  values = all_gene_ids,
  mart = ensembl
)

# Optional annotations (may be missing for some genes)
refseq_annotations <- getBM(
  attributes = c("ensembl_gene_id", "refseq_ncrna"),
  filters = "ensembl_gene_id",
  values = all_gene_ids,
  mart = ensembl
)

uniprot_annotations <- getBM(
  attributes = c("ensembl_gene_id", "uniprotswissprot"),
  filters = "ensembl_gene_id",
  values = all_gene_ids,
  mart = ensembl
)

uniprot_annotations <- uniprot_annotations %>%
  group_by(ensembl_gene_id) %>%
  slice(1) %>%
  ungroup()

# Combine all annotations (keeps all genes even if some fields are missing)
gene_annotations <- essential_annotations %>%
  left_join(refseq_annotations, by = "ensembl_gene_id") %>%
  left_join(uniprot_annotations, by = "ensembl_gene_id")

# Step 5a: Deduplicate (1 row per gene ID)
gene_annotations <- gene_annotations %>%
  group_by(ensembl_gene_id) %>%
  slice(1) %>%
  ungroup()

# Step 5b: Get transcript-level RNAcentral IDs
transcript_lookup <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "rnacentral"),
                           filters = "ensembl_gene_id",
                           values = all_gene_ids,
                           mart = ensembl)

# Step 5c: Deduplicate RNAcentral (1 transcript per gene if available)
transcript_lookup_clean <- transcript_lookup %>%
  filter(!is.na(rnacentral)) %>%
  group_by(ensembl_gene_id) %>%
  slice(1) %>%
  ungroup()

# Step 6: Merge gene and transcript annotations
full_annot <- gene_annotations %>%
  left_join(transcript_lookup_clean, by = "ensembl_gene_id")

# Step 7: Prepare Watson and Crick annotation tables
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

# Step 8: Merge into the BDP dataset
bdp_annotated <- bdp_pairs %>%
  left_join(watson_info, by = "watson_ensembl_id") %>%
  left_join(crick_info, by = "crick_ensembl_id")


#UNIPROT ANNOTATIONS
#head(na.omit(unique(c(bdp_annotated$watson_uniprot, bdp_annotated$crick_uniprot))), 500)

# Step 9: Combine all UniProt IDs and remove NA
all_uniprot_ids <- unique(c(bdp_annotated$watson_uniprot, bdp_annotated$crick_uniprot))
all_uniprot_ids <- na.omit(all_uniprot_ids)
all_uniprot_ids <- all_uniprot_ids[!is.na(all_uniprot_ids) & all_uniprot_ids != ""]


# Step 10: Function to query UniProt API for one ID
get_uniprot_keywords <- function(id) {
  
  url <- paste0("https://rest.uniprot.org/uniprotkb/", id, ".json")
  res <- tryCatch(GET(url), error = function(e) NULL)
  
  if (!is.null(res) && status_code(res) == 200) {
    data <- httr::content(res, encoding = "UTF-8")
    
    if (!is.null(data$keywords) && length(data$keywords) > 0) {
      valid_keywords <- Filter(function(k) !is.null(k$name) && nzchar(k$name), data$keywords)
      
      if (length(valid_keywords) > 0) {
        kw <- vapply(valid_keywords, function(k) k$name, character(1))
        return(data.frame(uniprot_id = rep(id, length(kw)), keyword = kw))
      }
    }
  }
  
  return(data.frame(uniprot_id = id, keyword = NA_character_))
}

# Step 11: Loop through IDs and fetch keywords
pb <- progress_bar$new(
  format = "Querying UniProt [:bar] :percent (:current/:total)",
  total = length(all_uniprot_ids), clear = FALSE, width = 60
)

keyword_list <- vector("list", length(all_uniprot_ids))

for (i in seq_along(all_uniprot_ids)) {
  keyword_list[[i]] <- get_uniprot_keywords(all_uniprot_ids[i])
  pb$tick()
  Sys.sleep(0.2)  # avoid rate limiting
}

# Step 12: Collapse multiple keyword rows into one per UniProt ID
uniprot_keywords <- bind_rows(keyword_list)

uniprot_keywords_collapsed <- uniprot_keywords %>%
  group_by(uniprot_id) %>%
  summarise(keyword = paste(na.omit(keyword), collapse = "; ")) %>%
  ungroup()

# Step 13: Join keywords back to Watson and Crick
watson_keywords <- uniprot_keywords_collapsed %>%
  rename(watson_uniprot = uniprot_id, watson_keyword = keyword)

crick_keywords <- uniprot_keywords_collapsed %>%
  rename(crick_uniprot = uniprot_id, crick_keyword = keyword)

bdp_annotated <- bdp_annotated %>%
  left_join(watson_keywords, by = "watson_uniprot") %>%
  left_join(crick_keywords, by = "crick_uniprot")


# Step 14: Extract all transcript IDs from Watson and Crick
transcript_ids <- unique(c(bdp_annotated$watson_transcript_id, bdp_annotated$crick_transcript_id))
transcript_ids <- na.omit(transcript_ids)

# Step 15: Query GENCODE transcript-level annotations via Ensembl BioMart
gencode_attrs <- c(
  "ensembl_transcript_id",
  "transcript_biotype",
  "transcript_gencode_basic",
  "transcript_tsl",
  "transcript_appris",
  "transcript_is_canonical"
)

gencode_data <- getBM(
  attributes = gencode_attrs,
  filters = "ensembl_transcript_id",
  values = transcript_ids,
  mart = ensembl
)

# Step 16: Clean and prioritise GENCODE entries (canonical, basic, TSL)
gencode_data_clean <- gencode_data %>%
  arrange(desc(transcript_is_canonical), desc(transcript_gencode_basic == "GENCODE basic"), transcript_tsl) %>%
  group_by(ensembl_transcript_id) %>%
  slice(1) %>%
  ungroup()

# Step 17: Prepare strand-specific GENCODE data
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

# Step 18: Merge GENCODE annotations into main table
bdp_annotated <- bdp_annotated %>%
  left_join(watson_gencode, by = "watson_transcript_id") %>%
  left_join(crick_gencode, by = "crick_transcript_id")

# Reorder columns to improve readability
core_cols <- c("BDP_ID", "Chromosome", "watson_start", "watson_end", "crick_start", "crick_end")

watson_cols <- names(bdp_annotated)[startsWith(names(bdp_annotated), "watson_") & !(names(bdp_annotated) %in% core_cols)]
crick_cols <- names(bdp_annotated)[startsWith(names(bdp_annotated), "crick_") & !(names(bdp_annotated) %in% core_cols)]

updated_col_order <- c(core_cols, sort(watson_cols), sort(crick_cols))
bdp_annotated <- bdp_annotated[, updated_col_order]

#Write Final Fully Annotated Table
write_csv(bdp_annotated, "BDP_Fully_Annotated.csv")
