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
gene_annotations <- getBM(attributes = gene_attributes,
                          filters = "ensembl_gene_id",
                          values = all_gene_ids,
                          mart = ensembl)

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

# Step 9: Combine all UniProt IDs and remove NA
all_uniprot_ids <- unique(c(bdp_annotated$watson_uniprot, bdp_annotated$crick_uniprot))
all_uniprot_ids <- na.omit(all_uniprot_ids)
all_uniprot_ids <- all_uniprot_ids[grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$", all_uniprot_ids)]


# Step 10: Function to query UniProt API for one ID
get_uniprot_keywords <- function(id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", id, ".json")
  res <- tryCatch(GET(url), error = function(e) NULL)
  
  if (!is.null(res) && status_code(res) == 200) {
    data <- httr::content(res, as = "parsed", encoding = "UTF-8")
    
    if (!is.null(data$keywords) && length(data$keywords) > 0) {
      # Use human-readable labels
      valid_keywords <- Filter(function(k) !is.null(k$label) && nzchar(k$label), data$keywords)
      
      if (length(valid_keywords) > 0) {
        kw <- vapply(valid_keywords, function(k) k$label, character(1))
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
  Sys.sleep(0.2)  # Pause to avoid rate limiting
}

uniprot_keywords <- bind_rows(keyword_list)
# Step 12: Join keywords to Watson and Crick columns
watson_keywords <- uniprot_keywords %>%
  rename(watson_uniprot = uniprot_id, watson_keyword = keyword)

crick_keywords <- uniprot_keywords %>%
  rename(crick_uniprot = uniprot_id, crick_keyword = keyword)

bdp_annotated <- bdp_annotated %>%
  left_join(watson_keywords, by = "watson_uniprot") %>%
  left_join(crick_keywords, by = "crick_uniprot")

#Write Final Annotated Table
write_csv(bdp_annotated, "BDP_Fully_Annotated.csv")
