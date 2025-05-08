# Load required libraries
library(dplyr)
library(readr)
library(biomaRt)

# Step 1: Load your existing annotated file
bdp <- read_csv("BDP_Fully_Annotated.csv")

# Step 2: Extract unique transcript IDs from both strands
transcript_ids <- unique(c(bdp$watson_transcript_id, bdp$crick_transcript_id))
transcript_ids <- na.omit(transcript_ids)
#transcript_ids

# Step 3: Connect to BioMart (Ensembl)
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# Step 4: Define GENCODE-relevant attributes
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
  values = transcript_ids,  # from your BDP_Fully_Annotated.csv
  mart = ensembl
)

# Step 6: Deduplicate per transcript (if multiple rows exist)
gencode_data_clean <- gencode_data %>%
  arrange(transcript_is_canonical == 1, transcript_gencode_basic == 1, transcript_tsl) %>%
  group_by(ensembl_transcript_id) %>%
  slice(1) %>%
  ungroup()

# Step 7: Split into Watson and Crick versions
watson_gencode <- gencode_data_clean %>%
  filter(ensembl_transcript_id %in% bdp$watson_transcript_id) %>%
  rename(watson_transcript_id = ensembl_transcript_id,
         watson_biotype_gencode = transcript_biotype,
         watson_gencode_basic = transcript_gencode_basic,
         watson_tsl = transcript_tsl,
         watson_appris = transcript_appris,
         watson_canonical = transcript_is_canonical)

crick_gencode <- gencode_data_clean %>%
  filter(ensembl_transcript_id %in% bdp$crick_transcript_id) %>%
  rename(crick_transcript_id = ensembl_transcript_id,
         crick_biotype_gencode = transcript_biotype,
         crick_gencode_basic = transcript_gencode_basic,
         crick_tsl = transcript_tsl,
         crick_appris = transcript_appris,
         crick_canonical = transcript_is_canonical)

# Step 8: Merge back into the main dataset
bdp_gencode <- bdp %>%
  left_join(watson_gencode, by = "watson_transcript_id") %>%
  left_join(crick_gencode, by = "crick_transcript_id")

# Step 9: Save results
write_csv(bdp_gencode, "BDP_Annotated_With_GENCODE.csv")
