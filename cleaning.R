library(dplyr)
library(readr)
library(stringr)

# Step 1: Load the Catalogue
catalogue <- read_csv("Catalogue.csv")

# Step 2: Filter to transcript-level entries only
transcripts <- catalogue %>%
  filter(Feature == "transcript")

# Step 3: Extract base BDP ID by removing '_2' or '_3' suffix from Annotation
transcripts <- transcripts %>%
  mutate(BDP_ID = str_replace(Annotation, "_[23]$", ""))

# Step 4: Separate Watson (+ strand) and Crick (− strand) entries
watson <- transcripts %>%
  filter(Strand_Pos == "+") %>%
  group_by(BDP_ID) %>%
  slice(1) %>%
  ungroup()

crick <- transcripts %>%
  filter(Strand_Neg == "-") %>%
  group_by(BDP_ID) %>%
  slice(1) %>%
  ungroup()

# Step 5: Join Watson and Crick into BDP pairs
bdp_pairs <- inner_join(watson, crick, by = "BDP_ID", suffix = c("_watson", "_crick"))

# Step 6: Select and rename useful fields
bdp_clean <- bdp_pairs %>%
  transmute(
    BDP_ID,
    Chromosome = Chromosome_watson, #Keeping it as singular as its shared between both strands
    watson_ensembl_id = ensembl_id_watson,
    watson_gene_name = external_gene_name_watson,
    watson_start = Start_watson,
    watson_end = End_watson,
    crick_ensembl_id = ensembl_id_crick,
    crick_gene_name = external_gene_name_crick,
    crick_start = Start_crick,
    crick_end = End_crick
  )

# Step 7: Output to CSV
write_csv(bdp_clean, "BDP_Cleaned_Pairs.csv")

# Step 8: Preview result
print(head(bdp_clean))
