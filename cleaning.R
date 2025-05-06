# Load required libraries
library(dplyr)
library(readr)
library(stringr)

# Step 1: Load the Catalogue
catalogue <- read_csv("Catalogue.csv")  # Adjust path if needed

# Step 2: Filter to transcript-level entries only
transcripts <- catalogue %>%
  filter(Feature == "transcript")

# Step 3: Extract base BDP ID by removing '_2' or '_3' suffix from Annotation
transcripts <- transcripts %>%
  mutate(BDP_ID = str_replace(Annotation, "_[23]$", ""))

# Step 4: Separate Watson (+ strand) and Crick (− strand) entries
watson_unique <- watson %>%
  group_by(BDP_ID) %>%
  slice(1) %>%
  ungroup()

crick_unique <- crick %>%
  group_by(BDP_ID) %>%
  slice(1) %>%
  ungroup()

# Now safely join without the many-to-many warning
bdp_pairs <- inner_join(watson_unique, crick_unique, by = "BDP_ID")

# Step 6: Output result to CSV
write_csv(bdp_pairs, "BDP_Cleaned_Pairs.csv")

# Step 7: Preview result
print(head(bdp_pairs))
