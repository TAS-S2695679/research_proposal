# ---- Libraries ----
library(dplyr)
library(readr)
library(tools)

# ---- Define input and output directories ----
input_dir <- "~/Desktop/bedGraphs_raw"  
output_dir <- "~/Desktop/bedGraphs_scaled"
dir.create(output_dir, showWarnings = FALSE)

# ---- List all .bedGraph files ----
bed_files <- list.files(input_dir, pattern = "\\.bedGraph$", full.names = TRUE)

# ---- Step 1: Compute total signal per file ----
cat("Calculating total signal per file...\n")
total_signals <- sapply(bed_files, function(file) {
  cat("  Reading:", basename(file), "...\n")
  df <- read_tsv(file, col_names = FALSE, progress = FALSE)
  sum(df$X4, na.rm = TRUE)
})
names(total_signals) <- basename(bed_files)

# ---- Step 2: Compute scale factors relative to 0h ----
ref_signal <- total_signals[grepl("0h", names(total_signals))]
cat("\n Reference total signal (0h):", ref_signal, "\n\n")

scale_factors <- ref_signal / total_signals
print(round(scale_factors, 3))

# ---- Step 3: Apply scaling and save new files ----
cat("\n Applying scale factors and saving files...\n")
for (i in seq_along(bed_files)) {
  file <- bed_files[i]
  cat("  Scaling:", basename(file), "\n")
  
  df <- read_tsv(file, col_names = FALSE, progress = FALSE)
  df$X4 <- df$X4 * scale_factors[i]
  
  out_file <- file.path(output_dir, paste0(file_path_sans_ext(basename(file)), "_scaled.bedGraph"))
  write_tsv(df, out_file, col_names = FALSE)
  cat("  Saved:", basename(out_file), "\n\n")
}

cat("All files scaled and saved to:", output_dir, "\n")
