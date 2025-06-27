# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("rtracklayer", "GenomicRanges"))

library(rtracklayer)
library(GenomicRanges)

# Set directories
input_dir <- "ttseq"
output_dir <- "ttseq_normalised_bw"
dir.create(output_dir, showWarnings = FALSE)

# Define timepoints and expected filename patterns
timepoints <- c("0h", "3h", "6h", "9h", "12h", "15h")

# Helper function to find the 4 relevant files per timepoint
get_bw_files <- function(tp) {
  list.files(input_dir, pattern = paste0(tp, ".*(plus|minus).*\\.bw$"), full.names = TRUE)
}

# Step 1: Calculate total signal per timepoint with progress
total_signal <- setNames(numeric(length(timepoints)), timepoints)

for (tp in timepoints) {
  cat("\n⏳ Processing timepoint:", tp, "\n")
  files <- get_bw_files(tp)
  
  tp_sum <- 0
  for (f in files) {
    cat("  📂 Reading:", basename(f), "... ")
    bw <- import(f, as = "RleList")
    s <- sum(as.numeric(unlist(bw)))
    cat("Done. Signal =", round(s), "\n")
    tp_sum <- tp_sum + s
  }
  
  total_signal[tp] <- tp_sum
  cat("✅ Total signal for", tp, "=", round(tp_sum), "\n")
}

cat("\n🚀 All timepoints processed. Scale factors next...\n")

# Step 2: Compute scale factors using 0h as reference
ref_total <- total_signal["0h"]
scale_factors <- ref_total / total_signal
print(round(scale_factors, 3))

# Step 3: Merge, normalise, and export .bw file per timepoint
for (tp in timepoints) {
  cat("\nProcessing:", tp)
  files <- get_bw_files(tp)
  sf <- scale_factors[tp]
  
  # Merge signal from 4 .bw files
  merged_bw <- Reduce(`+`, lapply(files, function(f) import(f, as = "RleList")))
  
  # Apply scaling
  norm_bw <- merged_bw * sf
  
  # Export
  out_file <- file.path(output_dir, paste0("GSE174774_Oct4_", tp, "_norm.bw"))
  export(norm_bw, out_file)
  cat(" → saved:", basename(out_file), "\n")
}

cat("\n All normalised files saved to:", output_dir, "\n")
