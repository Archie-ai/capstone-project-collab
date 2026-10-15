# BIOT 670I - Microbiome Abundance Consensus Machine
# Real-data differential abundance analysis: ALDEx2
# Dataset: cdi_vincent (CDI vs Healthy)

suppressPackageStartupMessages({
  library(ALDEx2)
})

# File locations
data_dir <- "data/processed/cdi_vincent"
output_dir <- "results/real/cdi_vincent/aldex2"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis settings
mc_samples <- 128
alpha <- 0.05

# Load processed dataset
counts <- read.csv(
  file.path(data_dir, "cdi_counts_filtered.csv"),
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  file.path(data_dir, "cdi_metadata.csv"),
  stringsAsFactors = FALSE
)

taxonomy <- read.csv(
  file.path(data_dir, "cdi_taxonomy_mapping.csv"),
  stringsAsFactors = FALSE
)

# Validate sample IDs
stopifnot(
  !anyDuplicated(rownames(counts)),
  !anyDuplicated(colnames(counts)),
  !anyDuplicated(metadata$SampleID),
  setequal(colnames(counts), metadata$SampleID)
)

# Align metadata to count columns
metadata <- metadata[
  match(colnames(counts), metadata$SampleID),
]

stopifnot(
  all(metadata$SampleID == colnames(counts)),
  all(metadata$Group %in% c("Case", "Control")),
  all(!is.na(counts)),
  all(counts >= 0),
  all(counts == floor(counts))
)

# Ensure Control is the reference group
groups <- factor(
  metadata$Group,
  levels = c("Control", "Case")
)

# Run ALDEx2
set.seed(42)

clr <- aldex.clr(
  counts,
  conds = as.character(groups),
  mc.samples = mc_samples,
  denom = "all",
  verbose = FALSE
)

test_results <- aldex.ttest(
  clr,
  paired.test = FALSE
)

effect_results <- aldex.effect(
  clr,
  verbose = FALSE
)

# Combine statistical tests and effect estimates
results <- data.frame(
  Taxon = rownames(test_results),
  test_results,
  effect_results,
  check.names = FALSE
)

# Match our simulated-data significance definition
results$Adjusted_PValue <- results$wi.eBH

results$Significant <- (
  !is.na(results$Adjusted_PValue) &
  results$Adjusted_PValue < alpha
)

# Add taxonomy annotations
results$Taxonomy <- taxonomy$Taxonomy[
  match(results$Taxon, taxonomy$OTU_ID)
]

# Sort by adjusted p-value
results <- results[
  order(results$Adjusted_PValue, na.last = TRUE),
]

# Save complete results
write.csv(
  results,
  file.path(output_dir, "aldex2_results.csv"),
  row.names = FALSE
)

# Save analysis summary
summary <- data.frame(
  Method = "ALDEx2",
  Dataset = "cdi_vincent",
  Samples = nrow(metadata),
  Case_Samples = sum(groups == "Case"),
  Control_Samples = sum(groups == "Control"),
  OTUs_Input = nrow(counts),
  OTUs_Tested = sum(!is.na(results$Adjusted_PValue)),
  Significant_OTUs = sum(results$Significant),
  Monte_Carlo_Samples = mc_samples,
  Alpha = alpha
)

write.csv(
  summary,
  file.path(output_dir, "analysis_summary.csv"),
  row.names = FALSE
)

cat("\nALDEx2 analysis completed!\n")
cat("Samples:", nrow(metadata), "\n")
cat("CDI samples:", sum(groups == "Case"), "\n")
cat("Healthy samples:", sum(groups == "Control"), "\n")
cat("OTUs input:", nrow(counts), "\n")
cat(
  "OTUs tested:",
  sum(!is.na(results$Adjusted_PValue)),
  "\n"
)
cat("Significant OTUs:", sum(results$Significant), "\n")
cat(
  "Minimum adjusted p-value:",
  min(results$Adjusted_PValue, na.rm = TRUE),
  "\n"
)
cat("Results saved to:", output_dir, "\n")
