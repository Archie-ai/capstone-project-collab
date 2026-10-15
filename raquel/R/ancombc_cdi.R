# BIOT 670I - Microbiome Abundance Consensus Machine
# Real-data differential abundance analysis: ANCOM-BC2
# Dataset: cdi_vincent (CDI vs Healthy)

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(TreeSummarizedExperiment)
  library(SummarizedExperiment)
})

# File locations
data_dir <- "data/processed/cdi_vincent"
output_dir <- "results/real/cdi_vincent/ancombc"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

# Validate and align samples
stopifnot(
  !anyDuplicated(rownames(counts)),
  !anyDuplicated(colnames(counts)),
  !anyDuplicated(metadata$SampleID),
  setequal(colnames(counts), metadata$SampleID),
  all(!is.na(counts)),
  all(counts >= 0),
  all(counts == floor(counts))
)

metadata <- metadata[
  match(colnames(counts), metadata$SampleID),
]

stopifnot(all(metadata$SampleID == colnames(counts)))

# Control is the reference group
metadata$Group <- factor(
  metadata$Group,
  levels = c("Control", "Case")
)

# Construct TreeSummarizedExperiment
count_matrix <- as.matrix(counts)
storage.mode(count_matrix) <- "integer"

sample_data <- S4Vectors::DataFrame(
  Group = metadata$Group,
  row.names = metadata$SampleID
)

tse <- TreeSummarizedExperiment(
  assays = list(counts = count_matrix),
  colData = sample_data
)

# Run ANCOM-BC2
set.seed(42)

fit <- ancombc2(
  data = tse,
  assay_name = "counts",
  tax_level = NULL,
  fix_formula = "Group",
  rand_formula = NULL,
  p_adj_method = "BH",
  prv_cut = 0,
  lib_cut = 0,
  s0_perc = 0.05,
  group = NULL,
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = alpha,
  n_cl = 1,
  verbose = FALSE
)

# Extract differential abundance results
res <- fit$res

required <- c(
  "taxon",
  "lfc_GroupCase",
  "q_GroupCase",
  "diff_GroupCase",
  "passed_ss_GroupCase"
)

if (!all(required %in% colnames(res))) {
  stop(
    "Expected ANCOM-BC2 columns not found: ",
    paste(setdiff(required, colnames(res)), collapse = ", ")
  )
}

results <- data.frame(
  Taxon = res$taxon,
  Effect = res$lfc_GroupCase,
  Adjusted_PValue = res$q_GroupCase,
  Passed_Sensitivity = res$passed_ss_GroupCase,
  stringsAsFactors = FALSE
)

# Standard adjusted-p-value significance
results$Significant <- (
  !is.na(results$Adjusted_PValue) &
  results$Adjusted_PValue < alpha
)

# Robust significance requires passing sensitivity analysis
results$Robust_Significant <- (
  results$Significant &
  !is.na(results$Passed_Sensitivity) &
  results$Passed_Sensitivity
)

# Add taxonomy
results$Taxonomy <- taxonomy$Taxonomy[
  match(results$Taxon, taxonomy$OTU_ID)
]

# Check that each taxon appears only once
stopifnot(!anyDuplicated(results$Taxon))

# Sort by adjusted p-value
results <- results[
  order(results$Adjusted_PValue, na.last = TRUE),
]

# Save complete results
write.csv(
  results,
  file.path(output_dir, "ancombc_results.csv"),
  row.names = FALSE
)

# Save summary
summary <- data.frame(
  Method = "ANCOMBC2",
  Dataset = "cdi_vincent",
  Samples = nrow(metadata),
  Case_Samples = sum(metadata$Group == "Case"),
  Control_Samples = sum(metadata$Group == "Control"),
  OTUs_Input = nrow(counts),
  OTUs_Reported = nrow(results),
  OTUs_Tested = sum(!is.na(results$Adjusted_PValue)),
  Significant_OTUs = sum(results$Significant),
  Robust_Significant_OTUs = sum(results$Robust_Significant),
  Alpha = alpha
)

write.csv(
  summary,
  file.path(output_dir, "analysis_summary.csv"),
  row.names = FALSE
)

cat("\nANCOM-BC2 analysis completed!\n")
cat("Samples:", nrow(metadata), "\n")
cat("CDI samples:", sum(metadata$Group == "Case"), "\n")
cat("Healthy samples:", sum(metadata$Group == "Control"), "\n")
cat("OTUs input:", nrow(counts), "\n")
cat("OTUs reported:", nrow(results), "\n")
cat(
  "OTUs tested:",
  sum(!is.na(results$Adjusted_PValue)),
  "\n"
)
cat(
  "Significant OTUs (adjusted p-value):",
  sum(results$Significant),
  "\n"
)
cat(
  "Robust significant OTUs:",
  sum(results$Robust_Significant),
  "\n"
)
cat("Results saved to:", output_dir, "\n")
