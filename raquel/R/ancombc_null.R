# Load ANCOMBC
library(ANCOMBC)

# Load null count data and metadata
counts <- read.csv(
  "data/null/null_counts.csv",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  "data/null/null_metadata.csv"
)

# Match metadata order to count table samples
metadata <- metadata[
  match(colnames(counts), metadata$SampleID),
]

# Remove taxa with zero counts across all samples
counts <- counts[rowSums(counts) > 0, ]

# Convert counts to samples x taxa for TreeSummarizedExperiment
count_matrix <- as.matrix(counts)

# Create TreeSummarizedExperiment object
tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
  assays = list(counts = count_matrix),
  colData = S4Vectors::DataFrame(metadata)
)

# Run ANCOM-BC2 using Group as the predictor
output <- ancombc2(
  data = tse,
  assay_name = "counts",
  fix_formula = "Group",
  p_adj_method = "BH",
  prv_cut = 0,
  lib_cut = 0,
  group = "Group",
  struc_zero = FALSE,
  neg_lb = FALSE,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  verbose = FALSE
)

# Extract primary results
results <- output$res

# Find the Case vs Control columns
lfc_col <- grep("^lfc_.*Group", names(results), value = TRUE)[1]
q_col <- grep("^q_.*Group", names(results), value = TRUE)[1]

# Create a standardized result table
final_results <- data.frame(
  Taxon = results$taxon,
  Effect = results[[lfc_col]],
  Adjusted_PValue = results[[q_col]]
)

# Mark significant taxa
final_results$Significant <- (
  !is.na(final_results$Adjusted_PValue) &
  final_results$Adjusted_PValue < 0.05
)

# Sort by adjusted p-value
final_results <- final_results[
  order(final_results$Adjusted_PValue),
]

# Create output folder if needed
dir.create(
  "results/null",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save results
write.csv(
  final_results,
  "results/null/ancombc_results.csv",
  row.names = FALSE
)

# Print summary
cat("ANCOM-BC analysis complete.\n")
cat("Taxa tested:", nrow(final_results), "\n")
cat(
  "Significant taxa:",
  sum(final_results$Significant),
  "\n"
)
cat(
  "Smallest adjusted p-value:",
  round(min(final_results$Adjusted_PValue, na.rm = TRUE), 4),
  "\n"
)
cat("Results saved to results/null/ancombc_results.csv\n")
