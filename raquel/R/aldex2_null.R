# Load ALDEx2
library(ALDEx2)

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

# Create the Control/Case condition vector
conditions <- metadata$Group

# Run ALDEx2 CLR transformation
# Monte Carlo sampling accounts for uncertainty in count data
clr <- aldex.clr(
  counts,
  conditions,
  mc.samples = 128,
  denom = "all",
  verbose = FALSE
)

# Compare Control and Case using ALDEx2
test <- aldex.ttest(
  clr,
  paired.test = FALSE,
  verbose = FALSE
)

# Calculate effect sizes
effect <- aldex.effect(
  clr,
  verbose = FALSE
)

# Combine ALDEx2 outputs
results <- cbind(test, effect)

# Add taxon names as a column
results$Taxon <- rownames(results)

# Use Wilcoxon Benjamini-Hochberg adjusted p-value
results$Adjusted_PValue <- results$wi.eBH

# Mark significant taxa
results$Significant <- results$Adjusted_PValue < 0.05

# Move Taxon to the first column
results <- results[
  c("Taxon", setdiff(names(results), "Taxon"))
]

# Create output folder if needed
dir.create(
  "results/null",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save results
write.csv(
  results,
  "results/null/aldex2_results.csv",
  row.names = FALSE
)

# Print summary
cat("ALDEx2 analysis complete.\n")
cat("Taxa tested:", nrow(results), "\n")
cat("Significant taxa:", sum(results$Significant), "\n")
cat(
  "Smallest adjusted p-value:",
  round(min(results$Adjusted_PValue, na.rm = TRUE), 4),
  "\n"
)
cat("Results saved to results/null/aldex2_results.csv\n")
