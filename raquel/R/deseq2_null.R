# Load DESeq2
library(DESeq2)

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

# Use sample IDs as metadata row names
rownames(metadata) <- metadata$SampleID

# Set Control as the reference group
metadata$Group <- factor(
  metadata$Group,
  levels = c("Control", "Case")
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = metadata,
  design = ~ Group
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Compare Case against Control
results <- results(
  dds,
  contrast = c("Group", "Case", "Control"),
  alpha = 0.05
)

# Convert results to a regular data frame
results <- as.data.frame(results)
results$Taxon <- rownames(results)

# Create standardized columns
results$Adjusted_PValue <- results$padj

results$Significant <- (
  !is.na(results$Adjusted_PValue) &
  results$Adjusted_PValue < 0.05
)

# Put Taxon first and sort by adjusted p-value
results <- results[
  c("Taxon", setdiff(names(results), "Taxon"))
]

results <- results[
  order(results$Adjusted_PValue, na.last = TRUE),
]

# Save results
dir.create(
  "results/null",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  results,
  "results/null/deseq2_results.csv",
  row.names = FALSE
)

# Count taxa with a usable adjusted p-value
taxa_tested <- sum(!is.na(results$Adjusted_PValue))

cat("DESeq2 analysis complete.\n")
cat("Taxa in input:", nrow(results), "\n")
cat("Taxa tested:", taxa_tested, "\n")
cat("Significant taxa:", sum(results$Significant), "\n")
cat(
  "Smallest adjusted p-value:",
  round(min(results$Adjusted_PValue, na.rm = TRUE), 4),
  "\n"
)
cat("Results saved to results/null/deseq2_results.csv\n")
