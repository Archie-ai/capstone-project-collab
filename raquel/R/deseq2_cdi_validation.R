# BIOT 670I - Microbiome Abundance Consensus Machine
# Real-data differential abundance analysis: DESeq2
# Dataset: cdi_vincent (CDI vs Healthy)

suppressPackageStartupMessages({
  library(DESeq2)
})

# File locations
data_dir <- "data/processed/cdi_vincent"
output_dir <- "results/real/cdi_vincent/deseq2_validation_no_filter"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis settings
alpha <- 0.05

# Load processed dataset
counts <- read.csv(
  file.path(data_dir, "cdi_counts.csv"),
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

# Validate input data
stopifnot(
  !anyDuplicated(rownames(counts)),
  !anyDuplicated(colnames(counts)),
  !anyDuplicated(metadata$SampleID),
  setequal(colnames(counts), metadata$SampleID),
  all(!is.na(counts)),
  all(counts >= 0),
  all(counts == floor(counts))
)

# Align metadata with count columns
metadata <- metadata[
  match(colnames(counts), metadata$SampleID),
]

stopifnot(
  all(metadata$SampleID == colnames(counts)),
  all(metadata$Group %in% c("Case", "Control"))
)

# Set Healthy (Control) as reference
metadata$Group <- factor(
  metadata$Group,
  levels = c("Control", "Case")
)

rownames(metadata) <- metadata$SampleID

# Prepare integer count matrix
count_matrix <- as.matrix(counts)
storage.mode(count_matrix) <- "integer"

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ Group
)

# Run differential abundance analysis
set.seed(42)

dds <- DESeq(
  dds,
  sfType = "poscounts",
  quiet = TRUE
)

# Extract Case vs Control results
# Keep DESeq2's default independent filtering
res <- results(
  dds,
  contrast = c("Group", "Case", "Control"),
  alpha = alpha,
  independentFiltering = FALSE
)

# Convert to data frame
results_df <- data.frame(
  Taxon = rownames(res),
  Effect = res$log2FoldChange,
  PValue = res$pvalue,
  Adjusted_PValue = res$padj,
  stringsAsFactors = FALSE
)

# Identify significant OTUs
results_df$Significant <- (
  !is.na(results_df$Adjusted_PValue) &
  results_df$Adjusted_PValue < alpha
)

# Attach taxonomy annotations
results_df$Taxonomy <- taxonomy$Taxonomy[
  match(results_df$Taxon, taxonomy$OTU_ID)
]

# Sort by adjusted p-value
results_df <- results_df[
  order(results_df$Adjusted_PValue, na.last = TRUE),
]

# Save complete results
write.csv(
  results_df,
  file.path(output_dir, "deseq2_results.csv"),
  row.names = FALSE
)

# Save analysis summary
summary <- data.frame(
  Method = "DESeq2",
  Dataset = "cdi_vincent",
  Samples = nrow(metadata),
  Case_Samples = sum(metadata$Group == "Case"),
  Control_Samples = sum(metadata$Group == "Control"),
  OTUs_Input = nrow(counts),
  OTUs_Tested = sum(!is.na(results_df$Adjusted_PValue)),
  OTUs_With_Missing_Adjusted_PValue = sum(
    is.na(results_df$Adjusted_PValue)
  ),
  Significant_OTUs = sum(results_df$Significant),
  Alpha = alpha
)

write.csv(
  summary,
  file.path(output_dir, "analysis_summary.csv"),
  row.names = FALSE
)

cat("\nDESeq2 analysis completed!\n")
cat("Samples:", nrow(metadata), "\n")
cat("CDI samples:", sum(metadata$Group == "Case"), "\n")
cat("Healthy samples:", sum(metadata$Group == "Control"), "\n")
cat("OTUs input:", nrow(counts), "\n")
cat(
  "OTUs tested:",
  sum(!is.na(results_df$Adjusted_PValue)),
  "\n"
)
cat(
  "OTUs with missing adjusted p-values:",
  sum(is.na(results_df$Adjusted_PValue)),
  "\n"
)
cat(
  "Significant OTUs:",
  sum(results_df$Significant),
  "\n"
)

if (any(!is.na(results_df$Adjusted_PValue))) {
  cat(
    "Minimum adjusted p-value:",
    min(results_df$Adjusted_PValue, na.rm = TRUE),
    "\n"
  )
}

cat("Results saved to:", output_dir, "\n")
