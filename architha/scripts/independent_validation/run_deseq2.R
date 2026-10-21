# ============================================================
# DESeq2 INDEPENDENT NULL VALIDATION
#
# Microbiome Abundance Consensus Machine
#
# Dataset:
#   100 taxa
#   50 samples
#   25 Control
#   25 Case
#   0 true differentially abundant taxa
#
# Purpose:
# Independently evaluate DESeq2 on the exact same null dataset
# previously analyzed using CLR/Wilcoxon, ALDEx2, and ANCOM-BC2.
#
# Primary DESeq2 analysis:
#   Negative-binomial model
#   Design = ~ Group
#   Reference = Control
#   Contrast = Case vs Control
#   Standard median-ratio size-factor estimation
#   Wald test
#   Benjamini-Hochberg adjustment
#   Default DESeq2 independent filtering retained
#
# Diagnostic analyses:
#   1. Independent filtering disabled
#   2. Independent filtering and Cook's cutoff disabled
# ============================================================


# ============================================================
# 1. CONFIGURATION
# ============================================================

COUNTS_FILE <- file.path(
  "data",
  "processed",
  "simulation_null",
  "null_counts.csv"
)

METADATA_FILE <- file.path(
  "data",
  "processed",
  "simulation_null",
  "null_metadata.csv"
)

TRUTH_FILE <- file.path(
  "data",
  "processed",
  "simulation_null",
  "null_ground_truth.csv"
)

OUTPUT_DIR <- file.path(
  "results",
  "independent_validation",
  "deseq2"
)

ALPHA <- 0.05
RANDOM_SEED <- 42
P_ADJ_METHOD <- "BH"


# ============================================================
# 2. LOAD DESEQ2
# ============================================================

suppressPackageStartupMessages(
  library(DESeq2)
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

separator <- paste(
  rep("=", 70),
  collapse = ""
)

cat("\n", separator, "\n", sep = "")
cat("DESeq2 INDEPENDENT NULL VALIDATION\n")
cat(separator, "\n", sep = "")

cat(
  "\nR version: ",
  R.version.string,
  "\n",
  sep = ""
)

cat(
  "DESeq2 version: ",
  as.character(packageVersion("DESeq2")),
  "\n",
  sep = ""
)


# ============================================================
# 3. LOAD INPUT DATA
# ============================================================

counts <- read.csv(
  COUNTS_FILE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata <- read.csv(
  METADATA_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

truth <- read.csv(
  TRUTH_FILE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\nInput files:\n")
cat("Counts:   ", COUNTS_FILE, "\n", sep = "")
cat("Metadata: ", METADATA_FILE, "\n", sep = "")
cat("Truth:    ", TRUTH_FILE, "\n", sep = "")

cat("\nInput dimensions:\n")

cat(
  "Counts:   ",
  nrow(counts),
  " taxa x ",
  ncol(counts),
  " samples\n",
  sep = ""
)

cat(
  "Metadata: ",
  nrow(metadata),
  " rows x ",
  ncol(metadata),
  " columns\n",
  sep = ""
)

cat(
  "Truth:    ",
  nrow(truth),
  " rows x ",
  ncol(truth),
  " columns\n",
  sep = ""
)


# ============================================================
# 4. VALIDATE METADATA
# ============================================================

required_metadata_columns <- c(
  "SampleID",
  "Group"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(metadata)
)

if (length(missing_metadata_columns) > 0) {

  stop(
    paste(
      "Missing metadata columns:",
      paste(
        missing_metadata_columns,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 5. VALIDATE SAMPLE IDENTIFIERS
# ============================================================

count_samples <- colnames(
  counts
)

metadata_samples <- metadata$SampleID

if (!setequal(
  count_samples,
  metadata_samples
)) {

  stop(
    "Count-table and metadata sample identifiers do not match."
  )
}

metadata <- metadata[
  match(
    count_samples,
    metadata$SampleID
  ),
  ,
  drop = FALSE
]

samples_aligned <- identical(
  colnames(counts),
  metadata$SampleID
)

cat("\nSample alignment:\n")

cat(
  "Counts and metadata aligned: ",
  samples_aligned,
  "\n",
  sep = ""
)

if (!samples_aligned) {

  stop(
    "Counts and metadata are not aligned."
  )
}


# ============================================================
# 6. DEFINE CASE / CONTROL GROUPS
# ============================================================

metadata$Group <- factor(
  metadata$Group,
  levels = c(
    "Control",
    "Case"
  )
)

if (any(is.na(metadata$Group))) {

  stop(
    "Unexpected group labels detected."
  )
}

control_count <- sum(
  metadata$Group == "Control"
)

case_count <- sum(
  metadata$Group == "Case"
)

cat("\nGroup counts:\n")

cat(
  "Control: ",
  control_count,
  "\n",
  sep = ""
)

cat(
  "Case:    ",
  case_count,
  "\n",
  sep = ""
)


# ============================================================
# 7. PREPARE COUNT MATRIX
# ============================================================

counts_matrix <- as.matrix(
  counts
)

storage.mode(
  counts_matrix
) <- "numeric"

if (any(!is.finite(counts_matrix))) {

  stop(
    "Non-finite values detected in the count matrix."
  )
}

negative_count_total <- sum(
  counts_matrix < 0
)

non_integer_total <- sum(
  abs(
    counts_matrix -
    round(counts_matrix)
  ) > 1e-8
)

if (negative_count_total > 0) {

  stop(
    "Negative counts detected."
  )
}

if (non_integer_total > 0) {

  stop(
    "Non-integer counts detected."
  )
}

counts_matrix <- round(
  counts_matrix
)

storage.mode(
  counts_matrix
) <- "integer"


# ============================================================
# 8. RAW DATA QC
# ============================================================

all_zero_taxa <- (
  rowSums(counts_matrix) == 0
)

present_samples <- rowSums(
  counts_matrix > 0
)

present_control <- rowSums(
  counts_matrix[
    ,
    metadata$Group == "Control",
    drop = FALSE
  ] > 0
)

present_case <- rowSums(
  counts_matrix[
    ,
    metadata$Group == "Case",
    drop = FALSE
  ] > 0
)

library_sizes <- colSums(
  counts_matrix
)

zero_cells <- sum(
  counts_matrix == 0
)

zero_fraction <- (
  zero_cells /
  length(counts_matrix)
)

cat("\nRaw-count QC:\n")

cat(
  "Total taxa: ",
  nrow(counts_matrix),
  "\n",
  sep = ""
)

cat(
  "Negative counts: ",
  negative_count_total,
  "\n",
  sep = ""
)

cat(
  "Non-integer counts: ",
  non_integer_total,
  "\n",
  sep = ""
)

cat(
  "All-zero taxa: ",
  sum(all_zero_taxa),
  "\n",
  sep = ""
)

cat(
  "Zero fraction: ",
  sprintf("%.4f", zero_fraction),
  "\n",
  sep = ""
)

cat(
  "Mean library size: ",
  round(mean(library_sizes), 2),
  "\n",
  sep = ""
)

cat(
  "Minimum library size: ",
  min(library_sizes),
  "\n",
  sep = ""
)

cat(
  "Maximum library size: ",
  max(library_sizes),
  "\n",
  sep = ""
)


# ============================================================
# 9. PREPARE DESEQ2 METADATA
# ============================================================

metadata_model <- data.frame(
  Group = metadata$Group,
  row.names = metadata$SampleID,
  stringsAsFactors = FALSE
)

metadata_model$Group <- factor(
  metadata_model$Group,
  levels = c(
    "Control",
    "Case"
  )
)

if (!identical(
  rownames(metadata_model),
  colnames(counts_matrix)
)) {

  stop(
    "DESeq2 metadata row names do not match count columns."
  )
}


# ============================================================
# 10. CREATE DESEQ2 DATASET
# ============================================================

cat("\nDESeq2 configuration:\n")

cat(
  "Design: ~ Group\n"
)

cat(
  "Reference group: Control\n"
)

cat(
  "Contrast: Case vs Control\n"
)

cat(
  "Size-factor method: ratio\n"
)

cat(
  "Statistical test: Wald\n"
)

cat(
  "Primary independent filtering: TRUE\n"
)

cat(
  "P-value adjustment: BH\n"
)

cat(
  "Alpha: ",
  ALPHA,
  "\n",
  sep = ""
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = metadata_model,
  design = ~ Group
)


# ============================================================
# 11. VERIFY DESIGN LEVELS
# ============================================================

cat("\nDESeq2 group levels:\n")

print(
  levels(
    dds$Group
  )
)

if (!identical(
  levels(dds$Group),
  c("Control", "Case")
)) {

  stop(
    "Unexpected DESeq2 group-level ordering."
  )
}


# ============================================================
# 12. RUN DESEQ2
# ============================================================

set.seed(
  RANDOM_SEED
)

cat(
  "\nRunning DESeq2 model fitting...\n"
)

dds <- DESeq(
  dds,
  test = "Wald",
  fitType = "parametric",
  sfType = "ratio",
  quiet = FALSE
)

cat(
  "\nDESeq2 fitting completed.\n"
)


# ============================================================
# 13. PRINT MODEL COEFFICIENTS
# ============================================================

cat(
  "\nDESeq2 result names:\n"
)

print(
  resultsNames(dds)
)


# ============================================================
# 14. EXTRACT PRIMARY DESEQ2 RESULTS
#
# Default independent filtering is retained.
# ============================================================

primary_result <- results(
  dds,
  contrast = c(
    "Group",
    "Case",
    "Control"
  ),
  alpha = ALPHA
)

primary_df <- as.data.frame(
  primary_result
)

primary_df$Taxon <- rownames(
  primary_df
)

primary_df <- primary_df[
  ,
  c(
    "Taxon",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj"
  ),
  drop = FALSE
]


# ============================================================
# 15. EXTRACT RESULT WITHOUT INDEPENDENT FILTERING
#
# This diagnostic tells us whether padj values disappeared
# because of DESeq2 independent filtering.
# ============================================================

no_filter_result <- results(
  dds,
  contrast = c(
    "Group",
    "Case",
    "Control"
  ),
  alpha = ALPHA,
  independentFiltering = FALSE
)

no_filter_df <- as.data.frame(
  no_filter_result
)

no_filter_df$Taxon <- rownames(
  no_filter_df
)

no_filter_df <- no_filter_df[
  ,
  c(
    "Taxon",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj"
  ),
  drop = FALSE
]


# ============================================================
# 16. EXTRACT FULL DIAGNOSTIC RESULT
#
# Independent filtering and Cook's cutoff are both disabled.
# This is diagnostic only and is NOT the primary result.
# ============================================================

full_diagnostic_result <- results(
  dds,
  contrast = c(
    "Group",
    "Case",
    "Control"
  ),
  alpha = ALPHA,
  independentFiltering = FALSE,
  cooksCutoff = FALSE
)

full_diagnostic_df <- as.data.frame(
  full_diagnostic_result
)

full_diagnostic_df$Taxon <- rownames(
  full_diagnostic_df
)

full_diagnostic_df <- full_diagnostic_df[
  ,
  c(
    "Taxon",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj"
  ),
  drop = FALSE
]


# ============================================================
# 17. VERIFY RESULT DIMENSIONS
# ============================================================

cat("\nDESeq2 output dimensions:\n")

cat(
  "Primary result taxa: ",
  nrow(primary_df),
  "\n",
  sep = ""
)

cat(
  "No-filter result taxa: ",
  nrow(no_filter_df),
  "\n",
  sep = ""
)

cat(
  "Full diagnostic result taxa: ",
  nrow(full_diagnostic_df),
  "\n",
  sep = ""
)

if (nrow(primary_df) != nrow(counts_matrix)) {

  stop(
    "Primary DESeq2 result does not contain the full taxon universe."
  )
}


# ============================================================
# 18. ALIGN RESULTS TO ORIGINAL TAXA
# ============================================================

all_taxa <- rownames(
  counts_matrix
)

primary_match <- match(
  all_taxa,
  primary_df$Taxon
)

no_filter_match <- match(
  all_taxa,
  no_filter_df$Taxon
)

full_match <- match(
  all_taxa,
  full_diagnostic_df$Taxon
)

if (any(is.na(primary_match))) {

  stop(
    "Primary DESeq2 results could not be aligned to all taxa."
  )
}

if (any(is.na(no_filter_match))) {

  stop(
    "No-filter DESeq2 results could not be aligned to all taxa."
  )
}

if (any(is.na(full_match))) {

  stop(
    "Full diagnostic DESeq2 results could not be aligned to all taxa."
  )
}


# ============================================================
# 19. BUILD STANDARDIZED RESULT TABLE
# ============================================================

results_standardized <- data.frame(
  Taxon = all_taxa,
  Method = rep(
    "DESeq2",
    length(all_taxa)
  ),
  stringsAsFactors = FALSE
)

results_standardized$Base_Mean <- primary_df$baseMean[
  primary_match
]

results_standardized$Log2_Fold_Change <- primary_df$log2FoldChange[
  primary_match
]

results_standardized$LFC_SE <- primary_df$lfcSE[
  primary_match
]

results_standardized$Wald_Statistic <- primary_df$stat[
  primary_match
]

results_standardized$P_Value <- primary_df$pvalue[
  primary_match
]

results_standardized$Adjusted_P_Value <- primary_df$padj[
  primary_match
]


# ============================================================
# 20. ADD DIAGNOSTIC P-VALUES
# ============================================================

results_standardized$No_Filter_P_Value <- no_filter_df$pvalue[
  no_filter_match
]

results_standardized$No_Filter_Adjusted_P_Value <- no_filter_df$padj[
  no_filter_match
]

results_standardized$Full_Diagnostic_P_Value <- full_diagnostic_df$pvalue[
  full_match
]

results_standardized$Full_Diagnostic_Adjusted_P_Value <- full_diagnostic_df$padj[
  full_match
]


# ============================================================
# 21. DEFINE TEST AND FILTER STATUS
# ============================================================

results_standardized$Tested <- !is.na(
  results_standardized$P_Value
)

results_standardized$Adjusted_P_Available <- !is.na(
  results_standardized$Adjusted_P_Value
)

results_standardized$Independent_Filtered <- (
  results_standardized$Tested &
  is.na(results_standardized$Adjusted_P_Value) &
  !is.na(results_standardized$No_Filter_Adjusted_P_Value)
)

results_standardized$Primary_Significant <- (
  !is.na(results_standardized$Adjusted_P_Value) &
  results_standardized$Adjusted_P_Value < ALPHA
)

results_standardized$No_Filter_Significant <- (
  !is.na(results_standardized$No_Filter_Adjusted_P_Value) &
  results_standardized$No_Filter_Adjusted_P_Value < ALPHA
)

results_standardized$Full_Diagnostic_Significant <- (
  !is.na(results_standardized$Full_Diagnostic_Adjusted_P_Value) &
  results_standardized$Full_Diagnostic_Adjusted_P_Value < ALPHA
)


# ============================================================
# 22. STANDARDIZE EFFECT DIRECTION
# ============================================================

results_standardized$Effect_Direction <- ifelse(
  is.na(results_standardized$Log2_Fold_Change),
  NA,
  ifelse(
    results_standardized$Log2_Fold_Change > 0,
    "Case_higher",
    ifelse(
      results_standardized$Log2_Fold_Change < 0,
      "Control_higher",
      "No_difference"
    )
  )
)


# ============================================================
# 23. ADD PREVALENCE INFORMATION
# ============================================================

taxon_match <- match(
  results_standardized$Taxon,
  names(present_samples)
)

results_standardized$Present_Samples <- as.integer(
  present_samples[
    taxon_match
  ]
)

results_standardized$Present_Control <- as.integer(
  present_control[
    taxon_match
  ]
)

results_standardized$Present_Case <- as.integer(
  present_case[
    taxon_match
  ]
)

results_standardized$All_Zero <- as.logical(
  all_zero_taxa[
    taxon_match
  ]
)


# ============================================================
# 24. MERGE GROUND TRUTH
# ============================================================

truth_match <- match(
  results_standardized$Taxon,
  rownames(truth)
)

if (any(is.na(truth_match))) {

  stop(
    "Ground truth could not be aligned."
  )
}

truth_aligned <- truth[
  truth_match,
  ,
  drop = FALSE
]

for (column_name in colnames(truth_aligned)) {

  results_standardized[
    ,
    column_name
  ] <- truth_aligned[
    ,
    column_name,
    drop = TRUE
  ]
}


# ============================================================
# 25. ADD ANALYSIS METADATA
# ============================================================

results_standardized$Alpha <- ALPHA
results_standardized$P_Adjust_Method <- P_ADJ_METHOD
results_standardized$Design <- "~ Group"
results_standardized$Reference_Group <- "Control"
results_standardized$Contrast <- "Case_vs_Control"
results_standardized$Size_Factor_Method <- "ratio"
results_standardized$Test_Type <- "Wald"
results_standardized$Scenario <- "Null_Baseline_Low_Sparsity"
results_standardized$Seed <- RANDOM_SEED

results_standardized$DESeq2_Version <- as.character(
  packageVersion("DESeq2")
)


# ============================================================
# 26. EXTRACT SIZE FACTORS
# ============================================================

sample_size_factors <- sizeFactors(
  dds
)

size_factor_table <- data.frame(
  SampleID = names(sample_size_factors),
  Group = as.character(
    metadata_model[
      names(sample_size_factors),
      "Group"
    ]
  ),
  Raw_Library_Size = as.numeric(
    library_sizes[
      names(sample_size_factors)
    ]
  ),
  Size_Factor = as.numeric(
    sample_size_factors
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 27. EXTRACT NORMALIZED COUNTS
# ============================================================

normalized_counts <- counts(
  dds,
  normalized = TRUE
)


# ============================================================
# 28. SAFE MINIMUM FUNCTION
# ============================================================

safe_min <- function(x) {

  usable <- x[
    is.finite(x)
  ]

  if (length(usable) == 0) {

    return(
      NA_real_
    )
  }

  return(
    min(usable)
  )
}


# ============================================================
# 29. SUMMARY STATISTICS
# ============================================================

n_available <- nrow(
  results_standardized
)

n_tested <- sum(
  results_standardized$Tested
)

n_primary_padj_available <- sum(
  results_standardized$Adjusted_P_Available
)

n_primary_padj_missing <- sum(
  !results_standardized$Adjusted_P_Available
)

n_independent_filtered <- sum(
  results_standardized$Independent_Filtered,
  na.rm = TRUE
)

n_primary_significant <- sum(
  results_standardized$Primary_Significant,
  na.rm = TRUE
)

n_no_filter_significant <- sum(
  results_standardized$No_Filter_Significant,
  na.rm = TRUE
)

n_full_diagnostic_significant <- sum(
  results_standardized$Full_Diagnostic_Significant,
  na.rm = TRUE
)

minimum_raw_p <- safe_min(
  results_standardized$P_Value
)

minimum_adjusted_p <- safe_min(
  results_standardized$Adjusted_P_Value
)

minimum_no_filter_adjusted_p <- safe_min(
  results_standardized$No_Filter_Adjusted_P_Value
)

minimum_full_diagnostic_adjusted_p <- safe_min(
  results_standardized$Full_Diagnostic_Adjusted_P_Value
)


# ============================================================
# 30. SORT STANDARDIZED RESULTS
# ============================================================

results_sorted <- results_standardized[
  order(
    results_standardized$Adjusted_P_Value,
    results_standardized$P_Value,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]


# ============================================================
# 31. SAVE OUTPUT FILES
# ============================================================

primary_file <- file.path(
  OUTPUT_DIR,
  "deseq2_primary_results.csv"
)

no_filter_file <- file.path(
  OUTPUT_DIR,
  "deseq2_no_independent_filter_results.csv"
)

full_diagnostic_file <- file.path(
  OUTPUT_DIR,
  "deseq2_full_diagnostic_results.csv"
)

standardized_file <- file.path(
  OUTPUT_DIR,
  "deseq2_standardized_results.csv"
)

size_factor_file <- file.path(
  OUTPUT_DIR,
  "deseq2_size_factors.csv"
)

normalized_counts_file <- file.path(
  OUTPUT_DIR,
  "deseq2_normalized_counts.csv"
)

write.csv(
  primary_df,
  primary_file,
  row.names = FALSE
)

write.csv(
  no_filter_df,
  no_filter_file,
  row.names = FALSE
)

write.csv(
  full_diagnostic_df,
  full_diagnostic_file,
  row.names = FALSE
)

write.csv(
  results_sorted,
  standardized_file,
  row.names = FALSE
)

write.csv(
  size_factor_table,
  size_factor_file,
  row.names = FALSE
)

write.csv(
  normalized_counts,
  normalized_counts_file,
  row.names = TRUE
)


# ============================================================
# 32. CREATE TEXT SUMMARY
# ============================================================

summary_file <- file.path(
  OUTPUT_DIR,
  "deseq2_summary.txt"
)

summary_lines <- c(
  "DESeq2 INDEPENDENT NULL VALIDATION",
  paste(rep("=", 60), collapse = ""),
  "",
  paste0(
    "DESeq2 version: ",
    packageVersion("DESeq2")
  ),
  "Design: ~ Group",
  "Reference group: Control",
  "Contrast: Case vs Control",
  "Size-factor method: ratio",
  "Test: Wald",
  "Primary independent filtering: TRUE",
  "P-value adjustment: BH",
  paste0(
    "Alpha: ",
    ALPHA
  ),
  "",
  paste0(
    "Taxa available: ",
    n_available
  ),
  paste0(
    "Taxa with raw p-values: ",
    n_tested
  ),
  paste0(
    "Taxa with primary adjusted p-values: ",
    n_primary_padj_available
  ),
  paste0(
    "Taxa with missing primary adjusted p-values: ",
    n_primary_padj_missing
  ),
  paste0(
    "Taxa identified as independently filtered: ",
    n_independent_filtered
  ),
  paste0(
    "Primary significant taxa: ",
    n_primary_significant
  ),
  paste0(
    "Significant without independent filtering: ",
    n_no_filter_significant
  ),
  paste0(
    "Significant with filtering and Cook's cutoff disabled: ",
    n_full_diagnostic_significant
  ),
  paste0(
    "Minimum primary raw p-value: ",
    signif(
      minimum_raw_p,
      10
    )
  ),
  paste0(
    "Minimum primary adjusted p-value: ",
    signif(
      minimum_adjusted_p,
      10
    )
  ),
  paste0(
    "Minimum no-filter adjusted p-value: ",
    signif(
      minimum_no_filter_adjusted_p,
      10
    )
  ),
  paste0(
    "Minimum full-diagnostic adjusted p-value: ",
    signif(
      minimum_full_diagnostic_adjusted_p,
      10
    )
  ),
  "",
  paste0(
    "Minimum size factor: ",
    signif(
      min(sample_size_factors),
      6
    )
  ),
  paste0(
    "Median size factor: ",
    signif(
      median(sample_size_factors),
      6
    )
  ),
  paste0(
    "Maximum size factor: ",
    signif(
      max(sample_size_factors),
      6
    )
  ),
  "",
  "Known true DA taxa: 0",
  "",
  paste0(
    "Under this complete-null simulation, ",
    "every significant taxon represents a ",
    "false-positive discovery."
  )
)

writeLines(
  summary_lines,
  summary_file
)


# ============================================================
# 33. PRINT MAIN SUMMARY
# ============================================================

cat("\n", separator, "\n", sep = "")
cat("RESULT SUMMARY\n")
cat(separator, "\n", sep = "")

cat(
  "Taxa available: ",
  n_available,
  "\n",
  sep = ""
)

cat(
  "Taxa with raw p-values: ",
  n_tested,
  "\n",
  sep = ""
)

cat(
  "Taxa with primary adjusted p-values: ",
  n_primary_padj_available,
  "\n",
  sep = ""
)

cat(
  "Taxa with missing primary adjusted p-values: ",
  n_primary_padj_missing,
  "\n",
  sep = ""
)

cat(
  "Taxa independently filtered: ",
  n_independent_filtered,
  "\n",
  sep = ""
)

cat(
  "Primary significant taxa at BH q < ",
  ALPHA,
  ": ",
  n_primary_significant,
  "\n",
  sep = ""
)

cat(
  "Significant without independent filtering: ",
  n_no_filter_significant,
  "\n",
  sep = ""
)

cat(
  "Significant with filtering and Cook's cutoff disabled: ",
  n_full_diagnostic_significant,
  "\n",
  sep = ""
)

cat(
  "Minimum primary raw p-value: ",
  signif(
    minimum_raw_p,
    10
  ),
  "\n",
  sep = ""
)

cat(
  "Minimum primary adjusted p-value: ",
  signif(
    minimum_adjusted_p,
    10
  ),
  "\n",
  sep = ""
)


# ============================================================
# 34. PRINT SIZE-FACTOR SUMMARY
# ============================================================

cat("\nSize-factor summary:\n")

cat(
  "Minimum: ",
  signif(
    min(sample_size_factors),
    6
  ),
  "\n",
  sep = ""
)

cat(
  "Median:  ",
  signif(
    median(sample_size_factors),
    6
  ),
  "\n",
  sep = ""
)

cat(
  "Maximum: ",
  signif(
    max(sample_size_factors),
    6
  ),
  "\n",
  sep = ""
)


# ============================================================
# 35. INSPECT TAXON_033
# ============================================================

display_columns <- c(
  "Taxon",
  "Present_Samples",
  "Present_Control",
  "Present_Case",
  "Base_Mean",
  "Log2_Fold_Change",
  "P_Value",
  "Adjusted_P_Value",
  "Adjusted_P_Available",
  "Independent_Filtered",
  "Primary_Significant",
  "No_Filter_Adjusted_P_Value",
  "Full_Diagnostic_Adjusted_P_Value",
  "Effect_Direction"
)

cat(
  "\nTaxon_033 comparison:\n"
)

print(
  results_standardized[
    results_standardized$Taxon == "Taxon_033",
    display_columns,
    drop = FALSE
  ],
  row.names = FALSE
)


# ============================================================
# 36. PRINT TOP 15 PRIMARY RESULTS
# ============================================================

cat(
  "\nTop 15 taxa by primary DESeq2 adjusted p-value:\n"
)

print(
  head(
    results_sorted[
      ,
      display_columns,
      drop = FALSE
    ],
    15
  ),
  row.names = FALSE
)


# ============================================================
# 37. PRINT RAREST TAXA
# ============================================================

rarest_results <- results_standardized[
  order(
    results_standardized$Present_Samples,
    results_standardized$Taxon
  ),
  ,
  drop = FALSE
]

cat(
  "\nRarest taxa and DESeq2 results:\n"
)

print(
  head(
    rarest_results[
      ,
      display_columns,
      drop = FALSE
    ],
    10
  ),
  row.names = FALSE
)


# ============================================================
# 38. PRINT OUTPUT FILES
# ============================================================

cat("\nOutput files:\n")
cat(primary_file, "\n")
cat(no_filter_file, "\n")
cat(full_diagnostic_file, "\n")
cat(standardized_file, "\n")
cat(size_factor_file, "\n")
cat(normalized_counts_file, "\n")
cat(summary_file, "\n")

cat(
  "\nDESeq2 independent null validation completed.\n"
)
