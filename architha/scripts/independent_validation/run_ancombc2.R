# ============================================================
# ANCOM-BC2 INDEPENDENT NULL VALIDATION
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
# ANCOMBC:
#   Version 2.15.2 bugfix branch
#
# Important workaround:
#   ANCOM-BC2 2.15.x has a known bug when metadata contains
#   exactly one column. Therefore a harmless Dummy column is
#   added to metadata but is NOT included in the model.
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
  "ancombc2"
)

ALPHA <- 0.05
RANDOM_SEED <- 42
P_ADJ_METHOD <- "BH"
PREVALENCE_CUT <- 0
LIBRARY_CUT <- 0
S0_PERC <- 0.05


# ============================================================
# 2. LOAD PACKAGE
# ============================================================

suppressPackageStartupMessages(
  library(ANCOMBC)
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

separator <- paste(rep("=", 70), collapse = "")

cat("\n", separator, "\n", sep = "")
cat("ANCOM-BC2 INDEPENDENT NULL VALIDATION\n")
cat(separator, "\n", sep = "")

cat(
  "\nR version: ",
  R.version.string,
  "\n",
  sep = ""
)

cat(
  "ANCOMBC version: ",
  as.character(packageVersion("ANCOMBC")),
  "\n",
  sep = ""
)


# ============================================================
# 3. LOAD DATA
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
# 5. ALIGN SAMPLES
# ============================================================

count_samples <- colnames(counts)

if (!setequal(count_samples, metadata$SampleID)) {

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
# 6. DEFINE GROUPS
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
    "Unexpected Case/Control group label."
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
    "Count matrix contains non-finite values."
  )
}

negative_count_total <- sum(
  counts_matrix < 0
)

non_integer_total <- sum(
  abs(
    counts_matrix - round(counts_matrix)
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
# 8. RAW COUNT QC
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

zero_cells <- sum(
  counts_matrix == 0
)

zero_fraction <- (
  zero_cells / length(counts_matrix)
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


# ============================================================
# 9. BUILD ANCOM-BC2 METADATA
#
# IMPORTANT:
# Adding Dummy avoids the known ANCOMBC single-column
# metadata bug. Dummy is NOT part of fix_formula.
# ============================================================

metadata_model <- data.frame(
  Group = metadata$Group,
  Dummy = rep(
    1L,
    nrow(metadata)
  ),
  stringsAsFactors = FALSE
)

metadata_model$Group <- factor(
  metadata_model$Group,
  levels = c(
    "Control",
    "Case"
  )
)

rownames(
  metadata_model
) <- metadata$SampleID

if (!identical(
  rownames(metadata_model),
  colnames(counts_matrix)
)) {

  stop(
    "Metadata row names do not match count-table columns."
  )
}

cat("\nANCOM-BC2 metadata columns:\n")

print(
  colnames(metadata_model)
)

cat(
  "\nMetadata preview:\n"
)

print(
  head(metadata_model)
)


# ============================================================
# 10. ANALYSIS CONFIGURATION
# ============================================================

cat("\nANCOM-BC2 configuration:\n")

cat("Fixed formula: Group\n")
cat("Reference group: Control\n")
cat("Target coefficient: GroupCase\n")

cat(
  "P-value adjustment: ",
  P_ADJ_METHOD,
  "\n",
  sep = ""
)

cat(
  "Prevalence cutoff: ",
  PREVALENCE_CUT,
  "\n",
  sep = ""
)

cat(
  "Library-size cutoff: ",
  LIBRARY_CUT,
  "\n",
  sep = ""
)

cat(
  "Pseudo-count sensitivity analysis: TRUE\n"
)

cat(
  "Conservative sensitivity mode: TRUE\n"
)

cat(
  "Alpha: ",
  ALPHA,
  "\n",
  sep = ""
)


# ============================================================
# 11. RUN ANCOM-BC2
# ============================================================

set.seed(
  RANDOM_SEED
)

cat(
  "\nRunning ANCOM-BC2...\n"
)

fit <- ANCOMBC::ancombc2(
  data = counts_matrix,
  taxa_are_rows = TRUE,
  meta_data = metadata_model,
  fix_formula = "Group",
  rand_formula = NULL,
  p_adj_method = P_ADJ_METHOD,
  pseudo_sens = TRUE,
  prv_cut = PREVALENCE_CUT,
  lib_cut = LIBRARY_CUT,
  s0_perc = S0_PERC,
  group = NULL,
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = ALPHA,
  n_cl = 1,
  verbose = TRUE,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  conservative = TRUE
)

cat(
  "\nANCOM-BC2 fitting completed.\n"
)


# ============================================================
# 12. CHECK PRIMARY RESULT
# ============================================================

if (is.null(fit$res)) {

  stop(
    "ANCOM-BC2 returned no primary result table."
  )
}

primary_result <- as.data.frame(
  fit$res,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\nPrimary-result dimensions:\n")

cat(
  nrow(primary_result),
  " rows x ",
  ncol(primary_result),
  " columns\n",
  sep = ""
)

cat(
  "\nPrimary-result columns:\n"
)

print(
  colnames(primary_result)
)


# ============================================================
# 13. SAVE RAW RESULT
# ============================================================

raw_result_file <- file.path(
  OUTPUT_DIR,
  "ancombc2_raw_primary_result.csv"
)

write.csv(
  primary_result,
  raw_result_file,
  row.names = FALSE
)


# ============================================================
# 14. CHECK EXPECTED CASE-vs-CONTROL COLUMNS
# ============================================================

expected_columns <- c(
  "taxon",
  "lfc_GroupCase",
  "se_GroupCase",
  "W_GroupCase",
  "p_GroupCase",
  "q_GroupCase",
  "diff_GroupCase",
  "passed_ss_GroupCase",
  "diff_robust_GroupCase"
)

missing_columns <- setdiff(
  expected_columns,
  colnames(primary_result)
)

if (length(missing_columns) > 0) {

  cat(
    "\nExpected columns that were not found:\n"
  )

  print(
    missing_columns
  )

  cat(
    "\nNo statistical interpretation will be forced.\n"
  )

  stop(
    "Inspect the printed ANCOM-BC2 column names before continuing."
  )
}


# ============================================================
# 15. ALIGN RESULTS TO ORIGINAL 100-TAXON UNIVERSE
# ============================================================

all_taxa <- rownames(
  counts_matrix
)

result_taxa <- as.character(
  primary_result[
    ,
    "taxon",
    drop = TRUE
  ]
)

result_match <- match(
  all_taxa,
  result_taxa
)

results <- data.frame(
  Taxon = all_taxa,
  Method = rep(
    "ANCOM_BC2",
    length(all_taxa)
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 16. HELPER FUNCTIONS
# ============================================================

extract_numeric <- function(column_name) {

  output <- rep(
    NA_real_,
    length(all_taxa)
  )

  matched <- !is.na(
    result_match
  )

  values <- as.numeric(
    primary_result[
      ,
      column_name,
      drop = TRUE
    ]
  )

  output[matched] <- values[
    result_match[matched]
  ]

  return(output)
}


extract_logical <- function(column_name) {

  output <- rep(
    NA,
    length(all_taxa)
  )

  matched <- !is.na(
    result_match
  )

  values <- as.logical(
    primary_result[
      ,
      column_name,
      drop = TRUE
    ]
  )

  output[matched] <- values[
    result_match[matched]
  ]

  return(output)
}


# ============================================================
# 17. EXTRACT ANCOM-BC2 STATISTICS
# ============================================================

results$Log_Fold_Change <- extract_numeric(
  "lfc_GroupCase"
)

results$Standard_Error <- extract_numeric(
  "se_GroupCase"
)

results$W_Statistic <- extract_numeric(
  "W_GroupCase"
)

results$P_Value <- extract_numeric(
  "p_GroupCase"
)

results$Adjusted_P_Value <- extract_numeric(
  "q_GroupCase"
)

results$Package_Diff_Call <- extract_logical(
  "diff_GroupCase"
)

results$Passed_Sensitivity <- extract_logical(
  "passed_ss_GroupCase"
)

results$Robust_Significant <- extract_logical(
  "diff_robust_GroupCase"
)


# ============================================================
# 18. DEFINE TEST STATUS
# ============================================================

results$Tested <- !is.na(
  results$Adjusted_P_Value
)

results$Standard_Significant <- (
  results$Tested &
  results$Adjusted_P_Value < ALPHA
)


# ============================================================
# 19. STANDARDIZE CASE-vs-CONTROL DIRECTION
# ============================================================

results$Effect_Direction <- ifelse(
  is.na(results$Log_Fold_Change),
  NA,
  ifelse(
    results$Log_Fold_Change > 0,
    "Case_higher",
    ifelse(
      results$Log_Fold_Change < 0,
      "Control_higher",
      "No_difference"
    )
  )
)


# ============================================================
# 20. ADD PREVALENCE
# ============================================================

taxon_match <- match(
  results$Taxon,
  names(present_samples)
)

results$Present_Samples <- as.integer(
  present_samples[
    taxon_match
  ]
)

results$Present_Control <- as.integer(
  present_control[
    taxon_match
  ]
)

results$Present_Case <- as.integer(
  present_case[
    taxon_match
  ]
)

results$All_Zero <- as.logical(
  all_zero_taxa[
    taxon_match
  ]
)


# ============================================================
# 21. MERGE GROUND TRUTH
# ============================================================

truth_match <- match(
  results$Taxon,
  rownames(truth)
)

if (any(is.na(truth_match))) {

  stop(
    "Ground-truth table could not be aligned."
  )
}

truth_aligned <- truth[
  truth_match,
  ,
  drop = FALSE
]

for (column_name in colnames(truth_aligned)) {

  results[
    ,
    column_name
  ] <- truth_aligned[
    ,
    column_name,
    drop = TRUE
  ]
}


# ============================================================
# 22. ADD ANALYSIS METADATA
# ============================================================

results$Alpha <- ALPHA
results$P_Adjust_Method <- P_ADJ_METHOD
results$Prevalence_Cut <- PREVALENCE_CUT
results$Library_Cut <- LIBRARY_CUT
results$Pseudo_Sensitivity <- TRUE
results$Conservative <- TRUE
results$S0_Percentile <- S0_PERC
results$Scenario <- "Null_Baseline_Low_Sparsity"
results$Seed <- RANDOM_SEED

results$ANCOMBC_Version <- as.character(
  packageVersion("ANCOMBC")
)


# ============================================================
# 23. SUMMARY STATISTICS
# ============================================================

n_available <- nrow(
  results
)

n_tested <- sum(
  results$Tested
)

n_not_tested <- sum(
  !results$Tested
)

n_standard_significant <- sum(
  results$Standard_Significant,
  na.rm = TRUE
)

n_robust_significant <- sum(
  results$Robust_Significant,
  na.rm = TRUE
)

minimum_raw_p <- min(
  results$P_Value,
  na.rm = TRUE
)

minimum_adjusted_p <- min(
  results$Adjusted_P_Value,
  na.rm = TRUE
)


# ============================================================
# 24. SORT RESULTS
# ============================================================

results_sorted <- results[
  order(
    results$Adjusted_P_Value,
    results$P_Value,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]


# ============================================================
# 25. SAVE STANDARDIZED RESULT
# ============================================================

standardized_file <- file.path(
  OUTPUT_DIR,
  "ancombc2_standardized_results.csv"
)

write.csv(
  results_sorted,
  standardized_file,
  row.names = FALSE
)


# ============================================================
# 26. SAVE SUMMARY
# ============================================================

summary_file <- file.path(
  OUTPUT_DIR,
  "ancombc2_summary.txt"
)

summary_lines <- c(
  "ANCOM-BC2 INDEPENDENT NULL VALIDATION",
  paste(rep("=", 60), collapse = ""),
  "",
  paste0(
    "ANCOMBC version: ",
    packageVersion("ANCOMBC")
  ),
  "Fixed formula: Group",
  "Reference group: Control",
  "Coefficient: GroupCase",
  "Dummy metadata column added: TRUE",
  paste0(
    "P-value adjustment: ",
    P_ADJ_METHOD
  ),
  paste0(
    "Prevalence cutoff: ",
    PREVALENCE_CUT
  ),
  paste0(
    "Library cutoff: ",
    LIBRARY_CUT
  ),
  "Pseudo-count sensitivity: TRUE",
  "Conservative sensitivity mode: TRUE",
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
    "Taxa tested: ",
    n_tested
  ),
  paste0(
    "Taxa not tested: ",
    n_not_tested
  ),
  paste0(
    "Standard BH-significant taxa: ",
    n_standard_significant
  ),
  paste0(
    "Robust significant taxa: ",
    n_robust_significant
  ),
  paste0(
    "Minimum raw p-value: ",
    signif(
      minimum_raw_p,
      10
    )
  ),
  paste0(
    "Minimum adjusted p-value: ",
    signif(
      minimum_adjusted_p,
      10
    )
  ),
  "",
  "Known true DA taxa: 0",
  "",
  paste0(
    "Any significant taxon in this dataset is ",
    "a false-positive discovery because the ",
    "simulation is a complete null."
  )
)

writeLines(
  summary_lines,
  summary_file
)


# ============================================================
# 27. PRINT SUMMARY
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
  "Taxa tested: ",
  n_tested,
  "\n",
  sep = ""
)

cat(
  "Taxa not tested: ",
  n_not_tested,
  "\n",
  sep = ""
)

cat(
  "Standard significant taxa at BH q < ",
  ALPHA,
  ": ",
  n_standard_significant,
  "\n",
  sep = ""
)

cat(
  "Robust significant taxa: ",
  n_robust_significant,
  "\n",
  sep = ""
)

cat(
  "Minimum raw p-value: ",
  signif(
    minimum_raw_p,
    10
  ),
  "\n",
  sep = ""
)

cat(
  "Minimum adjusted p-value: ",
  signif(
    minimum_adjusted_p,
    10
  ),
  "\n",
  sep = ""
)


# ============================================================
# 28. TAXON_033
# ============================================================

display_columns <- c(
  "Taxon",
  "Present_Samples",
  "Present_Control",
  "Present_Case",
  "Log_Fold_Change",
  "P_Value",
  "Adjusted_P_Value",
  "Standard_Significant",
  "Passed_Sensitivity",
  "Robust_Significant",
  "Effect_Direction"
)

cat(
  "\nTaxon_033 comparison:\n"
)

print(
  results[
    results$Taxon == "Taxon_033",
    display_columns,
    drop = FALSE
  ],
  row.names = FALSE
)


# ============================================================
# 29. TOP 15 TAXA
# ============================================================

cat(
  "\nTop 15 taxa by adjusted p-value:\n"
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
# 30. RAREST TAXA
# ============================================================

rarest_results <- results[
  order(
    results$Present_Samples,
    results$Taxon
  ),
  ,
  drop = FALSE
]

cat(
  "\nRarest taxa and ANCOM-BC2 results:\n"
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
# 31. OUTPUT FILES
# ============================================================

cat(
  "\nOutput files:\n"
)

cat(
  raw_result_file,
  "\n"
)

cat(
  standardized_file,
  "\n"
)

cat(
  summary_file,
  "\n"
)

cat(
  "\nANCOM-BC2 independent null validation completed.\n"
)
