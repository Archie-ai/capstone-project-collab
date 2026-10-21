# ============================================================
# ALDEx2 INDEPENDENT NULL VALIDATION
#
# Microbiome Abundance Consensus Machine
#
# Purpose:
# Independently analyze the validated simulated null dataset
# using ALDEx2 and compare its behavior with CLR/Wilcoxon.
#
# Dataset:
#   100 taxa
#   50 samples
#   25 Control
#   25 Case
#   0 true differentially abundant taxa
#
# ALDEx2 configuration:
#   Monte Carlo samples = 128
#   denominator = "all"
#   paired test = FALSE
#   random seed = 42
#   alpha = 0.05
#
# Both Welch and Wilcoxon results are retained.
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
  "aldex2"
)

MC_SAMPLES <- 128
DENOMINATOR <- "all"
ALPHA <- 0.05
RANDOM_SEED <- 42


# ============================================================
# 2. LOAD PACKAGE
# ============================================================

suppressPackageStartupMessages(
  library(ALDEx2)
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

separator <- paste(rep("=", 70), collapse = "")

cat("\n", separator, "\n", sep = "")
cat("ALDEx2 INDEPENDENT NULL VALIDATION\n")
cat(separator, "\n", sep = "")

cat("\nR version: ", R.version.string, "\n", sep = "")
cat(
  "ALDEx2 version: ",
  as.character(packageVersion("ALDEx2")),
  "\n",
  sep = ""
)


# ============================================================
# 3. LOAD INPUT FILES
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

required_columns <- c(
  "SampleID",
  "Group"
)

missing_columns <- setdiff(
  required_columns,
  colnames(metadata)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing metadata columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}


# ============================================================
# 5. VALIDATE SAMPLE IDENTIFIERS
# ============================================================

count_samples <- colnames(counts)
metadata_samples <- metadata$SampleID

if (!setequal(count_samples, metadata_samples)) {

  missing_from_metadata <- setdiff(
    count_samples,
    metadata_samples
  )

  missing_from_counts <- setdiff(
    metadata_samples,
    count_samples
  )

  stop(
    paste0(
      "Sample IDs do not match.\n",
      "Missing from metadata: ",
      paste(missing_from_metadata, collapse = ", "),
      "\n",
      "Missing from counts: ",
      paste(missing_from_counts, collapse = ", ")
    )
  )
}

metadata <- metadata[
  match(count_samples, metadata$SampleID),
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
  stop("Counts and metadata are not aligned.")
}


# ============================================================
# 6. DEFINE GROUPS
# ============================================================

conditions <- as.character(
  metadata$Group
)

control_count <- sum(
  conditions == "Control"
)

case_count <- sum(
  conditions == "Case"
)

cat("\nGroup counts:\n")
cat("Control: ", control_count, "\n", sep = "")
cat("Case:    ", case_count, "\n", sep = "")

if (control_count == 0) {
  stop("No Control samples found.")
}

if (case_count == 0) {
  stop("No Case samples found.")
}


# ============================================================
# 7. VALIDATE COUNT MATRIX
# ============================================================

counts_matrix <- as.matrix(
  counts
)

storage.mode(counts_matrix) <- "numeric"

if (any(!is.finite(counts_matrix))) {
  stop("Count matrix contains non-finite values.")
}

negative_count_total <- sum(
  counts_matrix < 0
)

non_integer_total <- sum(
  abs(
    counts_matrix - round(counts_matrix)
  ) > 1e-8
)

all_zero_taxa <- (
  rowSums(counts_matrix) == 0
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

if (negative_count_total > 0) {
  stop("Negative counts detected.")
}

if (non_integer_total > 0) {
  stop("Non-integer counts detected.")
}


# ALDEx2 requires integer counts.

counts_matrix <- round(
  counts_matrix
)

storage.mode(counts_matrix) <- "integer"


# ============================================================
# 8. CALCULATE TAXON PREVALENCE
# ============================================================

present_samples <- rowSums(
  counts_matrix > 0
)

present_control <- rowSums(
  counts_matrix[
    ,
    conditions == "Control",
    drop = FALSE
  ] > 0
)

present_case <- rowSums(
  counts_matrix[
    ,
    conditions == "Case",
    drop = FALSE
  ] > 0
)


# ============================================================
# 9. SET RANDOM SEED
# ============================================================

set.seed(
  RANDOM_SEED
)

cat("\nALDEx2 configuration:\n")
cat(
  "Monte Carlo samples: ",
  MC_SAMPLES,
  "\n",
  sep = ""
)

cat(
  "Denominator: ",
  DENOMINATOR,
  "\n",
  sep = ""
)

cat("Paired test: FALSE\n")

cat(
  "Random seed: ",
  RANDOM_SEED,
  "\n",
  sep = ""
)

cat(
  "Alpha: ",
  ALPHA,
  "\n",
  sep = ""
)


# ============================================================
# 10. RUN ALDEx2 CLR MONTE CARLO STEP
# ============================================================

cat(
  "\nGenerating ALDEx2 Monte Carlo CLR instances...\n"
)

aldex_clr <- aldex.clr(
  reads = counts_matrix,
  conds = conditions,
  mc.samples = MC_SAMPLES,
  denom = DENOMINATOR,
  verbose = FALSE,
  useMC = FALSE
)

cat(
  "Monte Carlo CLR generation completed.\n"
)


# ============================================================
# 11. RUN WELCH AND WILCOXON TESTS
# ============================================================

cat(
  "\nRunning ALDEx2 Welch and Wilcoxon tests...\n"
)

aldex_test <- aldex.ttest(
  aldex_clr,
  paired.test = FALSE,
  hist.plot = FALSE,
  verbose = FALSE
)

cat(
  "ALDEx2 testing completed.\n"
)


# ============================================================
# 12. CALCULATE ALDEx2 EFFECT SIZES
# ============================================================

cat(
  "\nCalculating ALDEx2 effect statistics...\n"
)

aldex_effect <- aldex.effect(
  aldex_clr,
  verbose = FALSE,
  include.sample.summary = FALSE,
  useMC = FALSE,
  CI = FALSE,
  paired.test = FALSE
)

cat(
  "ALDEx2 effect calculation completed.\n"
)


# ============================================================
# 13. CHECK OUTPUT DIMENSIONS
# ============================================================

cat("\nALDEx2 output dimensions:\n")

cat(
  "Test output taxa: ",
  nrow(aldex_test),
  "\n",
  sep = ""
)

cat(
  "Effect output taxa: ",
  nrow(aldex_effect),
  "\n",
  sep = ""
)

if (nrow(aldex_test) != nrow(counts_matrix)) {
  stop(
    paste0(
      "ALDEx2 test output has ",
      nrow(aldex_test),
      " taxa; expected ",
      nrow(counts_matrix),
      "."
    )
  )
}

if (nrow(aldex_effect) != nrow(counts_matrix)) {
  stop(
    paste0(
      "ALDEx2 effect output has ",
      nrow(aldex_effect),
      " taxa; expected ",
      nrow(counts_matrix),
      "."
    )
  )
}


# ============================================================
# 14. CREATE STANDARDIZED RESULT TABLE
# ============================================================

taxa <- rownames(
  aldex_test
)

results <- data.frame(
  Taxon = taxa,
  Method = "ALDEx2",
  Tested = TRUE,
  Welch_P_Value = aldex_test$we.ep,
  Welch_Adjusted_P_Value = aldex_test$we.eBH,
  Wilcoxon_P_Value = aldex_test$wi.ep,
  Wilcoxon_Adjusted_P_Value = aldex_test$wi.eBH,
  stringsAsFactors = FALSE
)


# ============================================================
# 15. ALIGN EFFECT-SIZE OUTPUT
# ============================================================

effect_order <- match(
  results$Taxon,
  rownames(aldex_effect)
)

if (any(is.na(effect_order))) {
  stop(
    "Unable to align ALDEx2 effect results with test results."
  )
}

effect_aligned <- aldex_effect[
  effect_order,
  ,
  drop = FALSE
]

results$ALDEx2_Effect <- effect_aligned$effect
results$ALDEx2_Diff_Between <- effect_aligned$diff.btw
results$ALDEx2_Diff_Within <- effect_aligned$diff.win
results$ALDEx2_Overlap <- effect_aligned$overlap


# ============================================================
# 16. DERIVE CASE-vs-CONTROL DIRECTION
# ============================================================

control_effect_column <- "rab.win.Control"
case_effect_column <- "rab.win.Case"

if (
  control_effect_column %in% colnames(effect_aligned) &&
  case_effect_column %in% colnames(effect_aligned)
) {

  results$Control_Median_CLR <- effect_aligned[
    ,
    control_effect_column
  ]

  results$Case_Median_CLR <- effect_aligned[
    ,
    case_effect_column
  ]

  results$Case_minus_Control_CLR <- (
    results$Case_Median_CLR -
    results$Control_Median_CLR
  )

  results$Effect_Direction <- ifelse(
    results$Case_minus_Control_CLR > 0,
    "Case_higher",
    ifelse(
      results$Case_minus_Control_CLR < 0,
      "Control_higher",
      "No_difference"
    )
  )

} else {

  results$Control_Median_CLR <- NA_real_
  results$Case_Median_CLR <- NA_real_
  results$Case_minus_Control_CLR <- NA_real_
  results$Effect_Direction <- NA_character_

  cat(
    "\nNOTE: Expected group-specific effect columns were not found.\n"
  )

  cat(
    "Available aldex.effect columns:\n"
  )

  print(
    colnames(effect_aligned)
  )
}


# ============================================================
# 17. ADD SIGNIFICANCE FLAGS
# ============================================================

results$Welch_Significant <- (
  !is.na(results$Welch_Adjusted_P_Value) &
  results$Welch_Adjusted_P_Value < ALPHA
)

results$Wilcoxon_Significant <- (
  !is.na(results$Wilcoxon_Adjusted_P_Value) &
  results$Wilcoxon_Adjusted_P_Value < ALPHA
)


# ============================================================
# 18. ADD PREVALENCE INFORMATION
# ============================================================

prevalence_order <- match(
  results$Taxon,
  rownames(counts_matrix)
)

if (any(is.na(prevalence_order))) {
  stop(
    "Unable to align taxon prevalence information."
  )
}

results$Present_Samples <- present_samples[
  prevalence_order
]

results$Present_Control <- present_control[
  prevalence_order
]

results$Present_Case <- present_case[
  prevalence_order
]

results$All_Zero <- all_zero_taxa[
  prevalence_order
]


# ============================================================
# 19. MERGE KNOWN GROUND TRUTH
# ============================================================

truth_order <- match(
  results$Taxon,
  rownames(truth)
)

if (any(is.na(truth_order))) {

  missing_truth_taxa <- results$Taxon[
    is.na(truth_order)
  ]

  stop(
    paste(
      "Missing ground truth for:",
      paste(
        missing_truth_taxa,
        collapse = ", "
      )
    )
  )
}

truth_aligned <- truth[
  truth_order,
  ,
  drop = FALSE
]

for (column_name in colnames(truth_aligned)) {
  results[[column_name]] <- truth_aligned[[column_name]]
}


# ============================================================
# 20. ADD ANALYSIS METADATA
# ============================================================

results$MC_Samples <- MC_SAMPLES
results$Denominator <- DENOMINATOR
results$Alpha <- ALPHA

results$Scenario <- (
  "Null_Baseline_Low_Sparsity"
)

results$Seed <- RANDOM_SEED

results$ALDEx2_Version <- as.character(
  packageVersion("ALDEx2")
)


# ============================================================
# 21. CALCULATE SUMMARY VALUES
# ============================================================

n_taxa <- nrow(
  results
)

welch_na <- sum(
  is.na(
    results$Welch_Adjusted_P_Value
  )
)

wilcoxon_na <- sum(
  is.na(
    results$Wilcoxon_Adjusted_P_Value
  )
)

effect_na <- sum(
  is.na(
    results$ALDEx2_Effect
  )
)

n_welch_tested <- sum(
  !is.na(
    results$Welch_Adjusted_P_Value
  )
)

n_wilcoxon_tested <- sum(
  !is.na(
    results$Wilcoxon_Adjusted_P_Value
  )
)

n_welch_significant <- sum(
  results$Welch_Significant,
  na.rm = TRUE
)

n_wilcoxon_significant <- sum(
  results$Wilcoxon_Significant,
  na.rm = TRUE
)

minimum_welch_p <- min(
  results$Welch_P_Value,
  na.rm = TRUE
)

minimum_welch_bh <- min(
  results$Welch_Adjusted_P_Value,
  na.rm = TRUE
)

minimum_wilcoxon_p <- min(
  results$Wilcoxon_P_Value,
  na.rm = TRUE
)

minimum_wilcoxon_bh <- min(
  results$Wilcoxon_Adjusted_P_Value,
  na.rm = TRUE
)


# ============================================================
# 22. SORT RESULTS
# ============================================================

welch_results <- results[
  order(
    results$Welch_Adjusted_P_Value,
    results$Welch_P_Value,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

wilcoxon_results <- results[
  order(
    results$Wilcoxon_Adjusted_P_Value,
    results$Wilcoxon_P_Value,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]


# ============================================================
# 23. SAVE RESULTS
# ============================================================

full_results_file <- file.path(
  OUTPUT_DIR,
  "aldex2_results.csv"
)

welch_results_file <- file.path(
  OUTPUT_DIR,
  "aldex2_welch_results.csv"
)

wilcoxon_results_file <- file.path(
  OUTPUT_DIR,
  "aldex2_wilcoxon_results.csv"
)

write.csv(
  results,
  full_results_file,
  row.names = FALSE
)

write.csv(
  welch_results,
  welch_results_file,
  row.names = FALSE
)

write.csv(
  wilcoxon_results,
  wilcoxon_results_file,
  row.names = FALSE
)


# ============================================================
# 24. CREATE TEXT SUMMARY
# ============================================================

summary_file <- file.path(
  OUTPUT_DIR,
  "aldex2_summary.txt"
)

summary_lines <- c(
  "ALDEx2 INDEPENDENT NULL VALIDATION",
  paste(rep("=", 60), collapse = ""),
  "",
  paste0(
    "ALDEx2 version: ",
    packageVersion("ALDEx2")
  ),
  paste0(
    "Monte Carlo samples: ",
    MC_SAMPLES
  ),
  paste0(
    "Denominator: ",
    DENOMINATOR
  ),
  paste0(
    "Random seed: ",
    RANDOM_SEED
  ),
  "",
  paste0(
    "Taxa available: ",
    n_taxa
  ),
  paste0(
    "Control samples: ",
    control_count
  ),
  paste0(
    "Case samples: ",
    case_count
  ),
  paste0(
    "All-zero taxa: ",
    sum(all_zero_taxa)
  ),
  "",
  "WELCH",
  paste0(
    "Taxa with usable BH values: ",
    n_welch_tested
  ),
  paste0(
    "Missing BH values: ",
    welch_na
  ),
  paste0(
    "Significant taxa: ",
    n_welch_significant
  ),
  paste0(
    "Minimum raw p-value: ",
    signif(minimum_welch_p, 10)
  ),
  paste0(
    "Minimum adjusted p-value: ",
    signif(minimum_welch_bh, 10)
  ),
  "",
  "WILCOXON",
  paste0(
    "Taxa with usable BH values: ",
    n_wilcoxon_tested
  ),
  paste0(
    "Missing BH values: ",
    wilcoxon_na
  ),
  paste0(
    "Significant taxa: ",
    n_wilcoxon_significant
  ),
  paste0(
    "Minimum raw p-value: ",
    signif(minimum_wilcoxon_p, 10)
  ),
  paste0(
    "Minimum adjusted p-value: ",
    signif(minimum_wilcoxon_bh, 10)
  ),
  "",
  paste0(
    "Missing effect estimates: ",
    effect_na
  ),
  "",
  "Known true DA taxa: 0",
  "",
  paste0(
    "Any significant result in this dataset ",
    "is a false-positive discovery because ",
    "the simulation is a global null."
  )
)

writeLines(
  summary_lines,
  summary_file
)


# ============================================================
# 25. PRINT MAIN SUMMARY
# ============================================================

cat("\n", separator, "\n", sep = "")
cat("RESULT SUMMARY\n")
cat(separator, "\n", sep = "")

cat(
  "Taxa available: ",
  n_taxa,
  "\n",
  sep = ""
)

cat("\nWELCH TEST\n")

cat(
  "Taxa with usable BH values: ",
  n_welch_tested,
  "\n",
  sep = ""
)

cat(
  "Missing BH values: ",
  welch_na,
  "\n",
  sep = ""
)

cat(
  "Significant taxa at BH q < ",
  ALPHA,
  ": ",
  n_welch_significant,
  "\n",
  sep = ""
)

cat(
  "Minimum raw p-value: ",
  signif(minimum_welch_p, 10),
  "\n",
  sep = ""
)

cat(
  "Minimum adjusted p-value: ",
  signif(minimum_welch_bh, 10),
  "\n",
  sep = ""
)


cat("\nWILCOXON TEST\n")

cat(
  "Taxa with usable BH values: ",
  n_wilcoxon_tested,
  "\n",
  sep = ""
)

cat(
  "Missing BH values: ",
  wilcoxon_na,
  "\n",
  sep = ""
)

cat(
  "Significant taxa at BH q < ",
  ALPHA,
  ": ",
  n_wilcoxon_significant,
  "\n",
  sep = ""
)

cat(
  "Minimum raw p-value: ",
  signif(minimum_wilcoxon_p, 10),
  "\n",
  sep = ""
)

cat(
  "Minimum adjusted p-value: ",
  signif(minimum_wilcoxon_bh, 10),
  "\n",
  sep = ""
)


# ============================================================
# 26. INSPECT TAXON_033
# ============================================================

cat(
  "\nTaxon_033 comparison:\n"
)

taxon_033_columns <- c(
  "Taxon",
  "Present_Samples",
  "Present_Control",
  "Present_Case",
  "Welch_P_Value",
  "Welch_Adjusted_P_Value",
  "Welch_Significant",
  "Wilcoxon_P_Value",
  "Wilcoxon_Adjusted_P_Value",
  "Wilcoxon_Significant",
  "ALDEx2_Effect",
  "Effect_Direction"
)

taxon_033 <- results[
  results$Taxon == "Taxon_033",
  taxon_033_columns,
  drop = FALSE
]

print(
  taxon_033,
  row.names = FALSE
)


# ============================================================
# 27. PRINT TOP 15 WELCH RESULTS
# ============================================================

cat(
  "\nTop 15 taxa by Welch adjusted p-value:\n"
)

welch_display_columns <- c(
  "Taxon",
  "Present_Samples",
  "Present_Control",
  "Present_Case",
  "Welch_P_Value",
  "Welch_Adjusted_P_Value",
  "Welch_Significant",
  "ALDEx2_Effect",
  "Effect_Direction"
)

print(
  head(
    welch_results[
      ,
      welch_display_columns,
      drop = FALSE
    ],
    15
  ),
  row.names = FALSE
)


# ============================================================
# 28. PRINT TOP 15 WILCOXON RESULTS
# ============================================================

cat(
  "\nTop 15 taxa by Wilcoxon adjusted p-value:\n"
)

wilcoxon_display_columns <- c(
  "Taxon",
  "Present_Samples",
  "Present_Control",
  "Present_Case",
  "Wilcoxon_P_Value",
  "Wilcoxon_Adjusted_P_Value",
  "Wilcoxon_Significant",
  "ALDEx2_Effect",
  "Effect_Direction"
)

print(
  head(
    wilcoxon_results[
      ,
      wilcoxon_display_columns,
      drop = FALSE
    ],
    15
  ),
  row.names = FALSE
)


# ============================================================
# 29. PRINT RAREST TAXA
# ============================================================

cat(
  "\nRarest taxa and ALDEx2 results:\n"
)

rarest_results <- results[
  order(
    results$Present_Samples,
    results$Taxon
  ),
  ,
  drop = FALSE
]

rare_display_columns <- c(
  "Taxon",
  "Present_Samples",
  "Present_Control",
  "Present_Case",
  "Welch_Adjusted_P_Value",
  "Wilcoxon_Adjusted_P_Value",
  "ALDEx2_Effect",
  "Effect_Direction"
)

print(
  head(
    rarest_results[
      ,
      rare_display_columns,
      drop = FALSE
    ],
    10
  ),
  row.names = FALSE
)


# ============================================================
# 30. PRINT OUTPUT LOCATIONS
# ============================================================

cat("\nOutput files:\n")
cat(full_results_file, "\n")
cat(welch_results_file, "\n")
cat(wilcoxon_results_file, "\n")
cat(summary_file, "\n")

cat(
  "\nALDEx2 independent null validation completed.\n"
)
