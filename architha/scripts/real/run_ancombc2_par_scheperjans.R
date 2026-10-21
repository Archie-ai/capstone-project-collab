#!/usr/bin/env Rscript

suppressPackageStartupMessages(
    library(ANCOMBC)
)

# ============================================================
# CONFIGURATION
# ============================================================

DATASET <- "par_scheperjans"

ALPHA <- 0.05
SEED <- 42

COUNTS_FILE <- "data/processed/par_scheperjans/par_scheperjans_counts.csv"
METADATA_FILE <- "data/processed/par_scheperjans/par_scheperjans_metadata.csv"
TAXONOMY_FILE <- "data/processed/par_scheperjans/par_scheperjans_taxonomy_mapping.csv"

OUTPUT_DIR <- "results/real/par_scheperjans/ancombc2"

dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

RESULT_FILE <- file.path(
    OUTPUT_DIR,
    "ancombc2_results.csv"
)

SIGNIFICANT_FILE <- file.path(
    OUTPUT_DIR,
    "ancombc2_significant_standard.csv"
)

ROBUST_FILE <- file.path(
    OUTPUT_DIR,
    "ancombc2_significant_robust.csv"
)

SUMMARY_FILE <- file.path(
    OUTPUT_DIR,
    "ancombc2_summary.txt"
)

RAW_RESULT_FILE <- file.path(
    OUTPUT_DIR,
    "ancombc2_raw_primary_result.csv"
)


# ============================================================
# HEADER
# ============================================================

cat(strrep("=", 72), "\n", sep = "")
cat("PAR_SCHEPERJANS ANCOM-BC2 ANALYSIS\n")
cat(strrep("=", 72), "\n", sep = "")

cat(
    "\nR version: ",
    as.character(getRversion()),
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
# LOAD CANONICAL DATA
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

taxonomy <- read.csv(
    TAXONOMY_FILE,
    check.names = FALSE,
    stringsAsFactors = FALSE
)


cat("\nCanonical input dimensions:\n")

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
    " samples x ",
    ncol(metadata),
    " columns\n",
    sep = ""
)

cat(
    "Taxonomy: ",
    nrow(taxonomy),
    " taxa x ",
    ncol(taxonomy),
    " columns\n",
    sep = ""
)


# ============================================================
# BASIC VALIDATION
# ============================================================

if (!all(c("SampleID", "Group") %in% colnames(metadata))) {
    stop("Metadata must contain SampleID and Group.")
}

if (!"Taxon" %in% colnames(taxonomy)) {
    stop("Taxonomy mapping must contain Taxon.")
}


metadata$SampleID <- trimws(
    as.character(metadata$SampleID)
)

metadata$Group <- trimws(
    as.character(metadata$Group)
)

colnames(counts) <- trimws(
    colnames(counts)
)

rownames(counts) <- trimws(
    rownames(counts)
)

taxonomy$Taxon <- trimws(
    as.character(taxonomy$Taxon)
)


if (anyDuplicated(metadata$SampleID) > 0) {
    stop("Duplicate metadata SampleIDs detected.")
}

if (anyDuplicated(colnames(counts)) > 0) {
    stop("Duplicate count-table sample IDs detected.")
}

if (anyDuplicated(rownames(counts)) > 0) {
    stop("Duplicate Taxon IDs detected.")
}


# ============================================================
# SAMPLE ALIGNMENT
# ============================================================

if (!setequal(
    colnames(counts),
    metadata$SampleID
)) {
    stop("Count-table and metadata sample sets do not match.")
}


metadata_index <- match(
    colnames(counts),
    metadata$SampleID
)

if (any(is.na(metadata_index))) {
    stop("Metadata alignment failed.")
}

metadata <- metadata[
    metadata_index,
    ,
    drop = FALSE
]


alignment_pass <- identical(
    colnames(counts),
    metadata$SampleID
)


cat("\nSample alignment:\n")

cat(
    "Counts and metadata aligned: ",
    alignment_pass,
    "\n",
    sep = ""
)

if (!alignment_pass) {
    stop("Sample alignment failed.")
}


# ============================================================
# GROUP VALIDATION
# ============================================================

n_control <- sum(
    metadata$Group == "Control"
)

n_case <- sum(
    metadata$Group == "Case"
)


cat("\nGroup counts:\n")

cat(
    "Control: ",
    n_control,
    "\n",
    sep = ""
)

cat(
    "Case:    ",
    n_case,
    "\n",
    sep = ""
)


if (n_control != 74) {
    stop("Expected 74 Control samples.")
}

if (n_case != 74) {
    stop("Expected 74 Case samples.")
}

if (!setequal(
    unique(metadata$Group),
    c("Control", "Case")
)) {
    stop("Unexpected group labels detected.")
}


# ============================================================
# COUNT MATRIX VALIDATION
# ============================================================

counts_matrix <- as.matrix(
    counts
)

storage.mode(
    counts_matrix
) <- "numeric"


if (any(is.na(counts_matrix))) {
    stop("Missing/non-numeric counts detected.")
}

if (any(!is.finite(counts_matrix))) {
    stop("Non-finite counts detected.")
}


negative_cells <- sum(
    counts_matrix < 0
)

non_integer_cells <- sum(
    abs(
        counts_matrix -
        round(counts_matrix)
    ) > 1e-8
)

all_zero_taxa <- sum(
    rowSums(counts_matrix) == 0
)

zero_library_samples <- sum(
    colSums(counts_matrix) == 0
)

zero_fraction <- mean(
    counts_matrix == 0
)


cat("\nCount QC:\n")

cat(
    "Negative cells: ",
    negative_cells,
    "\n",
    sep = ""
)

cat(
    "Non-integer cells: ",
    non_integer_cells,
    "\n",
    sep = ""
)

cat(
    "All-zero taxa: ",
    all_zero_taxa,
    "\n",
    sep = ""
)

cat(
    "Zero-library samples: ",
    zero_library_samples,
    "\n",
    sep = ""
)

cat(
    "Zero fraction: ",
    sprintf("%.6f", zero_fraction),
    "\n",
    sep = ""
)


if (negative_cells != 0) {
    stop("Negative counts detected.")
}

if (non_integer_cells != 0) {
    stop("Non-integer counts detected.")
}

if (all_zero_taxa != 0) {
    stop("All-zero taxa detected.")
}

if (zero_library_samples != 0) {
    stop("Zero-library samples detected.")
}


counts_matrix <- round(
    counts_matrix
)


# ============================================================
# ANCOM-BC2 METADATA
#
# Control is explicitly the reference group.
#
# Dummy is retained as a harmless constant metadata field.
# It is NOT included in the statistical formula.
# ============================================================

meta_ancom <- data.frame(

    Group = factor(
        metadata$Group,
        levels = c(
            "Control",
            "Case"
        )
    ),

    Dummy = rep(
        1,
        nrow(metadata)
    ),

    row.names = metadata$SampleID,

    stringsAsFactors = FALSE
)


if (!identical(
    colnames(counts_matrix),
    rownames(meta_ancom)
)) {
    stop("ANCOM-BC2 metadata rownames are not aligned.")
}


cat("\nANCOM-BC2 settings:\n")

cat("Reference group: Control\n")
cat("Target coefficient: GroupCase\n")
cat("Formula: Group\n")
cat("P-value adjustment: BH\n")
cat("External prevalence filter already applied: 10%\n")
cat("Internal prevalence filter: 0\n")
cat("Internal library-size filter: 0\n")
cat("Pseudo-count sensitivity analysis: TRUE\n")
cat("Structural-zero analysis: FALSE\n")
cat("s0 percentile: 0.05\n")

cat(
    "Alpha: ",
    ALPHA,
    "\n",
    sep = ""
)

cat(
    "Seed: ",
    SEED,
    "\n",
    sep = ""
)


# ============================================================
# RUN ANCOM-BC2
# ============================================================

set.seed(
    SEED
)

cat("\nRunning ancombc2()...\n")


ancom_out <- ancombc2(

    data = counts_matrix,

    taxa_are_rows = TRUE,

    aggregate_data = counts_matrix,

    meta_data = meta_ancom,

    fix_formula = "Group",

    rand_formula = NULL,

    p_adj_method = "BH",

    pseudo = 0,

    pseudo_sens = TRUE,

    prv_cut = 0,

    lib_cut = 0,

    s0_perc = 0.05,

    group = NULL,

    struc_zero = FALSE,

    neg_lb = FALSE,

    alpha = ALPHA,

    n_cl = 1,

    verbose = FALSE,

    global = FALSE,

    pairwise = FALSE,

    dunnet = FALSE,

    trend = FALSE,

    iter_control = list(
        tol = 1e-2,
        max_iter = 20,
        verbose = FALSE
    ),

    em_control = list(
        tol = 1e-5,
        max_iter = 100
    ),

    lme_control = NULL,

    mdfdr_control = list(
        fwer_ctrl_method = "holm",
        B = 100
    ),

    trend_control = NULL
)


cat("ancombc2() complete.\n")


# ============================================================
# EXTRACT PRIMARY RESULT
# ============================================================

if (is.null(ancom_out$res)) {
    stop("ANCOM-BC2 output does not contain res.")
}


res <- as.data.frame(
    ancom_out$res,
    stringsAsFactors = FALSE
)


cat("\nANCOM-BC2 primary result dimensions:\n")

cat(
    nrow(res),
    " rows x ",
    ncol(res),
    " columns\n",
    sep = ""
)


cat(
    "Result columns: ",
    paste(
        colnames(res),
        collapse = ", "
    ),
    "\n",
    sep = ""
)


# ============================================================
# TAXON IDENTIFIER
# ============================================================

if ("taxon" %in% colnames(res)) {

    res$Taxon <- as.character(
        res$taxon
    )

} else if (!is.null(rownames(res))) {

    res$Taxon <- rownames(
        res
    )

} else {

    stop(
        "Unable to determine Taxon IDs from ANCOM-BC2 result."
    )
}


# ============================================================
# REQUIRED GROUPCASE COLUMNS
# ============================================================

expected_columns <- c(
    "lfc_GroupCase",
    "se_GroupCase",
    "W_GroupCase",
    "p_GroupCase",
    "q_GroupCase",
    "diff_GroupCase",
    "passed_ss_GroupCase"
)


missing_result_columns <- setdiff(
    expected_columns,
    colnames(res)
)


if (length(missing_result_columns) > 0) {

    stop(
        paste0(
            "Missing expected GroupCase columns: ",
            paste(
                missing_result_columns,
                collapse = ", "
            ),
            "\nAvailable columns: ",
            paste(
                colnames(res),
                collapse = ", "
            )
        )
    )
}


# ============================================================
# ALIGN TO CANONICAL TAXA
# ============================================================

canonical_taxa <- rownames(
    counts_matrix
)


result_index <- match(
    canonical_taxa,
    res$Taxon
)

taxonomy_index <- match(
    canonical_taxa,
    taxonomy$Taxon
)


if (any(is.na(result_index))) {

    missing_taxa <- canonical_taxa[
        is.na(result_index)
    ]

    stop(
        paste0(
            "Taxa missing from ANCOM-BC2 result: ",
            paste(
                head(missing_taxa, 20),
                collapse = ", "
            )
        )
    )
}


if (any(is.na(taxonomy_index))) {
    stop("Taxa missing from taxonomy mapping.")
}


res <- res[
    result_index,
    ,
    drop = FALSE
]


taxonomy_aligned <- taxonomy[
    taxonomy_index,
    ,
    drop = FALSE
]


# ============================================================
# SIGNIFICANCE
#
# Standard:
#   q < 0.05
#
# Robust:
#   significant AND passed pseudo-count sensitivity
# ============================================================

standard_significant <- (
    !is.na(res$q_GroupCase) &
    res$q_GroupCase < ALPHA
)


robust_significant <- (
    standard_significant &
    !is.na(res$passed_ss_GroupCase) &
    res$passed_ss_GroupCase
)


effect_direction <- ifelse(
    res$lfc_GroupCase > 0,
    "Case_higher",
    ifelse(
        res$lfc_GroupCase < 0,
        "Control_higher",
        "No_difference"
    )
)


# ============================================================
# PREVALENCE
# ============================================================

control_ids <- metadata$SampleID[
    metadata$Group == "Control"
]

case_ids <- metadata$SampleID[
    metadata$Group == "Case"
]


present_total <- rowSums(
    counts_matrix > 0
)

present_control <- rowSums(
    counts_matrix[
        ,
        control_ids,
        drop = FALSE
    ] > 0
)

present_case <- rowSums(
    counts_matrix[
        ,
        case_ids,
        drop = FALSE
    ] > 0
)


# ============================================================
# STANDARDIZED OUTPUT
# ============================================================

results <- data.frame(

    Taxon = canonical_taxa,

    Method = "ANCOM_BC2",

    Tested = !is.na(
        res$p_GroupCase
    ),

    P_Value =
        res$p_GroupCase,

    Adjusted_P_Value =
        res$q_GroupCase,

    Significant =
        robust_significant,

    Standard_Significant =
        standard_significant,

    Passed_Pseudocount_Sensitivity =
        res$passed_ss_GroupCase,

    Robust_Significant =
        robust_significant,

    Effect_Estimate =
        res$lfc_GroupCase,

    Effect_Metric =
        "ANCOM-BC2 natural-log LFC (Case vs Control)",

    Effect_Direction =
        effect_direction,

    Standard_Error =
        res$se_GroupCase,

    Test_Statistic =
        res$W_GroupCase,

    ANCOM_Diff_Call =
        res$diff_GroupCase,

    Present_Samples =
        as.integer(present_total),

    Present_Control =
        as.integer(present_control),

    Present_Case =
        as.integer(present_case),

    Alpha = ALPHA,

    P_Adjustment = "BH",

    Pseudo_Sensitivity = TRUE,

    Internal_Prevalence_Cutoff = 0,

    External_Prevalence_Cutoff = 0.10,

    Dataset = DATASET,

    Scenario = "Real_Data",

    Denovo_ID =
        taxonomy_aligned$Denovo_ID,

    Kingdom =
        taxonomy_aligned$Kingdom,

    Phylum =
        taxonomy_aligned$Phylum,

    Class =
        taxonomy_aligned$Class,

    Order =
        taxonomy_aligned$Order,

    Family =
        taxonomy_aligned$Family,

    Genus =
        taxonomy_aligned$Genus,

    Species =
        taxonomy_aligned$Species,

    Source_Taxon =
        taxonomy_aligned$Source_Taxon,

    stringsAsFactors = FALSE,
    check.names = FALSE
)


# ============================================================
# SORT
# ============================================================

sort_q <- results$Adjusted_P_Value
sort_p <- results$P_Value

sort_q[
    is.na(sort_q)
] <- Inf

sort_p[
    is.na(sort_p)
] <- Inf


sort_index <- order(
    sort_q,
    sort_p,
    results$Taxon
)


results <- results[
    sort_index,
    ,
    drop = FALSE
]

row.names(results) <- NULL


standard_results <- results[
    results$Standard_Significant,
    ,
    drop = FALSE
]


robust_results <- results[
    results$Robust_Significant,
    ,
    drop = FALSE
]


# ============================================================
# SUMMARY STATISTICS
# ============================================================

n_available <- nrow(
    results
)

n_tested <- sum(
    results$Tested,
    na.rm = TRUE
)

n_standard <- sum(
    results$Standard_Significant,
    na.rm = TRUE
)

n_passed_ss <- sum(
    results$Passed_Pseudocount_Sensitivity,
    na.rm = TRUE
)

n_robust <- sum(
    results$Robust_Significant,
    na.rm = TRUE
)


robust_case_higher <- sum(
    results$Robust_Significant &
    results$Effect_Direction == "Case_higher",
    na.rm = TRUE
)

robust_control_higher <- sum(
    results$Robust_Significant &
    results$Effect_Direction == "Control_higher",
    na.rm = TRUE
)


min_p <- min(
    results$P_Value,
    na.rm = TRUE
)

min_q <- min(
    results$Adjusted_P_Value,
    na.rm = TRUE
)


# ============================================================
# PRINT SUMMARY
# ============================================================

cat(
    "\n",
    strrep("=", 72),
    "\n",
    sep = ""
)

cat("RESULT SUMMARY\n")

cat(
    strrep("=", 72),
    "\n",
    sep = ""
)


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
    "Standard significant taxa (q < 0.05): ",
    n_standard,
    "\n",
    sep = ""
)

cat(
    "Taxa passing pseudo-count sensitivity: ",
    n_passed_ss,
    "\n",
    sep = ""
)

cat(
    "Robust significant taxa: ",
    n_robust,
    "\n",
    sep = ""
)

cat(
    "Robust Case-higher taxa: ",
    robust_case_higher,
    "\n",
    sep = ""
)

cat(
    "Robust Control-higher taxa: ",
    robust_control_higher,
    "\n",
    sep = ""
)

cat(
    "Minimum raw p-value: ",
    format(
        min_p,
        digits = 10
    ),
    "\n",
    sep = ""
)

cat(
    "Minimum adjusted p-value: ",
    format(
        min_q,
        digits = 10
    ),
    "\n",
    sep = ""
)


# ============================================================
# TOP 20
# ============================================================

cat(
    "\nTop 20 taxa by adjusted p-value:\n"
)


print(
    head(
        results[
            ,
            c(
                "Taxon",
                "Denovo_ID",
                "Present_Samples",
                "Present_Control",
                "Present_Case",
                "Effect_Estimate",
                "Effect_Direction",
                "P_Value",
                "Adjusted_P_Value",
                "Standard_Significant",
                "Passed_Pseudocount_Sensitivity",
                "Robust_Significant",
                "Phylum",
                "Family",
                "Genus"
            ),
            drop = FALSE
        ],
        20
    ),
    row.names = FALSE
)


# ============================================================
# SAVE RESULTS
# ============================================================

write.csv(
    res,
    RAW_RESULT_FILE,
    row.names = FALSE,
    na = ""
)

write.csv(
    results,
    RESULT_FILE,
    row.names = FALSE,
    na = ""
)

write.csv(
    standard_results,
    SIGNIFICANT_FILE,
    row.names = FALSE,
    na = ""
)

write.csv(
    robust_results,
    ROBUST_FILE,
    row.names = FALSE,
    na = ""
)


# ============================================================
# SUMMARY FILE
# ============================================================

summary_lines <- c(

    "PAR_SCHEPERJANS ANCOM-BC2 SUMMARY",

    strrep("=", 68),

    "",

    paste0(
        "Taxa analyzed: ",
        n_available
    ),

    paste0(
        "Control samples: ",
        n_control
    ),

    paste0(
        "Case samples: ",
        n_case
    ),

    paste0(
        "R version: ",
        as.character(getRversion())
    ),

    paste0(
        "ANCOMBC version: ",
        as.character(
            packageVersion("ANCOMBC")
        )
    ),

    "",

    "SETTINGS",

    "Formula: Group",

    "Reference group: Control",

    "Target coefficient: GroupCase",

    "P-value adjustment: BH",

    "Pseudo-count sensitivity: TRUE",

    "Internal prevalence cutoff: 0",

    "External prevalence cutoff: 0.10",

    "Internal library cutoff: 0",

    "Structural-zero detection: FALSE",

    "s0 percentile: 0.05",

    paste0(
        "Alpha: ",
        ALPHA
    ),

    "",

    "RESULTS",

    paste0(
        "Taxa tested: ",
        n_tested
    ),

    paste0(
        "Standard significant taxa: ",
        n_standard
    ),

    paste0(
        "Robust significant taxa: ",
        n_robust
    ),

    paste0(
        "Robust Case-higher taxa: ",
        robust_case_higher
    ),

    paste0(
        "Robust Control-higher taxa: ",
        robust_control_higher
    ),

    paste0(
        "Minimum raw p-value: ",
        format(
            min_p,
            digits = 10
        )
    ),

    paste0(
        "Minimum adjusted p-value: ",
        format(
            min_q,
            digits = 10
        )
    ),

    "",

    "PRIMARY CROSS-METHOD CALL",

    paste0(
        "Robust significance = q < 0.05 ",
        "AND passed pseudo-count sensitivity."
    ),

    "",

    "EFFECT",

    paste0(
        "ANCOM-BC2 LFC is on the natural-log scale. ",
        "Positive GroupCase LFC = Case higher."
    ),

    "",

    "INTERPRETATION",

    paste0(
        "Real-data findings have unknown ground truth ",
        "and are not classified as true positives."
    )
)


writeLines(
    summary_lines,
    SUMMARY_FILE
)


# ============================================================
# COMPLETE
# ============================================================

cat("\nOutput files:\n")

cat(
    RESULT_FILE,
    "\n"
)

cat(
    SIGNIFICANT_FILE,
    "\n"
)

cat(
    ROBUST_FILE,
    "\n"
)

cat(
    RAW_RESULT_FILE,
    "\n"
)

cat(
    SUMMARY_FILE,
    "\n"
)


cat(
    "\n",
    strrep("=", 72),
    "\n",
    sep = ""
)

cat(
    "ANCOM-BC2 REAL-DATA ANALYSIS: PASS\n"
)

cat(
    strrep("=", 72),
    "\n",
    sep = ""
)
