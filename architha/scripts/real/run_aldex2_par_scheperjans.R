#!/usr/bin/env Rscript

suppressPackageStartupMessages(
    library(ALDEx2)
)

# ============================================================
# CONFIGURATION
# ============================================================

DATASET <- "par_scheperjans"
ALPHA <- 0.05
MC_SAMPLES <- 128
SEED <- 42

COUNTS_FILE <- "data/processed/par_scheperjans/par_scheperjans_counts.csv"
METADATA_FILE <- "data/processed/par_scheperjans/par_scheperjans_metadata.csv"
TAXONOMY_FILE <- "data/processed/par_scheperjans/par_scheperjans_taxonomy_mapping.csv"

OUTPUT_DIR <- "results/real/par_scheperjans/aldex2"

dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

RESULT_FILE <- file.path(
    OUTPUT_DIR,
    "aldex2_results.csv"
)

WILCOXON_FILE <- file.path(
    OUTPUT_DIR,
    "aldex2_wilcoxon_results.csv"
)

WELCH_FILE <- file.path(
    OUTPUT_DIR,
    "aldex2_welch_results.csv"
)

SIGNIFICANT_FILE <- file.path(
    OUTPUT_DIR,
    "aldex2_significant_wilcoxon.csv"
)

SUMMARY_FILE <- file.path(
    OUTPUT_DIR,
    "aldex2_summary.txt"
)

# ============================================================
# HEADER
# ============================================================

cat(strrep("=", 72), "\n", sep = "")
cat("PAR_SCHEPERJANS ALDEX2 ANALYSIS\n")
cat(strrep("=", 72), "\n", sep = "")

cat(
    "\nR version: ",
    as.character(getRversion()),
    "\n",
    sep = ""
)

cat(
    "ALDEx2 version: ",
    as.character(packageVersion("ALDEx2")),
    "\n",
    sep = ""
)

# ============================================================
# LOAD DATA
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
# VALIDATE REQUIRED COLUMNS
# ============================================================

if (!all(c("SampleID", "Group") %in% colnames(metadata))) {
    stop("Metadata must contain SampleID and Group columns.")
}

required_taxonomy <- c(
    "Taxon",
    "Denovo_ID",
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "Source_Taxon"
)

if (!all(required_taxonomy %in% colnames(taxonomy))) {
    stop("Taxonomy mapping is missing required columns.")
}

# ============================================================
# CLEAN IDS
# ============================================================

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

# ============================================================
# DUPLICATE CHECKS
# ============================================================

if (anyDuplicated(metadata$SampleID) > 0) {
    stop("Duplicate metadata SampleIDs detected.")
}

if (anyDuplicated(colnames(counts)) > 0) {
    stop("Duplicate count-table sample IDs detected.")
}

if (anyDuplicated(rownames(counts)) > 0) {
    stop("Duplicate Taxon IDs detected.")
}

if (anyDuplicated(taxonomy$Taxon) > 0) {
    stop("Duplicate taxonomy Taxon IDs detected.")
}

# ============================================================
# ALIGN METADATA
# ============================================================

if (!setequal(colnames(counts), metadata$SampleID)) {
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
# GROUP CHECK
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
# COUNT MATRIX
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
    "Matrix class: ",
    class(counts_matrix)[1],
    "\n",
    sep = ""
)

cat(
    "Matrix storage mode: ",
    storage.mode(counts_matrix),
    "\n",
    sep = ""
)

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
# CONDITIONS
#
# IMPORTANT:
# ALDEx2 receives a CHARACTER VECTOR, not a factor.
# ============================================================

conditions <- as.character(
    metadata$Group
)

if (!is.character(conditions)) {
    stop("ALDEx2 conditions are not character.")
}

if (length(conditions) != ncol(counts_matrix)) {
    stop("Condition vector length does not match sample count.")
}

cat("\nALDEx2 condition QC:\n")

cat(
    "Conditions class: ",
    class(conditions),
    "\n",
    sep = ""
)

cat(
    "Condition count: ",
    length(conditions),
    "\n",
    sep = ""
)

cat(
    "Unique conditions: ",
    paste(
        sort(unique(conditions)),
        collapse = ", "
    ),
    "\n",
    sep = ""
)

# ============================================================
# SETTINGS
# ============================================================

cat("\nALDEx2 settings:\n")

cat(
    "Monte Carlo samples: ",
    MC_SAMPLES,
    "\n",
    sep = ""
)

cat("Denominator: all\n")
cat("Paired test: FALSE\n")
cat("Primary statistic: wi.eBH\n")

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
# RUN ALDEX2
# ============================================================

set.seed(
    SEED
)

cat("\nRunning aldex.clr()...\n")

aldex_clr <- aldex.clr(
    reads = counts_matrix,
    conds = conditions,
    mc.samples = MC_SAMPLES,
    denom = "all",
    verbose = FALSE
)

cat("aldex.clr() complete.\n")

cat("Running aldex.ttest()...\n")

aldex_test <- aldex.ttest(
    aldex_clr,
    paired.test = FALSE,
    hist.plot = FALSE,
    verbose = FALSE
)

cat("aldex.ttest() complete.\n")

cat("Running aldex.effect()...\n")

aldex_effect <- aldex.effect(
    aldex_clr,
    paired.test = FALSE,
    verbose = FALSE
)

cat("aldex.effect() complete.\n")

# ============================================================
# OUTPUT COLUMN CHECK
# ============================================================

required_test_columns <- c(
    "we.ep",
    "we.eBH",
    "wi.ep",
    "wi.eBH"
)

if (!all(
    required_test_columns %in%
    colnames(aldex_test)
)) {

    stop(
        paste(
            "Unexpected aldex.ttest columns:",
            paste(
                colnames(aldex_test),
                collapse = ", "
            )
        )
    )
}

required_effect_columns <- c(
    "rab.win.Control",
    "rab.win.Case",
    "diff.btw",
    "diff.win",
    "effect",
    "overlap"
)

if (!all(
    required_effect_columns %in%
    colnames(aldex_effect)
)) {

    stop(
        paste(
            "Unexpected aldex.effect columns:",
            paste(
                colnames(aldex_effect),
                collapse = ", "
            )
        )
    )
}

# ============================================================
# ALIGN RESULTS
# ============================================================

canonical_taxa <- rownames(
    counts_matrix
)

test_index <- match(
    canonical_taxa,
    rownames(aldex_test)
)

effect_index <- match(
    canonical_taxa,
    rownames(aldex_effect)
)

taxonomy_index <- match(
    canonical_taxa,
    taxonomy$Taxon
)

if (any(is.na(test_index))) {
    stop("Taxa missing from aldex.ttest output.")
}

if (any(is.na(effect_index))) {
    stop("Taxa missing from aldex.effect output.")
}

if (any(is.na(taxonomy_index))) {
    stop("Taxa missing from taxonomy mapping.")
}

aldex_test <- aldex_test[
    test_index,
    ,
    drop = FALSE
]

aldex_effect <- aldex_effect[
    effect_index,
    ,
    drop = FALSE
]

taxonomy_aligned <- taxonomy[
    taxonomy_index,
    ,
    drop = FALSE
]

# ============================================================
# EFFECT DIRECTION
#
# Explicitly calculate:
# Case median CLR - Control median CLR
# ============================================================

control_median_clr <- (
    aldex_effect$rab.win.Control
)

case_median_clr <- (
    aldex_effect$rab.win.Case
)

case_minus_control <- (
    case_median_clr -
    control_median_clr
)

effect_direction <- ifelse(
    case_minus_control > 0,
    "Case_higher",
    ifelse(
        case_minus_control < 0,
        "Control_higher",
        "No_difference"
    )
)

# ============================================================
# SIGNIFICANCE
# ============================================================

wilcoxon_significant <- (
    !is.na(aldex_test$wi.eBH) &
    aldex_test$wi.eBH < ALPHA
)

welch_significant <- (
    !is.na(aldex_test$we.eBH) &
    aldex_test$we.eBH < ALPHA
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
# STANDARDIZED RESULT
# ============================================================

results <- data.frame(

    Taxon = canonical_taxa,

    Method = "ALDEx2",

    Tested = !is.na(
        aldex_test$wi.ep
    ),

    P_Value = aldex_test$wi.ep,

    Adjusted_P_Value = aldex_test$wi.eBH,

    Significant = wilcoxon_significant,

    Effect_Estimate = case_minus_control,

    Effect_Metric =
        "ALDEx2 median CLR difference (Case-Control)",

    Effect_Direction = effect_direction,

    Control_Median_CLR =
        control_median_clr,

    Case_Median_CLR =
        case_median_clr,

    Wilcoxon_P_Value =
        aldex_test$wi.ep,

    Wilcoxon_Adjusted_P_Value =
        aldex_test$wi.eBH,

    Wilcoxon_Significant =
        wilcoxon_significant,

    Welch_P_Value =
        aldex_test$we.ep,

    Welch_Adjusted_P_Value =
        aldex_test$we.eBH,

    Welch_Significant =
        welch_significant,

    ALDEx2_Native_Diff_Btw =
        aldex_effect$diff.btw,

    ALDEx2_Native_Diff_Win =
        aldex_effect$diff.win,

    ALDEx2_Native_Effect =
        aldex_effect$effect,

    ALDEx2_Overlap =
        aldex_effect$overlap,

    Present_Samples =
        as.integer(present_total),

    Present_Control =
        as.integer(present_control),

    Present_Case =
        as.integer(present_case),

    Alpha = ALPHA,

    MC_Samples = MC_SAMPLES,

    Denominator = "all",

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

# ============================================================
# SUBTABLES
# ============================================================

wilcoxon_results <- results[
    ,
    c(
        "Taxon",
        "Method",
        "Tested",
        "Wilcoxon_P_Value",
        "Wilcoxon_Adjusted_P_Value",
        "Wilcoxon_Significant",
        "Effect_Estimate",
        "Effect_Direction",
        "Present_Samples",
        "Present_Control",
        "Present_Case",
        "Denovo_ID",
        "Phylum",
        "Family",
        "Genus",
        "Source_Taxon"
    ),
    drop = FALSE
]

welch_results <- results[
    ,
    c(
        "Taxon",
        "Method",
        "Welch_P_Value",
        "Welch_Adjusted_P_Value",
        "Welch_Significant",
        "Effect_Estimate",
        "Effect_Direction",
        "Present_Samples",
        "Present_Control",
        "Present_Case",
        "Denovo_ID",
        "Phylum",
        "Family",
        "Genus",
        "Source_Taxon"
    ),
    drop = FALSE
]

significant_wilcoxon <- results[
    results$Wilcoxon_Significant,
    ,
    drop = FALSE
]

# ============================================================
# SUMMARY VALUES
# ============================================================

n_available <- nrow(results)

n_tested_wilcoxon <- sum(
    results$Tested,
    na.rm = TRUE
)

n_sig_wilcoxon <- sum(
    results$Wilcoxon_Significant,
    na.rm = TRUE
)

n_sig_welch <- sum(
    results$Welch_Significant,
    na.rm = TRUE
)

case_higher_sig <- sum(
    results$Wilcoxon_Significant &
    results$Effect_Direction == "Case_higher",
    na.rm = TRUE
)

control_higher_sig <- sum(
    results$Wilcoxon_Significant &
    results$Effect_Direction == "Control_higher",
    na.rm = TRUE
)

min_wilcoxon_p <- min(
    results$Wilcoxon_P_Value,
    na.rm = TRUE
)

min_wilcoxon_q <- min(
    results$Wilcoxon_Adjusted_P_Value,
    na.rm = TRUE
)

min_welch_p <- min(
    results$Welch_P_Value,
    na.rm = TRUE
)

min_welch_q <- min(
    results$Welch_Adjusted_P_Value,
    na.rm = TRUE
)

# ============================================================
# PRINT RESULT SUMMARY
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

cat("\nWILCOXON PRIMARY RESULT\n")

cat(
    "Taxa tested: ",
    n_tested_wilcoxon,
    "\n",
    sep = ""
)

cat(
    "Significant taxa at BH q < 0.05: ",
    n_sig_wilcoxon,
    "\n",
    sep = ""
)

cat(
    "Case-higher significant taxa: ",
    case_higher_sig,
    "\n",
    sep = ""
)

cat(
    "Control-higher significant taxa: ",
    control_higher_sig,
    "\n",
    sep = ""
)

cat(
    "Minimum raw p-value: ",
    format(
        min_wilcoxon_p,
        digits = 10
    ),
    "\n",
    sep = ""
)

cat(
    "Minimum adjusted p-value: ",
    format(
        min_wilcoxon_q,
        digits = 10
    ),
    "\n",
    sep = ""
)

cat("\nWELCH DIAGNOSTIC RESULT\n")

cat(
    "Significant taxa at BH q < 0.05: ",
    n_sig_welch,
    "\n",
    sep = ""
)

cat(
    "Minimum raw p-value: ",
    format(
        min_welch_p,
        digits = 10
    ),
    "\n",
    sep = ""
)

cat(
    "Minimum adjusted p-value: ",
    format(
        min_welch_q,
        digits = 10
    ),
    "\n",
    sep = ""
)

# ============================================================
# TOP 20
# ============================================================

cat(
    "\nTop 20 taxa by Wilcoxon adjusted p-value:\n"
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
                "Wilcoxon_P_Value",
                "Wilcoxon_Adjusted_P_Value",
                "Wilcoxon_Significant",
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
# SAVE
# ============================================================

write.csv(
    results,
    RESULT_FILE,
    row.names = FALSE,
    na = ""
)

write.csv(
    wilcoxon_results,
    WILCOXON_FILE,
    row.names = FALSE,
    na = ""
)

write.csv(
    welch_results,
    WELCH_FILE,
    row.names = FALSE,
    na = ""
)

write.csv(
    significant_wilcoxon,
    SIGNIFICANT_FILE,
    row.names = FALSE,
    na = ""
)

# ============================================================
# SUMMARY FILE
# ============================================================

summary_lines <- c(

    "PAR_SCHEPERJANS ALDEX2 SUMMARY",

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
        "ALDEx2 version: ",
        as.character(
            packageVersion("ALDEx2")
        )
    ),

    paste0(
        "Monte Carlo samples: ",
        MC_SAMPLES
    ),

    "Denominator: all",

    "Conditions supplied as character vector",

    "Paired test: FALSE",

    paste0(
        "Seed: ",
        SEED
    ),

    "",

    "PRIMARY WILCOXON RESULT",

    paste0(
        "Taxa tested: ",
        n_tested_wilcoxon
    ),

    paste0(
        "Significant taxa: ",
        n_sig_wilcoxon
    ),

    paste0(
        "Case-higher significant taxa: ",
        case_higher_sig
    ),

    paste0(
        "Control-higher significant taxa: ",
        control_higher_sig
    ),

    paste0(
        "Minimum raw p-value: ",
        format(
            min_wilcoxon_p,
            digits = 10
        )
    ),

    paste0(
        "Minimum adjusted p-value: ",
        format(
            min_wilcoxon_q,
            digits = 10
        )
    ),

    "",

    "WELCH DIAGNOSTIC RESULT",

    paste0(
        "Significant taxa: ",
        n_sig_welch
    ),

    paste0(
        "Minimum raw p-value: ",
        format(
            min_welch_p,
            digits = 10
        )
    ),

    paste0(
        "Minimum adjusted p-value: ",
        format(
            min_welch_q,
            digits = 10
        )
    ),

    "",

    "EFFECT DIRECTION",

    paste0(
        "Standardized effect = ",
        "ALDEx2 median CLR(Case) - ",
        "median CLR(Control)."
    ),

    "",

    "INTERPRETATION",

    paste0(
        "Real-data discoveries have unknown ",
        "ground truth and are not classified ",
        "as true positives."
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
    WILCOXON_FILE,
    "\n"
)

cat(
    WELCH_FILE,
    "\n"
)

cat(
    SIGNIFICANT_FILE,
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
    "ALDEX2 REAL-DATA ANALYSIS: PASS\n"
)

cat(
    strrep("=", 72),
    "\n",
    sep = ""
)
