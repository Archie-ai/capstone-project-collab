#!/usr/bin/env Rscript

suppressPackageStartupMessages(
    library(DESeq2)
)

# ============================================================
# DESeq2 REAL-DATA ANALYSIS
# par_scheperjans
#
# Microbiome Abundance Consensus Machine
#
# Canonical input:
# 691 taxa x 148 samples
# 74 Control / 74 Case
#
# Size-factor estimation:
# poscounts
#
# Primary result:
# DESeq2 Wald test
# Case vs Control
# independent filtering ON
# Cook's filtering ON
#
# Diagnostics:
# 1. independent filtering OFF
# 2. independent filtering OFF + Cook's filtering OFF
# ============================================================


# ============================================================
# CONFIGURATION
# ============================================================

DATASET <- "par_scheperjans"

ALPHA <- 0.05

COUNTS_FILE <- "data/processed/par_scheperjans/par_scheperjans_counts.csv"

METADATA_FILE <- "data/processed/par_scheperjans/par_scheperjans_metadata.csv"

TAXONOMY_FILE <- "data/processed/par_scheperjans/par_scheperjans_taxonomy_mapping.csv"

OUTPUT_DIR <- "results/real/par_scheperjans/deseq2"


dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


RESULT_FILE <- file.path(
    OUTPUT_DIR,
    "deseq2_results.csv"
)

SIGNIFICANT_FILE <- file.path(
    OUTPUT_DIR,
    "deseq2_significant_taxa.csv"
)

NO_IF_FILE <- file.path(
    OUTPUT_DIR,
    "deseq2_no_independent_filter.csv"
)

FULL_DIAGNOSTIC_FILE <- file.path(
    OUTPUT_DIR,
    "deseq2_no_filter_no_cooks.csv"
)

SIZE_FACTOR_FILE <- file.path(
    OUTPUT_DIR,
    "deseq2_size_factors.csv"
)

SUMMARY_FILE <- file.path(
    OUTPUT_DIR,
    "deseq2_summary.txt"
)


# ============================================================
# HEADER
# ============================================================

cat(
    strrep("=", 72),
    "\n",
    sep = ""
)

cat(
    "PAR_SCHEPERJANS DESEQ2 ANALYSIS\n"
)

cat(
    strrep("=", 72),
    "\n",
    sep = ""
)


cat(
    "\nR version: ",
    as.character(getRversion()),
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


cat(
    "\nCanonical input dimensions:\n"
)

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
# REQUIRED COLUMNS
# ============================================================

if (!all(
    c(
        "SampleID",
        "Group"
    ) %in%
    colnames(metadata)
)) {

    stop(
        "Metadata must contain SampleID and Group."
    )
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


if (!all(
    required_taxonomy %in%
    colnames(taxonomy)
)) {

    stop(
        "Taxonomy mapping is missing required columns."
    )
}


# ============================================================
# CLEAN IDENTIFIERS
# ============================================================

metadata$SampleID <- trimws(
    as.character(
        metadata$SampleID
    )
)

metadata$Group <- trimws(
    as.character(
        metadata$Group
    )
)

colnames(counts) <- trimws(
    colnames(counts)
)

rownames(counts) <- trimws(
    rownames(counts)
)

taxonomy$Taxon <- trimws(
    as.character(
        taxonomy$Taxon
    )
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
# SAMPLE ALIGNMENT
# ============================================================

if (!setequal(
    colnames(counts),
    metadata$SampleID
)) {

    stop(
        "Count-table and metadata sample sets do not match."
    )
}


metadata_index <- match(
    colnames(counts),
    metadata$SampleID
)


if (any(is.na(metadata_index))) {

    stop(
        "Metadata alignment failed."
    )
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


cat(
    "\nSample alignment:\n"
)

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
# GROUP CHECKS
# ============================================================

n_control <- sum(
    metadata$Group == "Control"
)

n_case <- sum(
    metadata$Group == "Case"
)


cat(
    "\nGroup counts:\n"
)

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
    c(
        "Control",
        "Case"
    )
)) {

    stop(
        "Unexpected group labels detected."
    )
}


# ============================================================
# COUNT MATRIX QC
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


cat(
    "\nCount QC:\n"
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
    sprintf(
        "%.6f",
        zero_fraction
    ),
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

storage.mode(
    counts_matrix
) <- "integer"


# ============================================================
# DESEQ2 METADATA
#
# Control is explicit reference.
# ============================================================

col_data <- data.frame(

    Group = factor(
        metadata$Group,
        levels = c(
            "Control",
            "Case"
        )
    ),

    row.names = metadata$SampleID,

    stringsAsFactors = FALSE
)


if (!identical(
    colnames(counts_matrix),
    rownames(col_data)
)) {

    stop(
        "DESeq2 metadata alignment failed."
    )
}


# ============================================================
# CREATE DESEQ DATASET
# ============================================================

dds <- DESeqDataSetFromMatrix(

    countData = counts_matrix,

    colData = col_data,

    design = ~ Group
)


# ============================================================
# SIZE FACTORS
#
# poscounts handles sparse count data with many zeros.
# ============================================================

cat(
    "\nEstimating DESeq2 size factors using poscounts...\n"
)


dds <- estimateSizeFactors(
    dds,
    type = "poscounts"
)


sf <- sizeFactors(
    dds
)


if (any(is.na(sf))) {
    stop("Missing DESeq2 size factors detected.")
}

if (any(!is.finite(sf))) {
    stop("Non-finite DESeq2 size factors detected.")
}

if (any(sf <= 0)) {
    stop("Non-positive DESeq2 size factors detected.")
}


cat(
    "Minimum size factor: ",
    format(
        min(sf),
        digits = 8
    ),
    "\n",
    sep = ""
)

cat(
    "Median size factor: ",
    format(
        median(sf),
        digits = 8
    ),
    "\n",
    sep = ""
)

cat(
    "Maximum size factor: ",
    format(
        max(sf),
        digits = 8
    ),
    "\n",
    sep = ""
)


# ============================================================
# RUN DESEQ2
# ============================================================

cat(
    "\nRunning DESeq2 Wald model...\n"
)


dds <- DESeq(
    dds,
    test = "Wald",
    fitType = "parametric",
    quiet = FALSE
)


cat(
    "DESeq2 model complete.\n"
)


# ============================================================
# ACTUAL DISPERSION FIT
# ============================================================

fit_type_used <- attr(
    dispersionFunction(dds),
    "fitType"
)


if (is.null(fit_type_used)) {

    fit_type_used <- "unknown"
}


cat(
    "\nDispersion fit used: ",
    fit_type_used,
    "\n",
    sep = ""
)


# ============================================================
# PRIMARY RESULTS
#
# Default:
# independent filtering ON
# Cook's filtering ON
# ============================================================

res_primary <- results(

    dds,

    contrast = c(
        "Group",
        "Case",
        "Control"
    ),

    alpha = ALPHA
)


# ============================================================
# DIAGNOSTIC 1
#
# Independent filtering OFF
# ============================================================

res_no_if <- results(

    dds,

    contrast = c(
        "Group",
        "Case",
        "Control"
    ),

    alpha = ALPHA,

    independentFiltering = FALSE
)


# ============================================================
# DIAGNOSTIC 2
#
# Independent filtering OFF
# Cook's filtering OFF
# ============================================================

res_full <- results(

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


# ============================================================
# ALIGN RESULTS
# ============================================================

canonical_taxa <- rownames(
    counts_matrix
)


primary_df <- as.data.frame(
    res_primary
)

no_if_df <- as.data.frame(
    res_no_if
)

full_df <- as.data.frame(
    res_full
)


primary_index <- match(
    canonical_taxa,
    rownames(primary_df)
)

no_if_index <- match(
    canonical_taxa,
    rownames(no_if_df)
)

full_index <- match(
    canonical_taxa,
    rownames(full_df)
)

taxonomy_index <- match(
    canonical_taxa,
    taxonomy$Taxon
)


if (any(is.na(primary_index))) {
    stop("Taxa missing from primary DESeq2 result.")
}

if (any(is.na(no_if_index))) {
    stop("Taxa missing from no-filter DESeq2 result.")
}

if (any(is.na(full_index))) {
    stop("Taxa missing from full diagnostic result.")
}

if (any(is.na(taxonomy_index))) {
    stop("Taxa missing from taxonomy mapping.")
}


primary_df <- primary_df[
    primary_index,
    ,
    drop = FALSE
]

no_if_df <- no_if_df[
    no_if_index,
    ,
    drop = FALSE
]

full_df <- full_df[
    full_index,
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
# ============================================================

primary_significant <- (
    !is.na(primary_df$padj) &
    primary_df$padj < ALPHA
)


no_if_significant <- (
    !is.na(no_if_df$padj) &
    no_if_df$padj < ALPHA
)


full_significant <- (
    !is.na(full_df$padj) &
    full_df$padj < ALPHA
)


effect_direction <- ifelse(
    primary_df$log2FoldChange > 0,
    "Case_higher",
    ifelse(
        primary_df$log2FoldChange < 0,
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
# BUILD PRIMARY STANDARDIZED TABLE
# ============================================================

deseq_results <- data.frame(

    Taxon = canonical_taxa,

    Method = "DESeq2",

    Tested = !is.na(
        primary_df$pvalue
    ),

    P_Value =
        primary_df$pvalue,

    Adjusted_P_Value =
        primary_df$padj,

    Significant =
        primary_significant,

    Effect_Estimate =
        primary_df$log2FoldChange,

    Effect_Metric =
        "DESeq2 log2 fold change (Case vs Control)",

    Effect_Direction =
        effect_direction,

    Base_Mean =
        primary_df$baseMean,

    Log2_Fold_Change =
        primary_df$log2FoldChange,

    LFC_Standard_Error =
        primary_df$lfcSE,

    Wald_Statistic =
        primary_df$stat,

    NoIF_P_Value =
        no_if_df$pvalue,

    NoIF_Adjusted_P_Value =
        no_if_df$padj,

    NoIF_Significant =
        no_if_significant,

    FullDiagnostic_P_Value =
        full_df$pvalue,

    FullDiagnostic_Adjusted_P_Value =
        full_df$padj,

    FullDiagnostic_Significant =
        full_significant,

    Present_Samples =
        as.integer(present_total),

    Present_Control =
        as.integer(present_control),

    Present_Case =
        as.integer(present_case),

    Alpha = ALPHA,

    Size_Factor_Method = "poscounts",

    Dispersion_Fit =
        as.character(fit_type_used),

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
# BUILD DIAGNOSTIC TABLES
# ============================================================

no_if_results <- data.frame(

    Taxon = canonical_taxa,

    P_Value =
        no_if_df$pvalue,

    Adjusted_P_Value =
        no_if_df$padj,

    Significant =
        no_if_significant,

    Log2_Fold_Change =
        no_if_df$log2FoldChange,

    stringsAsFactors = FALSE
)


full_results <- data.frame(

    Taxon = canonical_taxa,

    P_Value =
        full_df$pvalue,

    Adjusted_P_Value =
        full_df$padj,

    Significant =
        full_significant,

    Log2_Fold_Change =
        full_df$log2FoldChange,

    stringsAsFactors = FALSE
)


# ============================================================
# SORT PRIMARY RESULTS
# ============================================================

sort_q <- deseq_results$Adjusted_P_Value

sort_p <- deseq_results$P_Value


sort_q[
    is.na(sort_q)
] <- Inf

sort_p[
    is.na(sort_p)
] <- Inf


sort_index <- order(
    sort_q,
    sort_p,
    deseq_results$Taxon
)


deseq_results <- deseq_results[
    sort_index,
    ,
    drop = FALSE
]


row.names(
    deseq_results
) <- NULL


significant_results <- deseq_results[
    deseq_results$Significant,
    ,
    drop = FALSE
]


# ============================================================
# SIZE FACTOR TABLE
# ============================================================

size_factor_table <- data.frame(

    SampleID = names(sf),

    Size_Factor = as.numeric(sf),

    Group = metadata$Group,

    stringsAsFactors = FALSE
)


# ============================================================
# SUMMARY VALUES
# ============================================================

n_available <- nrow(
    deseq_results
)


n_raw_p <- sum(
    !is.na(
        deseq_results$P_Value
    )
)


n_primary_q <- sum(
    !is.na(
        deseq_results$Adjusted_P_Value
    )
)


n_primary_sig <- sum(
    deseq_results$Significant,
    na.rm = TRUE
)


n_no_if_sig <- sum(
    no_if_results$Significant,
    na.rm = TRUE
)


n_full_sig <- sum(
    full_results$Significant,
    na.rm = TRUE
)


n_primary_na_q <- sum(
    is.na(
        deseq_results$Adjusted_P_Value
    )
)


case_higher_sig <- sum(
    deseq_results$Significant &
    deseq_results$Effect_Direction ==
        "Case_higher",
    na.rm = TRUE
)


control_higher_sig <- sum(
    deseq_results$Significant &
    deseq_results$Effect_Direction ==
        "Control_higher",
    na.rm = TRUE
)


min_p <- min(
    deseq_results$P_Value,
    na.rm = TRUE
)


min_q <- min(
    deseq_results$Adjusted_P_Value,
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

cat(
    "RESULT SUMMARY\n"
)

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
    "Taxa with raw p-values: ",
    n_raw_p,
    "\n",
    sep = ""
)

cat(
    "Taxa with primary adjusted p-values: ",
    n_primary_q,
    "\n",
    sep = ""
)

cat(
    "Taxa with NA primary adjusted p-values: ",
    n_primary_na_q,
    "\n",
    sep = ""
)

cat(
    "Primary significant taxa: ",
    n_primary_sig,
    "\n",
    sep = ""
)

cat(
    "Primary Case-higher significant taxa: ",
    case_higher_sig,
    "\n",
    sep = ""
)

cat(
    "Primary Control-higher significant taxa: ",
    control_higher_sig,
    "\n",
    sep = ""
)

cat(
    "No-independent-filter significant taxa: ",
    n_no_if_sig,
    "\n",
    sep = ""
)

cat(
    "No-filter/no-Cook's significant taxa: ",
    n_full_sig,
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
    "\nTop 20 taxa by primary adjusted p-value:\n"
)


print(
    head(
        deseq_results[
            ,
            c(
                "Taxon",
                "Denovo_ID",
                "Present_Samples",
                "Present_Control",
                "Present_Case",
                "Base_Mean",
                "Effect_Estimate",
                "Effect_Direction",
                "P_Value",
                "Adjusted_P_Value",
                "Significant",
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
    deseq_results,
    RESULT_FILE,
    row.names = FALSE,
    na = ""
)


write.csv(
    significant_results,
    SIGNIFICANT_FILE,
    row.names = FALSE,
    na = ""
)


write.csv(
    no_if_results,
    NO_IF_FILE,
    row.names = FALSE,
    na = ""
)


write.csv(
    full_results,
    FULL_DIAGNOSTIC_FILE,
    row.names = FALSE,
    na = ""
)


write.csv(
    size_factor_table,
    SIZE_FACTOR_FILE,
    row.names = FALSE,
    na = ""
)


# ============================================================
# SUMMARY TEXT
# ============================================================

summary_lines <- c(

    "PAR_SCHEPERJANS DESEQ2 SUMMARY",

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
        "DESeq2 version: ",
        as.character(
            packageVersion("DESeq2")
        )
    ),

    "",

    "SETTINGS",

    "Design: ~ Group",

    "Reference group: Control",

    "Contrast: Case vs Control",

    "Test: Wald",

    "Size-factor method: poscounts",

    paste0(
        "Dispersion fit used: ",
        fit_type_used
    ),

    paste0(
        "Alpha: ",
        ALPHA
    ),

    "",

    "PRIMARY RESULTS",

    paste0(
        "Taxa with raw p-values: ",
        n_raw_p
    ),

    paste0(
        "Taxa with adjusted p-values: ",
        n_primary_q
    ),

    paste0(
        "Taxa with NA adjusted p-values: ",
        n_primary_na_q
    ),

    paste0(
        "Significant taxa: ",
        n_primary_sig
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

    "DIAGNOSTICS",

    paste0(
        "No-independent-filter significant taxa: ",
        n_no_if_sig
    ),

    paste0(
        "No-independent-filter/no-Cook's significant taxa: ",
        n_full_sig
    ),

    "",

    "EFFECT",

    paste0(
        "Effect estimate is DESeq2 log2 fold change ",
        "for Case versus Control."
    ),

    "",

    "INTERPRETATION",

    paste0(
        "Real-data discoveries have unknown ground truth ",
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

cat(
    "\nOutput files:\n"
)

cat(
    RESULT_FILE,
    "\n"
)

cat(
    SIGNIFICANT_FILE,
    "\n"
)

cat(
    NO_IF_FILE,
    "\n"
)

cat(
    FULL_DIAGNOSTIC_FILE,
    "\n"
)

cat(
    SIZE_FACTOR_FILE,
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
    "DESEQ2 REAL-DATA ANALYSIS: PASS\n"
)

cat(
    strrep("=", 72),
    "\n",
    sep = ""
)
