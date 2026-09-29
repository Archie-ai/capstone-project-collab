#!/usr/bin/env Rscript

# ============================================================
# ANCOM-BC 
# ============================================================
#
# This script runs ANCOM-BC on microbiome count data that passed shared QC.
#
# This script was developed with assistance from OpenAI's ChatGPT
#

suppressPackageStartupMessages({
    library(ANCOMBC)
    library(phyloseq)
})


# ============================================================
# COMMAND-LINE INPUTS
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

# Expected arguments:
#
# 1. QC count table
# 2. QC metadata
# 3. output directory
# 4. dataset label
#
# Example:
#
# Rscript ancombc.R \
#   erinmaley/data/qc_2/qc_null_counts.csv \
#   erinmaley/data/qc_2/qc_null_metadata.csv \
#   erinmaley/data/method_outputs_3/ancombc \
#   null

if (length(args) != 4) {
    stop(
        paste(
            "Usage:",
            "Rscript ancombc.R",
            "<counts.csv>",
            "<metadata.csv>",
            "<output_directory>",
            "<dataset_label>"
        )
    )
}

count_file <- args[1]
metadata_file <- args[2]
output_dir <- args[3]
dataset_label <- args[4]


# ============================================================
# CHECK INPUT FILES
# ============================================================

if (!file.exists(count_file)) {
    stop(paste("Count file not found:", count_file))
}

if (!file.exists(metadata_file)) {
    stop(paste("Metadata file not found:", metadata_file))
}

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# LOAD QC OUTPUTS
# ============================================================

counts <- read.csv(
    count_file,
    row.names = 1,
    check.names = FALSE
)

metadata <- read.csv(
    metadata_file,
    check.names = FALSE
)

cat("\nRunning ANCOM-BC for dataset:", dataset_label, "\n\n")

cat("Input taxa:", nrow(counts), "\n")
cat("Input samples:", ncol(counts), "\n")


# ============================================================
# CHECK REQUIRED METADATA
# ============================================================

required_columns <- c("SampleID", "Group")

missing_columns <- setdiff(
    required_columns,
    colnames(metadata)
)

if (length(missing_columns) > 0) {
    stop(
        paste(
            "Metadata is missing required column(s):",
            paste(missing_columns, collapse = ", ")
        )
    )
}


# ============================================================
# MATCH SAMPLES
# ============================================================

rownames(metadata) <- metadata$SampleID

if (!setequal(colnames(counts), rownames(metadata))) {
    stop(
        "Sample IDs do not match between the QC count table and metadata."
    )
}

# Put metadata in exactly the same order as the count table.
metadata <- metadata[
    colnames(counts),
    ,
    drop = FALSE
]


# ============================================================
# DEFINE CASE / CONTROL COMPARISON
# ============================================================

# Control is intentionally set as the reference group.
#
# Therefore:
#   positive log fold change = higher abundance in Case
#   negative log fold change = lower abundance in Case

metadata$Group <- factor(
    metadata$Group,
    levels = c("Control", "Case")
)

if (any(is.na(metadata$Group))) {
    stop(
        "Group must contain the expected labels: Control and Case."
    )
}

cat(
    "Groups:",
    paste(
        names(table(metadata$Group)),
        table(metadata$Group),
        collapse = ", "
    ),
    "\n\n"
)


# ============================================================
# BUILD PHYLOSEQ OBJECT
# ============================================================

# ANCOM-BC accepts a phyloseq object containing:
# - the taxon count table
# - the sample metadata

count_matrix <- as.matrix(counts)

storage.mode(count_matrix) <- "numeric"

otu <- otu_table(
    count_matrix,
    taxa_are_rows = TRUE
)

sample_info <- sample_data(
    metadata
)

ps <- phyloseq(
    otu,
    sample_info
)


# ============================================================
# RUN ANCOM-BC
# ============================================================

set.seed(42)

ancombc_result <- ancombc(
    data = ps,

    # Group is the only predictor in this test.
    formula = "Group",

    # Use BH/FDR adjustment so adjusted p-values can be
    # compared consistently with the other DA methods later.
    p_adj_method = "BH",

    # Shared filtering already occurred during QC stage 2.
    # Do not apply another prevalence filter here.
    prv_cut = 0,

    # QC already checked library sizes.
    lib_cut = 0,

    # Two-group comparison; structural-zero testing is not
    # required for this benchmark run.
    group = NULL,
    struc_zero = FALSE,

    neg_lb = FALSE,
    tol = 1e-5,
    max_iter = 100,
    conserve = FALSE,

    alpha = 0.05,

    # No multi-group global test is needed for Control vs Case.
    global = FALSE,

    n_cl = 1,
    verbose = TRUE
)


# ============================================================
# EXTRACT CASE VS CONTROL RESULT
# ============================================================

result <- ancombc_result$res

# ANCOM-BC includes the taxon ID as a column in the result table.
# Exclude it so only Case vs Control model coefficient remains.

coefficient_names <- setdiff(
    colnames(result$lfc),
    c("taxon", "(Intercept)")
)

if (length(coefficient_names) != 1) {
    stop(
        paste(
            "Expected one Group comparison coefficient but found:",
            paste(coefficient_names, collapse = ", ")
        )
    )
}

group_coefficient <- coefficient_names[1]

cat(
    "\nANCOM-BC comparison coefficient:",
    group_coefficient,
    "\n"
)


# ============================================================
# CREATE METHOD OUTPUT TABLE
# ============================================================

ancombc_table <- data.frame(

    Taxon = result$lfc$taxon,

    LogFoldChange =
        result$lfc[, group_coefficient],

    StandardError =
        result$se[, group_coefficient],

    Statistic =
        result$W[, group_coefficient],

    PValue =
        result$p_val[, group_coefficient],

    AdjustedPValue =
        result$q_val[, group_coefficient],

    Significant =
        result$diff_abn[, group_coefficient],

    stringsAsFactors = FALSE
)


# Add direction for easier interpretation later.

ancombc_table$Direction <- ifelse(
    ancombc_table$LogFoldChange > 0,
    "Higher_in_Case",
    ifelse(
        ancombc_table$LogFoldChange < 0,
        "Higher_in_Control",
        "No_change"
    )
)


# ============================================================
# SAVE METHOD OUTPUT
# ============================================================

output_file <- file.path(
    output_dir,
    paste0(
        "ancombc_",
        dataset_label,
        "_results.csv"
    )
)

write.csv(
    ancombc_table,
    output_file,
    row.names = FALSE
)


# ============================================================
# SAVE RUN INFORMATION
# ============================================================

summary_file <- file.path(
    output_dir,
    paste0(
        "ancombc_",
        dataset_label,
        "_summary.txt"
    )
)

significant_count <- sum(
    ancombc_table$Significant,
    na.rm = TRUE
)

summary_text <- c(
    "ANCOM-BC METHOD SUMMARY",
    "=======================",
    "",
    paste("Dataset:", dataset_label),
    paste("ANCOMBC version:", packageVersion("ANCOMBC")),
    paste("Samples:", ncol(counts)),
    paste("Taxa analyzed:", nrow(ancombc_table)),
    "",
    "Comparison:",
    "Case vs Control",
    "",
    "Reference group:",
    "Control",
    "",
    "P-value adjustment:",
    "Benjamini-Hochberg (BH)",
    "",
    "Significance threshold:",
    "Adjusted p-value < 0.05",
    "",
    paste(
        "Significant taxa:",
        significant_count
    ),
    "",
    "Direction:",
    "Positive log fold change = higher in Case",
    "Negative log fold change = higher in Control",
    "",
    "QC filtering was completed before ANCOM-BC.",
    "No additional prevalence or library-size filtering was applied."
)

writeLines(
    summary_text,
    summary_file
)


# ============================================================
# TERMINAL OUTPUT
# ============================================================

cat("\nANCOM-BC COMPLETED\n")
cat("==================\n")

cat(
    "Taxa analyzed:",
    nrow(ancombc_table),
    "\n"
)

cat(
    "Significant taxa:",
    significant_count,
    "\n"
)

cat("\nOutput:\n")
cat(" ", output_file, "\n")

cat("\nSummary:\n")
cat(" ", summary_file, "\n\n")
