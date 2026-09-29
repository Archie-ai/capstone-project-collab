#!/usr/bin/env Rscript

# ============================================================
# ALDEx2
# ============================================================
#
# This script runs ALDEx2 on microbiome count data that passed shared QC.
#
# This script was developed with assistance from OpenAI's ChatGPT
#

suppressPackageStartupMessages({
    library(ALDEx2)
})


# ============================================================
# COMMAND-LINE INPUTS
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

# Expected arguments:
# 1. QC count table
# 2. QC metadata
# 3. output directory
# 4. dataset label

if (length(args) != 4) {
    stop(
        paste(
            "Usage:",
            "Rscript aldex2.R",
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

cat("\nRunning ALDEx2 for dataset:", dataset_label, "\n\n")

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

# Put metadata in the same order as the count-table columns.
metadata <- metadata[
    colnames(counts),
    ,
    drop = FALSE
]


# ============================================================
# DEFINE CASE / CONTROL COMPARISON
# ============================================================

# Keep the same direction convention used by the other methods.
#
# Control = reference
# Positive difference/effect = higher in Case
# Negative difference/effect = higher in Control

metadata$Group <- factor(
    metadata$Group,
    levels = c("Control", "Case")
)

if (any(is.na(metadata$Group))) {
    stop(
        "Group must contain the expected labels: Control and Case."
    )
}

conditions <- as.character(metadata$Group)

cat(
    "Groups:",
    paste(
        names(table(conditions)),
        table(conditions),
        collapse = ", "
    ),
    "\n\n"
)


# ============================================================
# PREPARE COUNT MATRIX
# ============================================================

count_matrix <- as.matrix(counts)

storage.mode(count_matrix) <- "integer"


# ============================================================
# RUN ALDEx2 CLR MONTE CARLO MODEL
# ============================================================

# ALDEx2 uses Monte Carlo sampling to represent uncertainty
# in the compositional count data before statistical testing.

set.seed(42)

clr_result <- aldex.clr(
    reads = count_matrix,
    conds = conditions,
    mc.samples = 128,
    denom = "all",
    verbose = FALSE
)


# ============================================================
# RUN ALDEx2 STATISTICAL TESTS
# ============================================================

# aldex.ttest reports both:
# - Welch's t-test
# - Wilcoxon rank-sum test
#
# Both are preserved in the output.

test_result <- aldex.ttest(
    clr_result,
    paired.test = FALSE,
    verbose = FALSE
)


# ============================================================
# CALCULATE ALDEx2 EFFECT SIZES
# ============================================================

effect_result <- aldex.effect(
    clr_result,
    verbose = FALSE
)


# ============================================================
# CREATE METHOD OUTPUT TABLE
# ============================================================

aldex2_table <- data.frame(

    Taxon = rownames(test_result),

    # Expected Welch test results.
    WelchPValue =
        test_result$we.ep,

    WelchAdjustedPValue =
        test_result$we.eBH,

    # Expected Wilcoxon test results.
    WilcoxonPValue =
        test_result$wi.ep,

    WilcoxonAdjustedPValue =
        test_result$wi.eBH,

    # ALDEx2 effect information.
    DifferenceBetween =
        effect_result$diff.btw,

    DifferenceWithin =
        effect_result$diff.win,

    Effect =
        effect_result$effect,

    Overlap =
        effect_result$overlap,

    stringsAsFactors = FALSE
)


# ============================================================
# PRIMARY RESULT FIELDS
# ============================================================

# Use ALDEx2's Wilcoxon expected p-value as the primary
# p-value for the later standardized-results step.

aldex2_table$PValue <- (
    aldex2_table$WilcoxonPValue
)

aldex2_table$AdjustedPValue <- (
    aldex2_table$WilcoxonAdjustedPValue
)

aldex2_table$Significant <- (
    !is.na(aldex2_table$AdjustedPValue)
    & aldex2_table$AdjustedPValue < 0.05
)


# ============================================================
# ADD DIRECTION
# ============================================================

# diff.btw is the median difference between the two groups.
# With Control first and Case second:
#
# positive = higher in Case
# negative = higher in Control

aldex2_table$Direction <- ifelse(
    is.na(aldex2_table$DifferenceBetween),
    "Not_available",
    ifelse(
        aldex2_table$DifferenceBetween > 0,
        "Higher_in_Case",
        ifelse(
            aldex2_table$DifferenceBetween < 0,
            "Higher_in_Control",
            "No_change"
        )
    )
)


# ============================================================
# SAVE METHOD OUTPUT
# ============================================================

output_file <- file.path(
    output_dir,
    paste0(
        "aldex2_",
        dataset_label,
        "_results.csv"
    )
)

write.csv(
    aldex2_table,
    output_file,
    row.names = FALSE
)


# ============================================================
# SAVE RUN SUMMARY
# ============================================================

summary_file <- file.path(
    output_dir,
    paste0(
        "aldex2_",
        dataset_label,
        "_summary.txt"
    )
)

significant_count <- sum(
    aldex2_table$Significant,
    na.rm = TRUE
)

summary_text <- c(
    "ALDEx2 METHOD SUMMARY",
    "=====================",
    "",
    paste("Dataset:", dataset_label),
    paste("ALDEx2 version:", packageVersion("ALDEx2")),
    paste("Samples:", ncol(counts)),
    paste("Taxa analyzed:", nrow(aldex2_table)),
    "",
    "Comparison:",
    "Case vs Control",
    "",
    "Reference group:",
    "Control",
    "",
    "Monte Carlo samples:",
    "128",
    "",
    "CLR denominator:",
    "all",
    "",
    "Primary statistical test:",
    "Wilcoxon rank-sum",
    "",
    "P-value adjustment:",
    "Benjamini-Hochberg (BH)",
    "",
    "Significance threshold:",
    "Adjusted p-value < 0.05",
    "",
    paste("Significant taxa:", significant_count),
    "",
    "Direction:",
    "Positive difference = higher in Case",
    "Negative difference = higher in Control",
    "",
    "Welch and Wilcoxon ALDEx2 results are both retained in the output.",
    "QC filtering was completed before ALDEx2."
)

writeLines(
    summary_text,
    summary_file
)


# ============================================================
# TERMINAL OUTPUT
# ============================================================

cat("\nALDEx2 COMPLETED\n")
cat("================\n")

cat(
    "Taxa analyzed:",
    nrow(aldex2_table),
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
