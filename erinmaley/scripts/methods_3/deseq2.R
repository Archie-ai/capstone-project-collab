#!/usr/bin/env Rscript

# ============================================================
# DESeq2
# ============================================================
#
# This script runs DESeq2 on microbiome count data that passed shared QC.
#
# This script was developed with assistance from OpenAI's ChatGPT
#

suppressPackageStartupMessages({
    library(DESeq2)
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
            "Rscript deseq2.R",
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

cat("\nRunning DESeq2 for dataset:", dataset_label, "\n\n")

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

# Match metadata order to the count-table columns.
metadata <- metadata[
    colnames(counts),
    ,
    drop = FALSE
]


# ============================================================
# DEFINE CASE / CONTROL COMPARISON
# ============================================================

# Keep the same direction convention used for ANCOM-BC.
#
# Control = reference group
# Positive log2 fold change = higher in Case
# Negative log2 fold change = higher in Control

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
# PREPARE COUNTS
# ============================================================

# DESeq2 requires a matrix of non-negative integer counts.
count_matrix <- as.matrix(counts)

storage.mode(count_matrix) <- "integer"


# ============================================================
# BUILD DESEQ2 DATASET
# ============================================================

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = ~ Group
)


# ============================================================
# RUN DESEQ2
# ============================================================

# DESeq2 estimates size factors, dispersion, and the model internally.
dds <- DESeq(
    dds,
    quiet = TRUE
)


# ============================================================
# EXTRACT CASE VS CONTROL RESULTS
# ============================================================

result <- results(
    dds,
    contrast = c(
        "Group",
        "Case",
        "Control"
    ),
    alpha = 0.05
)


# ============================================================
# CREATE METHOD OUTPUT TABLE
# ============================================================

deseq2_table <- data.frame(

    Taxon = rownames(result),

    BaseMean =
        result$baseMean,

    Log2FoldChange =
        result$log2FoldChange,

    StandardError =
        result$lfcSE,

    Statistic =
        result$stat,

    PValue =
        result$pvalue,

    AdjustedPValue =
        result$padj,

    stringsAsFactors = FALSE
)


# A taxon is significant when its adjusted p-value is below 0.05.
# NA adjusted p-values are not treated as significant.

deseq2_table$Significant <- (
    !is.na(deseq2_table$AdjustedPValue)
    & deseq2_table$AdjustedPValue < 0.05
)


# Add direction for easier interpretation later.

deseq2_table$Direction <- ifelse(
    is.na(deseq2_table$Log2FoldChange),
    "Not_available",
    ifelse(
        deseq2_table$Log2FoldChange > 0,
        "Higher_in_Case",
        ifelse(
            deseq2_table$Log2FoldChange < 0,
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
        "deseq2_",
        dataset_label,
        "_results.csv"
    )
)

write.csv(
    deseq2_table,
    output_file,
    row.names = FALSE
)


# ============================================================
# SAVE RUN SUMMARY
# ============================================================

summary_file <- file.path(
    output_dir,
    paste0(
        "deseq2_",
        dataset_label,
        "_summary.txt"
    )
)

significant_count <- sum(
    deseq2_table$Significant,
    na.rm = TRUE
)

na_adjusted_count <- sum(
    is.na(deseq2_table$AdjustedPValue)
)

summary_text <- c(
    "DESeq2 METHOD SUMMARY",
    "=====================",
    "",
    paste("Dataset:", dataset_label),
    paste("DESeq2 version:", packageVersion("DESeq2")),
    paste("Samples:", ncol(counts)),
    paste("Taxa analyzed:", nrow(deseq2_table)),
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
    paste("Significant taxa:", significant_count),
    paste("Taxa with NA adjusted p-value:", na_adjusted_count),
    "",
    "Direction:",
    "Positive log2 fold change = higher in Case",
    "Negative log2 fold change = higher in Control",
    "",
    "QC filtering was completed before DESeq2.",
    "DESeq2 performed its own method-specific size-factor normalization."
)

writeLines(
    summary_text,
    summary_file
)


# ============================================================
# TERMINAL OUTPUT
# ============================================================

cat("\nDESeq2 COMPLETED\n")
cat("================\n")

cat(
    "Taxa analyzed:",
    nrow(deseq2_table),
    "\n"
)

cat(
    "Significant taxa:",
    significant_count,
    "\n"
)

cat(
    "NA adjusted p-values:",
    na_adjusted_count,
    "\n"
)

cat("\nOutput:\n")
cat(" ", output_file, "\n")

cat("\nSummary:\n")
cat(" ", summary_file, "\n\n")
