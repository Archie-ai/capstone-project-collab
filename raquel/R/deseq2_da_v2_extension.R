# DESeq2: Version 2 DA simulations
# BIOT 670I - Microbiome Abundance Consensus Machine

suppressPackageStartupMessages(library(DESeq2))

n_simulations <- 100
alpha <- 0.05

input_dir <- "data/da/simulations_v2"
output_dir <- "results/da/deseq2_v2_extension"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

performance_results <- list()

for (sim in 101:200) {

  sim_name <- sprintf("simulation_%03d", sim)
  sim_dir <- file.path(input_dir, sim_name)

  # Load count data, metadata, and ground truth
  counts <- read.csv(
    file.path(sim_dir, "da_counts.csv"),
    row.names = 1,
    check.names = FALSE
  )

  metadata <- read.csv(
    file.path(sim_dir, "da_metadata.csv")
  )

  truth <- read.csv(
    file.path(sim_dir, "da_ground_truth.csv")
  )

  # Align metadata with count columns
  metadata <- metadata[
    match(colnames(counts), metadata$SampleID),
    ,
    drop = FALSE
  ]

  if (anyNA(metadata$SampleID)) {
    stop("Sample metadata mismatch in ", sim_name)
  }

  rownames(metadata) <- metadata$SampleID

  # Set Control as reference group
  metadata$Group <- factor(
    metadata$Group,
    levels = c("Control", "Case")
  )

  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = metadata,
    design = ~ Group
  )

  # Run DESeq2
  dds <- DESeq(dds, quiet = TRUE)

  # Compare Case against Control
  res <- results(
    dds,
    contrast = c("Group", "Case", "Control"),
    alpha = alpha
  )

  # Standardize results
  res <- as.data.frame(res)
  res$Taxon <- rownames(res)

  final_results <- data.frame(
    Taxon = res$Taxon,
    Effect = res$log2FoldChange,
    Adjusted_PValue = res$padj
  )

  final_results$Significant <- (
    !is.na(final_results$Adjusted_PValue) &
    final_results$Adjusted_PValue < alpha
  )

  # Match results with ground truth
  final_results <- merge(
    truth,
    final_results,
    by = "Taxon",
    all.x = TRUE,
    sort = FALSE
  )

  final_results$Significant[
    is.na(final_results$Significant)
  ] <- FALSE

  true_da <- as.logical(final_results$True_DA)
  detected <- final_results$Significant

  # Calculate performance metrics
  TP <- sum(detected & true_da)
  FP <- sum(detected & !true_da)
  FN <- sum(!detected & true_da)
  TN <- sum(!detected & !true_da)

  discoveries <- TP + FP

  FDP <- if (discoveries > 0) {
    FP / discoveries
  } else {
    0
  }

  power <- TP / (TP + FN)
  specificity <- TN / (TN + FP)

  performance_results[[sim - 100]] <- data.frame(
    Simulation = sim,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    Discoveries = discoveries,
    FDP = FDP,
    Power = power,
    Sensitivity = power,
    Specificity = specificity,
    Taxa_Tested = sum(!is.na(res$padj))
  )

  # Save individual simulation results
  write.csv(
    final_results,
    file.path(
      output_dir,
      paste0(sim_name, "_results.csv")
    ),
    row.names = FALSE
  )

  cat(
  "Simulation", sim, "/",
  200, "analyzed.\n"
)
}

# Combine performance results
performance <- do.call(
  rbind,
  performance_results
)

write.csv(
  performance,
  file.path(
    output_dir,
    "performance_by_simulation.csv"
  ),
  row.names = FALSE
)

# Create summary statistics
metrics <- setdiff(
  names(performance),
  "Simulation"
)

summary_table <- data.frame(
  Metric = metrics,
  Mean = sapply(
    performance[metrics],
    mean,
    na.rm = TRUE
  ),
  SD = sapply(
    performance[metrics],
    sd,
    na.rm = TRUE
  ),
  Min = sapply(
    performance[metrics],
    min,
    na.rm = TRUE
  ),
  Max = sapply(
    performance[metrics],
    max,
    na.rm = TRUE
  ),
  row.names = NULL
)

write.csv(
  summary_table,
  file.path(
    output_dir,
    "performance_summary.csv"
  ),
  row.names = FALSE
)

cat("\nAll DESeq2 analyses completed!\n")
print(summary_table, row.names = FALSE)
cat("\nResults saved to:", output_dir, "\n")
