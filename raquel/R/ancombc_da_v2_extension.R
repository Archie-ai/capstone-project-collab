# ANCOM-BC2: Version 2 DA simulations
# BIOT 670I - Microbiome Abundance Consensus Machine

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(TreeSummarizedExperiment)
})

n_simulations <- 100
alpha <- 0.05

input_dir <- "data/da/simulations_v2"
output_dir <- "results/da/ancombc_v2_extension"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

performance_results <- list()

for (sim in 101:200) {

  sim_name <- sprintf("simulation_%03d", sim)
  sim_dir <- file.path(input_dir, sim_name)

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

  # Align sample metadata with count columns
  metadata <- metadata[
    match(colnames(counts), metadata$SampleID),
    ,
    drop = FALSE
  ]

  if (anyNA(metadata$SampleID)) {
    stop("Sample metadata mismatch in ", sim_name)
  }

  # Ensure Control is the reference group
  metadata$Group <- factor(
    metadata$Group,
    levels = c("Control", "Case")
  )

  # Remove taxa with zero counts in all samples
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]

  # Create TreeSummarizedExperiment
  tse <- TreeSummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    colData = S4Vectors::DataFrame(metadata)
  )

  # Run ANCOM-BC2
  output <- ancombc2(
    data = tse,
    assay_name = "counts",
    fix_formula = "Group",
    p_adj_method = "BH",
    prv_cut = 0,
    lib_cut = 0,
    group = "Group",
    struc_zero = FALSE,
    neg_lb = FALSE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE,
    verbose = FALSE
  )

  results <- output$res

  # Identify Case vs Control coefficient
  lfc_col <- "lfc_GroupCase"
  q_col <- "q_GroupCase"

  if (!all(c("taxon", lfc_col, q_col) %in% names(results))) {
    stop(
      "Expected ANCOM-BC2 result columns missing in ",
      sim_name,
      ": ",
      paste(names(results), collapse = ", ")
    )
  }

robust_col <- "diff_robust_GroupCase"
passed_col <- "passed_ss_GroupCase"

if (!all(c(robust_col, passed_col) %in% names(results))) {
  stop(
    "ANCOM-BC2 robustness columns missing: ",
    paste(names(results), collapse = ", ")
  )
}

final_results <- data.frame(
  Taxon = results$taxon,
  Effect = results[[lfc_col]],
  Adjusted_PValue = results[[q_col]],
  Passed_Sensitivity = results[[passed_col]],
  Robust_Significant = results[[robust_col]]
)

  final_results$Significant <- (
    !is.na(final_results$Adjusted_PValue) &
    final_results$Adjusted_PValue < alpha
  )

  # Compare findings with planted ground truth
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

  # Ground-truth column should identify planted DA taxa
  truth_col <- setdiff(
    names(truth),
    "Taxon"
  )

  if (length(truth_col) != 1) {
    stop("Unexpected ground-truth columns in ", sim_name)
  }

  true_da <- as.logical(final_results[[truth_col]])
  detected <- final_results$Significant

  if (anyNA(true_da)) {
    stop("Missing ground-truth values in ", sim_name)
  }

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

  power <- if ((TP + FN) > 0) {
    TP / (TP + FN)
  } else {
    NA_real_
  }

  specificity <- if ((TN + FP) > 0) {
    TN / (TN + FP)
  } else {
    NA_real_
  }

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
    Taxa_Tested = nrow(results)
  )

  write.csv(
    final_results,
    file.path(output_dir, paste0(sim_name, "_results.csv")),
    row.names = FALSE
  )

  cat("Simulation", sim, "/ 200 analyzed.\n")
}

# Save performance metrics for all simulations
performance <- do.call(rbind, performance_results)

write.csv(
  performance,
  file.path(output_dir, "performance_by_simulation.csv"),
  row.names = FALSE
)

metrics <- setdiff(names(performance), "Simulation")

summary_table <- data.frame(
  Metric = metrics,
  Mean = sapply(performance[metrics], mean, na.rm = TRUE),
  SD = sapply(performance[metrics], sd, na.rm = TRUE),
  Min = sapply(performance[metrics], min, na.rm = TRUE),
  Max = sapply(performance[metrics], max, na.rm = TRUE),
  row.names = NULL
)

write.csv(
  summary_table,
  file.path(output_dir, "performance_summary.csv"),
  row.names = FALSE
)

cat("\nAll ANCOM-BC2 analyses completed!\n")
print(summary_table, row.names = FALSE)
cat("\nResults saved to:", output_dir, "\n")
