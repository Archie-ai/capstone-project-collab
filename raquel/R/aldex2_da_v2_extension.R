library(ALDEx2)

# Simulation settings
n_simulations <- 100
alpha <- 0.05
mc_samples <- 128

# Input and output directories
input_dir <- "data/da/simulations_v2"
output_dir <- "results/da/aldex2_v2_extension"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Store performance metrics
performance_results <- list()

for (sim in 101:200) {

  # Reproducible Monte Carlo sampling
  set.seed(1000 + sim)

  sim_dir <- file.path(
    input_dir,
    sprintf("simulation_%03d", sim)
  )

  # Load counts, metadata, and ground truth
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

  # Match metadata to count table
  metadata <- metadata[
    match(colnames(counts), metadata$SampleID),
  ]

  if (anyNA(metadata$SampleID)) {
    stop("Metadata mismatch in simulation ", sim)
  }

  conditions <- metadata$Group

  # Run ALDEx2 CLR transformation
  clr <- aldex.clr(
    counts,
    conditions,
    mc.samples = mc_samples,
    denom = "all",
    verbose = FALSE
  )

  # Run Wilcoxon test
  test <- aldex.ttest(
    clr,
    paired.test = FALSE,
    verbose = FALSE
  )

  # Calculate effect sizes
  effect <- aldex.effect(
    clr,
    verbose = FALSE
  )

  # Combine results
  results <- cbind(test, effect)

  results$Taxon <- rownames(results)
  results$Adjusted_PValue <- results$wi.eBH

  results$Significant <- (
    !is.na(results$Adjusted_PValue) &
    results$Adjusted_PValue < alpha
  )

  # Merge with all ground-truth taxa
  # Excluded taxa are treated as not significant
  results <- merge(
    truth,
    results,
    by = "Taxon",
    all.x = TRUE,
    sort = FALSE
  )

  results$Significant[
    is.na(results$Significant)
  ] <- FALSE

  if (anyNA(results$True_DA)) {
    stop("Missing ground truth in simulation ", sim)
  }

  # Confusion matrix
  tp <- sum(
    results$Significant &
    results$True_DA == 1
  )

  fp <- sum(
    results$Significant &
    results$True_DA == 0
  )

  fn <- sum(
    !results$Significant &
    results$True_DA == 1
  )

  tn <- sum(
    !results$Significant &
    results$True_DA == 0
  )

  # Performance metrics
  discoveries <- tp + fp

  fdp <- if (discoveries > 0) {
    fp / discoveries
  } else {
    0
  }

  power <- tp / (tp + fn)
  sensitivity <- power
  specificity <- tn / (tn + fp)

  performance_results[[sim - 100]] <- data.frame(
    Simulation = sim,
    TP = tp,
    FP = fp,
    FN = fn,
    TN = tn,
    Discoveries = discoveries,
    FDP = fdp,
    Power = power,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Taxa_Tested = sum(!is.na(results$Adjusted_PValue))
  )

  # Save individual simulation results
  write.csv(
    results,
    file.path(
      output_dir,
      sprintf("simulation_%03d_results.csv", sim)
    ),
    row.names = FALSE
  )

  cat(
    sprintf("Simulation %d/200 analyzed.\n", sim)
  )
}

# Combine performance results
performance_df <- do.call(
  rbind,
  performance_results
)

write.csv(
  performance_df,
  file.path(output_dir, "performance_by_simulation.csv"),
  row.names = FALSE
)

# Calculate summary statistics
metrics <- c(
  "TP", "FP", "FN", "TN",
  "Discoveries", "FDP", "Power",
  "Sensitivity", "Specificity",
  "Taxa_Tested"
)

summary_df <- data.frame(
  Metric = metrics,
  Mean = sapply(
    performance_df[metrics], mean
  ),
  SD = sapply(
    performance_df[metrics], sd
  ),
  Min = sapply(
    performance_df[metrics], min
  ),
  Max = sapply(
    performance_df[metrics], max
  ),
  row.names = NULL
)

write.csv(
  summary_df,
  file.path(output_dir, "performance_summary.csv"),
  row.names = FALSE
)

cat("\nAll 100 ALDEx2 analyses completed!\n")
print(summary_df, row.names = FALSE)
cat("\nResults saved to:", output_dir, "\n")
