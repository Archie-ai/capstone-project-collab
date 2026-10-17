#!/usr/bin/env bash
# Developed with assistance from OpenAI's ChatGPT.
set -euo pipefail

BASE=erinmaley/data/benchmark_11

for i in $(seq -w 1 5); do
    SIM="$BASE/simulation_00${i}"
    QC="$SIM/qc"
    OUT="$SIM/method_outputs"

    echo "===== Processing simulation ${i}/5 ====="

    python3 erinmaley/scripts/quality_control.py \
      --counts "$SIM/counts.csv" \
      --metadata "$SIM/metadata.csv" \
      --output-dir "$QC" \
      --dataset-label signal

    mkdir -p "$OUT"/{aldex2,ancombc,clr_wilcoxon,deseq2}

    Rscript erinmaley/scripts/methods_3/aldex2.R \
      "$QC/qc_signal_counts.csv" \
      "$QC/qc_signal_metadata.csv" "$OUT/aldex2" signal

    Rscript erinmaley/scripts/methods_3/ancombc.R \
      "$QC/qc_signal_counts.csv" \
      "$QC/qc_signal_metadata.csv" "$OUT/ancombc" signal

    python3 erinmaley/scripts/methods_3/clr_wilcoxon.py \
      --counts "$QC/qc_signal_counts.csv" \
      --metadata "$QC/qc_signal_metadata.csv" \
      --output-dir "$OUT/clr_wilcoxon" --dataset-label signal

    Rscript erinmaley/scripts/methods_3/deseq2.R \
      "$QC/qc_signal_counts.csv" \
      "$QC/qc_signal_metadata.csv" "$OUT/deseq2" signal

    python3 erinmaley/scripts/standard_results_4/standardize_results.py \
      --input-dir "$OUT" \
      --output-dir "$SIM/standard_results" \
      --dataset-label signal

    python3 erinmaley/scripts/analysis_5/analyze_results.py \
      --standard-results "$SIM/standard_results/standardized_signal_all_methods.csv" \
      --output-dir "$SIM/analysis" \
      --dataset-label signal \
      --ground-truth "$SIM/ground_truth.csv"

    echo "===== Simulation ${i} complete ====="
done

echo "ALL FIVE BENCHMARKS COMPLETE"
