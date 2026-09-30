# BIOT670I Microbiome Group Project

This repository contains the code, analyses, documentation, and results for the BIOT670I microbiome group project (Capstone Project).

## Project Status

- Project planning: Complete
- Simulation framework: Complete
- Phase 1 – Null calibration: Complete
- Phase 2 – Signal simulation: Complete
- Phase 3 – Extended/real-data analysis: Next phase
- Final integration and interpretation: Planned

## Repository Structure

```text
Microbiome_Group_Project/
├── data/
├── notebooks/
│   ├── BIOT670I_Microbiome_Project.ipynb
│   └── README.md
├── results/
│   ├── README.md
│   ├── phase1_final_calibration_summary.csv
│   ├── phase1_raw_type1_error.png
│   ├── phase1_any_bh_false_discovery.png
│   ├── phase2_performance_summary_8F.csv
│   ├── phase2_sensitivity_plot_8G.png
│   ├── phase2_empirical_FDR_plot_8H.png
│   └── phase2_combined_performance_plot_8I.png
├── scripts/
├── .gitignore
└── README.md
```

## Phase 1 Results

Phase 1 evaluates the statistical analysis pipeline under the null hypothesis using repeated microbiome simulations.

Selected Phase 1 results, including empirical Type I error, Benjamini-Hochberg false-discovery behavior, and the final calibration summary, are available in the [`results`](results/) directory.


## Phase 2 Results

Phase 2 evaluates differential-abundance method performance under simulated microbiome datasets containing known signal effects.

The analysis compares CLR/Wilcoxon, ALDEx2 Welch, ANCOM-BC2, and DESeq2 across increasing signal strengths (1.5×, 2×, and 4×). Performance is evaluated primarily using sensitivity and empirical false discovery rate (FDR).

Selected Phase 2 outputs are available in the [`results`](results/) directory:

- [`phase2_performance_summary_8F.csv`](results/phase2_performance_summary_8F.csv) – summary performance metrics across methods and signal strengths
- [`phase2_sensitivity_plot_8G.png`](results/phase2_sensitivity_plot_8G.png) – sensitivity across signal strengths
- [`phase2_empirical_FDR_plot_8H.png`](results/phase2_empirical_FDR_plot_8H.png) – empirical FDR across signal strengths
- [`phase2_combined_performance_plot_8I.png`](results/phase2_combined_performance_plot_8I.png) – combined sensitivity and empirical FDR comparison

Additional methodological details and interpretation are documented in [`results/README.md`](results/README.md).
