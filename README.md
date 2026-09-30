# BIOT670I Microbiome Group Project

This repository contains the code, analyses, documentation, and results for the BIOT670I microbiome capstone project.

The project evaluates differential-abundance analysis methods using simulated microbiome data, beginning with null calibration and progressing to planted-signal simulations before extended and real-data analysis.

## Project Status

- Project planning: Complete
- Simulation framework: Complete
- Phase 1 - Null calibration: Complete
- Phase 2 - Signal simulation and performance evaluation: Complete
- Phase 3 - Extended/real-data analysis: Next phase
- Final integration and interpretation: Planned

## Repository Structure

```text
Microbiome_Group_Project/
├── data/
│   ├── simulated/   # optional saved simulation datasets; most simulations are generated in the notebook
│   └── real/        # real microbiome datasets for Phase 3
├── notebooks/
│   ├── biot670i_microbiome_project.py
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
├── scripts/         # reusable standalone analysis code as the workflow is modularized
├── .gitignore
└── README.md
```

The primary analysis is maintained in the Jupyter notebook, with a Python export included in the `notebooks/` directory for easier code review and version tracking.

The `data/simulated/` and `scripts/` directories are intentionally lightweight at this stage. Phase 1 and Phase 2 simulations are generated directly in the notebook, while reusable standalone functions can be moved into `scripts/` later as the project workflow is modularized. Simulation-derived summary tables and figures belong in `results/`, not `data/simulated/`.

## Phase 1: Null Calibration

Phase 1 evaluates the statistical analysis pipeline under the null hypothesis, where no true differential-abundance signal is present.

Repeated microbiome simulations were used to assess calibration behavior, including:

- empirical raw Type I error relative to the nominal 0.05 level
- Benjamini-Hochberg false-discovery behavior
- robustness and diagnostic summaries across repeated simulation runs

Selected Phase 1 tables and figures are available in the [`results`](results/) directory.

## Phase 2: Signal Simulation

Phase 2 evaluates method performance after known differential-abundance signals are introduced into the simulated microbiome datasets.

Four methods were compared:

- CLR/Wilcoxon
- ALDEx2 Welch
- ANCOM-BC2
- DESeq2

Signal strengths of 1.5x, 2x, and 4x were evaluated across repeated simulations. Performance was summarized primarily using sensitivity and empirical false discovery rate (FDR).

Across the simulated conditions, sensitivity generally increased as signal strength increased. The methods also differed in their false-discovery behavior, illustrating a trade-off between detection sensitivity and FDR control rather than a single method performing best under every condition.

Selected Phase 2 outputs include:

- [`phase2_performance_summary_8F.csv`](results/phase2_performance_summary_8F.csv) - summary performance metrics across methods and signal strengths
- [`phase2_sensitivity_plot_8G.png`](results/phase2_sensitivity_plot_8G.png) - sensitivity across signal strengths
- [`phase2_empirical_FDR_plot_8H.png`](results/phase2_empirical_FDR_plot_8H.png) - empirical FDR across signal strengths
- [`phase2_combined_performance_plot_8I.png`](results/phase2_combined_performance_plot_8I.png) - combined sensitivity and empirical FDR comparison

Additional methodological details and interpretation are documented in [`results/README.md`](results/README.md).

## Next Phase

Phase 3 will extend the completed simulation framework toward additional analyses and/or real microbiome data. The existing Phase 1 and Phase 2 workflow provides the calibration and performance baseline for interpreting those later results.

## Reproducibility

The repository is organized so that analysis code, selected outputs, and documentation remain version controlled while large or temporary intermediate files can be excluded through `.gitignore`.
