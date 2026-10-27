# Peggy Hall — Microbiome Differential-Abundance Analysis

**BIOT 670I Bioinformatics Capstone | University of Maryland Global Campus**

Part of the collaborative [Microbiome Abundance Consensus Machine (MACM)](../../) project.

## Overview

My capstone work examines how statistical choices affect microbiome differential-abundance findings. I contributed two complementary analyses:

- **Real-data case study:** comparison of four differential-abundance methods in a colorectal cancer (CRC) versus healthy-control microbiome dataset (`crc_baxter`).
- **Simulation-based evaluation:** assessment of method performance and consensus rules when the true differential-abundance status of taxa is known.

The emphasis is on **data quality, statistical reliability, reproducibility, and cautious biological interpretation**—not treating statistical significance or agreement between methods as proof of biological truth.

## Methods and tools

- **Methods:** CLR/Wilcoxon, ALDEx2 (Welch-primary in my workflow), ANCOM-BC2, and DESeq2
- **Computational tools:** Python and R
- **Evaluation:** sensitivity, empirical false discovery rate (FDR), specificity, F1 score, agreement, and consensus across methods
- **Workflow practices:** input validation, preprocessing documentation, standardized result tables, simulation checkpoints, and visual summaries

## 1. CRC microbiome case study

**Dataset:** MicrobiomeHD `crc_baxter`, comparing colorectal cancer samples with healthy controls.

| Dataset / analysis characteristic | Value |
| --- | ---: |
| Samples | 292 (120 CRC; 172 healthy) |
| OTUs analyzed | 1,133 |
| Significance threshold | BH-adjusted p < 0.05 |

### Differential-abundance results

| Method | Significant OTUs |
| --- | ---: |
| CLR/Wilcoxon | 24 |
| ALDEx2 (Welch) | 0 |
| ANCOM-BC2 (standard significance) | 486 |
| DESeq2 | 20 |

Despite using the same dataset, the four methods produced substantially different discovery sets.

**Cross-method consensus:**

| Support threshold | OTUs |
| --- | ---: |
| At least 1 method | 494 |
| At least 2 methods | 34 |
| At least 3 methods | 2 |
| All 4 methods | 0 |

These counts describe agreement, **not** ground-truth accuracy. The real dataset does not establish which detected taxa are truly differential; the analysis also does not establish effect-direction agreement or adjust for clinical covariates. ANCOM-BC2 counts use BH-adjusted significance rather than confirmed pseudo-count-robust calls.

![Pairwise method overlap (Jaccard similarity)](data/real/crc_baxter/DA_results/crc_baxter_jaccard_heatmap.png)

**Explore the work:** [CRC results and consensus tables](data/real/crc_baxter/DA_results/) · [Final CRC interpretation](data/real/crc_baxter/DA_results/crc_baxter_final_interpretation.txt)

## 2. Simulation and consensus evaluation

Simulated microbiome data allow method outputs to be compared against known truth. My Phase 3 outputs evaluate performance across **27 simulated conditions** and compare consensus strategies requiring agreement from two, three, or all four methods.

| Consensus rule | Mean sensitivity | Mean empirical FDR |
| --- | ---: | ---: |
| 2 of 4 | 25.01% | 2.58% |
| 3 of 4 | 20.33% | 0.47% |
| 4 of 4 | 14.16% | 0.05% |

*Source: Phase 3 consensus summary; values are averages within the evaluated simulation design, not universal estimates.*

Stricter consensus lowered the empirical false discovery rate while also lowering sensitivity. The Phase 3 output recommends **3-of-4** overall for its tested conditions; that recommendation should **not** be generalized to every microbiome dataset.

![Sensitivity versus false-discovery rate across consensus approaches](data/simulated/Phase3/Consensus/Phase3_Final_Figures/Figure_9N1_Sensitivity_vs_FDR.png)

**Explore the work:** [Phase 3 consensus summaries, figures, and checkpoints](data/simulated/Phase3/Consensus/) · [Final consensus recommendation (CSV)](data/simulated/Phase3/Consensus/phase3_final_consensus_recommendation.csv)

## Key takeaways

1. Method selection can substantially change which microbial taxa are labeled significant.
2. Real-data concordance is not the same as accuracy against independently established biological truth.
3. Consensus offers a way to prioritize findings, but stricter agreement trades sensitivity for fewer false discoveries.
4. Simulations and real-data studies answer different, complementary questions.
5. Biological context, transparent preprocessing, and careful interpretation remain central to trustworthy bioinformatics.

## Directory guide

| Location | What it contains |
| --- | --- |
| [`data/real/crc_baxter/`](data/real/crc_baxter/) | CRC dataset metadata, filtered data, validation summaries, method outputs, and agreement analyses |
| [`data/simulated/`](data/simulated/) | Simulation calibration outputs, method-performance tables, and figures |
| [`data/simulated/Phase3/Consensus/`](data/simulated/Phase3/Consensus/) | Consensus evaluation, condition-level summaries, and final figures |
| [`scripts/`](scripts/) | Current script-directory notes; standalone source scripts are not included here |

**Reproducibility note:** This directory currently emphasizes retained analysis outputs and documentation. The `scripts/` folder does not yet provide a complete executable workflow. Results should be interpreted alongside the preprocessing details, simulation conditions, software choices, and limitations of the [team README](../../).

## Collaboration and scope

This analysis was completed as part of a four-person UMGC capstone team. Other contributors analyzed additional real microbiome datasets and simulation scenarios. See the [main MACM project README](../../) for team-wide methodology, comparisons, contributors, and project limitations.

This repository is academic research; its results are not clinical or causal conclusions.
