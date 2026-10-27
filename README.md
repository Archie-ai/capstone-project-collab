# Microbiome Abundance Consensus Machine

### Comparative Evaluation of Differential-Abundance Methods in Simulated and Real Microbiome Data

**University of Maryland Global Campus**  
**BIOT 670I — Capstone Project**

**Team Members:**  
Architha Sankar • Erin Maley • Peggy Hall • Raquel De Leon

---

## Project Overview

Differential-abundance (DA) analysis is widely used in microbiome research to identify microbial taxa whose abundance differs between biological groups, such as disease and control populations.

However, microbiome datasets present several statistical challenges. They are typically **compositional**, often highly **sparse**, contain many zero observations, and can vary substantially in sequencing depth. Different DA methods address these characteristics using different transformations, normalization procedures, probability models, and statistical assumptions.

As a result, applying different statistical methods to the same microbiome dataset can produce substantially different sets of significant taxa.

The **Microbiome Abundance Consensus Machine (MACM)** was developed to systematically investigate this problem.

The project compares four differential-abundance approaches:

- **CLR/Wilcoxon**
- **ALDEx2**
- **ANCOM-BC2**
- **DESeq2**

The analysis combines:

1. **Simulated microbiome datasets with known ground truth**
2. **Repeated null simulations**
3. **Known-signal simulations**
4. **Four real microbiome datasets**
5. **Cross-method agreement and consensus analysis**

The objective is not to identify one universally superior DA method. Instead, the project examines:

> **How data characteristics and statistical assumptions influence method behavior, when methods agree or disagree, and whether consensus among multiple methods can improve the reliability of microbiome differential-abundance interpretation.**

---

# Research Question

### Central Research Question

**Under what microbiome data conditions do differential-abundance methods agree, disagree, or fail, and can a consensus strategy produce more reliable biological conclusions?**

---

# Project Objectives

The project was organized around six major objectives:

1. Develop a reproducible computational framework integrating four differential-abundance approaches.

2. Generate simulated microbiome datasets in which the truly differential taxa are explicitly known.

3. Evaluate method performance using metrics such as:

   - sensitivity;
   - statistical power;
   - false-positive behavior;
   - Type I error;
   - false discovery proportion;
   - false discovery rate;
   - specificity; and
   - related performance measures.

4. Apply the four methods to real microbiome datasets and compare:

   - significant-taxon discovery burden;
   - pairwise overlap;
   - Jaccard similarity;
   - effect-direction agreement; and
   - cross-method concordance.

5. Evaluate consensus strategies requiring agreement from:

   - **2 of 4 methods**
   - **3 of 4 methods**
   - **4 of 4 methods**

6. Examine how microbiome data properties such as sample size, sparsity, prevalence, sequencing depth, and effect size influence method behavior.

---

# Why This Problem Matters

A statistically significant microbial association can influence biological interpretation, biomarker development, downstream laboratory experiments, or disease-mechanism hypotheses.

However, if the result changes substantially depending on which DA method is selected, then statistical significance alone may not be sufficient evidence of robustness.

Microbiome DA methods differ in:

- normalization;
- compositional treatment;
- zero handling;
- variance modeling;
- statistical testing;
- filtering;
- effect-size estimation; and
- multiple-testing behavior.

Therefore, the project focuses not only on **which taxa are significant**, but also on **how reproducible those findings are across analytical approaches**.

---

# Two Complementary Experimental Arms

## 1. Simulated Data — Ground Truth Known

Simulations provide a controlled environment in which the true DA status of every taxon is known before statistical analysis.

This allows direct classification of results as:

- True Positive (TP)
- False Positive (FP)
- True Negative (TN)
- False Negative (FN)

Simulation therefore permits direct measurement of statistical accuracy.

Examples of measurable outcomes include:

- Sensitivity
- Power
- Specificity
- Type I error
- False-positive rate
- False discovery proportion (FDP)
- False discovery rate (FDR)

### Key Principle

**Simulation evaluates statistical performance because the truth is known.**

---

## 2. Real Microbiome Data — Ground Truth Unknown

For real microbiome datasets, the complete biological truth is not independently known.

Therefore, real-data discoveries cannot automatically be classified as true positives or false positives.

Instead, the analysis focuses on:

- discovery burden;
- discovery rate;
- significant-taxon overlap;
- Jaccard similarity;
- effect-direction agreement; and
- consensus support.

### Key Principle

**Real-data analysis evaluates concordance and reproducibility rather than direct statistical accuracy.**

---

# Differential-Abundance Methods

| Method | General Approach | Role in This Project |
|---|---|---|
| **CLR/Wilcoxon** | Centered log-ratio transformation followed by a nonparametric rank test | Transparent compositional baseline |
| **ALDEx2** | Dirichlet Monte Carlo sampling and log-ratio compositional analysis | Explicitly propagates uncertainty in compositional count data |
| **ANCOM-BC2** | Bias-corrected regression framework | Estimates abundance differences while correcting sampling/compositional bias |
| **DESeq2** | Negative-binomial generalized linear model | Count-modeling contrast originally developed for RNA-seq |

The effect estimates produced by these methods are **not directly interchangeable**.

For example:

- CLR differences;
- ALDEx2 effect estimates;
- ANCOM log fold changes; and
- DESeq2 log2 fold changes

represent different statistical quantities.

Cross-method comparisons therefore emphasize:

**significance + effect direction**

rather than direct comparison of effect-size magnitude.

---

# Overall Analytical Workflow

```text
SIMULATED DATA                         REAL MICROBIOME DATA
Known ground truth                     Unknown biological truth
       |                                         |
       +----------------+------------------------+
                        |
                Input validation / QC
                        |
       +----------------+----------------+
       |                |               |
 CLR/Wilcoxon        ALDEx2         ANCOM-BC2
                                         |
                                      DESeq2
       +----------------+----------------+
                        |
              Standardized outputs
                        |
       +----------------+----------------+
       |                                 |
Simulation performance             Real-data agreement
TP / FP / TN / FN                  Discovery burden
Sensitivity / power                Jaccard similarity
Type I error / FDP / FDR           Effect direction
       |                                 |
       +----------------+----------------+
                        |
               CONSENSUS ANALYSIS
             2-of-4 / 3-of-4 / 4-of-4
                        |
               Integrated interpretation
```

---

# Real Microbiome Datasets

Four real datasets were evaluated.

| Dataset | Biological Condition | Primary Team Analysis |
|---|---|---|
| `par_scheperjans` | Parkinson's disease | Architha |
| `cdi_schubert` | *Clostridioides difficile* infection | Erin |
| `crc_baxter` | Colorectal cancer | Peggy |
| `cdi_vincent` | *Clostridioides difficile* infection | Raquel |

Because biological ground truth is unknown, these datasets were used primarily to examine **method disagreement and concordance**.

---

# Repository Organization

The repository is organized as a collaborative team workspace.

```text
capstone-project-collab/
│
├── architha/
│   ├── scripts/
│   └── results/
│
├── erinmaley/
│   └── data/
│
├── peggyhall/
│   └── data/
│
├── raquel/
│   └── data/
│
└── README.md
```

Each contributor's directory contains the code, intermediate outputs, simulation results, and/or real-data analyses associated with that contributor's workflow.

Important consolidated and benchmark locations include:

```text
architha/results/team_consolidation/null_vs_real/

erinmaley/data/null_benchmark_12/

erinmaley/data/benchmark_11/

raquel/data/da/consensus_v2_200/

peggyhall/data/simulated/Phase3/Consensus/
```

The project was developed collaboratively, so individual analyses may contain their own supporting directory structures and scripts.

---

# Computational Environment

The project combines Python and R.

## Core Software

- Python 3.10+
- R 4.x
- Git
- GitHub

## Python Ecosystem

Python was used for activities including:

- simulation;
- data manipulation;
- quality control;
- statistical evaluation;
- result consolidation; and
- visualization.

Common tools include:

- pandas
- NumPy
- SciPy

## R Ecosystem

R was used for the established DA packages:

- ALDEx2
- ANCOMBC / ANCOM-BC2
- DESeq2

Because analyses were developed across multiple team environments, exact package versions should be interpreted together with the scripts and environment information associated with each workflow.

---

# Reproducing the Repository

Clone the shared repository:

```bash
git clone https://github.com/Archie-ai/capstone-project-collab.git
cd capstone-project-collab
```

Verify Python:

```bash
python --version
```

Verify R:

```bash
Rscript --version
```

The collaborative repository contains multiple analysis workflows rather than one verified top-level executable.

Users should therefore navigate to the relevant contributor directory before reproducing a specific experiment.

For example, the team consolidation analysis includes:

```text
architha/scripts/team_consolidation/compare_null_vs_real.py
```

with outputs under:

```text
architha/results/team_consolidation/null_vs_real/
```

Before reproducing any analysis, confirm:

1. count-table and metadata sample IDs are correctly aligned;
2. filtering and prevalence criteria match the original analysis;
3. the same significance threshold is used;
4. method-specific preprocessing requirements are preserved;
5. software/package versions are documented; and
6. fixed random seeds are used for reproducible simulation.

---

# Canonical Null Simulation

The canonical null dataset was designed as an initial controlled validation dataset.

### Simulation Design

- **50 samples**
  - 25 Control
  - 25 Case
- **100 taxa**
- **0 truly differential taxa**
- Mean library size: **9,994.22**
- Observed library-size range: **9,753–10,191**
- Zero fraction: **4.68%**
- Random seed: **42**

Because there were **no true DA taxa**, every statistically significant result is known to be a false discovery.

## Canonical Null Results

| Method | Significant Taxa |
|---|---:|
| CLR/Wilcoxon | 1 |
| ALDEx2 | 0 |
| ANCOM-BC2 | 1 |
| DESeq2 | 0 |

CLR/Wilcoxon and ANCOM-BC2 both identified:

```text
Taxon_033
```

Because the simulation contains no true group differences, `Taxon_033` is a known false discovery in this particular replicate.

## Canonical Null Consensus

| Threshold | Consensus-Positive Taxa |
|---|---:|
| 2-of-4 | 1 |
| 3-of-4 | 0 |
| 4-of-4 | 0 |

This provides an important example:

> **Agreement among multiple statistical methods does not automatically establish biological or statistical truth.**

---

# Repeated-Null Benchmark

A single simulation cannot characterize expected method behavior.

Therefore, repeated null simulations were used to evaluate false-positive calibration.

The current 200-run benchmark produced:

| Method | Probability of ≥1 False Discovery | Mean False Positives per Run |
|---|---:|---:|
| CLR/Wilcoxon | **3.5%** | 0.035 |
| ALDEx2 | **3.0%** | 0.030 |
| ANCOM-BC | **9.5%** | 0.135 |
| DESeq2 | **14.0%** | 0.180 |

These results describe behavior under this specific simulation framework.

They should **not** be interpreted as a universal ranking of the four methods across all microbiome datasets.

---

# Repeated-Null Consensus Calibration

Consensus was also evaluated across repeated null simulations.

| Consensus Threshold | Mean Per-Taxon Type I Error | Replicates With ≥1 False Consensus |
|---|---:|---:|
| **2-of-4** | 0.080% | 7.0% |
| **3-of-4** | 0.030% | 3.0% |
| **4-of-4** | 0.020% | 2.0% |

Increasing the consensus requirement reduced false consensus.

However:

> **Even 4-of-4 agreement did not completely eliminate false-positive consensus.**

---

# Known-Signal Consensus Evaluation

Known-signal simulations contain deliberately planted DA taxa, allowing both sensitivity and false-discovery proportion to be evaluated.

| Analysis | Consensus | Sensitivity | Mean FDP |
|---|---|---:|---:|
| Erin — 5 replicates | 2-of-4 | 57.00% | 4.33% |
| | 3-of-4 | 50.00% | 0.00% |
| | 4-of-4 | 45.00% | 0.00% |
| Raquel — 200 replicates | 2-of-4 | 46.50% | 4.69% |
| | 3-of-4 | 42.55% | 1.02% |
| | 4-of-4 | 32.75% | 0.23% |
| Peggy — 5,400 replicates | 2-of-4 | 25.01% | 2.58% |
| | 3-of-4 | 20.33% | 0.47% |
| | 4-of-4 | 14.16% | 0.05% |

These results were produced under different simulation conditions, numbers of replicates, and consensus implementations.

They are therefore presented separately rather than pooled.

### Overall Pattern

Across the available known-signal analyses:

> **Stricter consensus reduced false discovery proportion but also reduced sensitivity.**

This illustrates the central consensus tradeoff:

```text
More stringent agreement
        ↓
Fewer false discoveries
        BUT
More missed true signals
```

---

# Real-Data Differential-Abundance Results

For harmonized cross-dataset comparison, **standard ANCOM significance** was used.

| Dataset | CLR/Wilcoxon | ALDEx2 | ANCOM Standard | DESeq2 |
|---|---:|---:|---:|---:|
| `cdi_schubert` | **96.76%** | 6.48% | 36.69% | 19.93% |
| `crc_baxter` | 2.12% | 0.00% | **42.89%** | 1.77% |
| `cdi_vincent` | 7.18% | 1.54% | **31.28%** | 17.44% |
| `par_scheperjans` | 0.145% | 0.00% | **52.82%** | 0.434% |

These results demonstrate substantial **method-by-dataset interaction**.

For example:

- CLR/Wilcoxon identified **96.76%** of taxa in `cdi_schubert`;
- the same approach identified only **0.145%** in `par_scheperjans`.

Likewise, standard ANCOM produced relatively large discovery sets across the four datasets.

These observations describe **method behavior**.

They do not establish which method was biologically correct because real-data ground truth is unknown.

---

# ANCOM Robustness Analysis

A notable finding was the difference between standard ANCOM significance and robust/sensitivity-supported findings.

| Dataset | Standard ANCOM | Robust ANCOM |
|---|---:|---:|
| `cdi_vincent` | 61 | 1 |
| `par_scheperjans` | 365 | 1 |

This substantial reduction indicates that standard ANCOM calls in these datasets were highly sensitive to the robustness criterion.

However:

> **Non-robust real-data taxa should not be described as false positives.**

Robustness evaluates statistical stability under the sensitivity procedure. It does not independently establish biological truth.

---

# Consensus Interpretation

The consensus framework evaluates whether taxa are supported by multiple analytical approaches.

### 2-of-4 Consensus

At least two methods identify the taxon.

- Most permissive consensus threshold
- Generally preserves more sensitivity
- May retain more false or unstable signals

### 3-of-4 Consensus

At least three methods identify the taxon.

- Intermediate threshold
- Stronger cross-method support
- Intermediate sensitivity/FDP tradeoff

### 4-of-4 Consensus

All four methods identify the taxon.

- Most stringent threshold
- Typically produces the smallest candidate set
- Can substantially reduce sensitivity

### Important Principle

Consensus should be interpreted as:

**cross-method reproducibility**

not:

**proof of biological correctness**

All four methods analyze the same underlying data and may therefore share biases or respond to the same dataset artifact.

---

# Main Findings

The combined analyses support the following conclusions:

1. **Different differential-abundance methods can produce substantially different results from the same microbiome dataset.**

2. **Method behavior depends strongly on dataset characteristics and statistical assumptions.**

3. **A single simulation replicate is insufficient to characterize method reliability.**

4. **Repeated null simulations revealed meaningful differences in false-positive calibration among methods.**

5. **Consensus reduced false discoveries as the agreement threshold became stricter.**

6. **Stricter consensus also reduced sensitivity to true signals.**

7. **Real-data discovery burden varied substantially across both methods and datasets.**

8. **Agreement among methods should be interpreted as reproducibility rather than biological truth.**

9. **No single DA method was consistently most favorable across all analyses examined.**

---

# Important Interpretation Rules

The following rules were used throughout the project to prevent statistical overinterpretation.

### Do not:

- label real-data discoveries as true positives or false positives without independent ground truth;
- call non-robust ANCOM findings false positives;
- compare effect-size magnitudes directly across methods;
- assume 4-of-4 consensus proves biological correctness;
- assume 3-of-4 consensus is universally optimal;
- pool simulation summaries generated under substantially different designs without qualification;
- interpret discovery count alone as method accuracy.

### Instead:

- distinguish simulated truth from real-data concordance;
- report the analysis universe;
- normalize discovery counts when comparing datasets;
- examine effect direction;
- document filtering and preprocessing;
- report method-specific limitations; and
- interpret consensus as evidence supporting prioritization for further validation.

---

# Project Limitations

Several limitations should be considered when interpreting the project.

### Simulation differences

Simulation designs and numbers of replicates differed across team contributors.

Therefore, not all simulation summaries can be directly pooled.

### Preprocessing differences

Real datasets required different filtering and preprocessing because of differences in:

- sample number;
- sparsity;
- prevalence;
- taxonomic structure; and
- available data.

### ALDEx2 implementation

Peggy's workflow used a Welch-primary ALDEx2 analysis, whereas other team workflows used Wilcoxon-primary ALDEx2.

### ANCOM robustness

Robustness/sensitivity reporting was not available in exactly the same form for every dataset.

### DESeq2

DESeq2 was developed primarily for RNA-seq rather than microbiome data and was included as a methodological contrast.

### Unknown real-data ground truth

The true set of biologically differential taxa in the real datasets is unknown.

### Shared statistical biases

Agreement among methods is not equivalent to independent confirmation because the methods analyze the same underlying data.

### Project scope

The study evaluates four DA approaches and four real datasets rather than every available microbiome DA method.

---

# Future Work

Future extensions could include:

- standardizing one common simulation design across all team members;
- increasing independent simulation replicates;
- evaluating additional real microbiome datasets;
- testing additional DA approaches;
- expanding sample-size, sparsity, prevalence, and effect-size scenario grids;
- developing weighted rather than simple voting consensus;
- predicting method reliability from measurable dataset properties;
- validating high-consensus taxa in independent biological cohorts; and
- integrating external biological or functional evidence.

---

# Reproducibility Principles

The project emphasizes reproducible bioinformatics practice through:

- version-controlled scripts;
- Git/GitHub collaboration;
- fixed random seeds for simulation;
- documented simulation parameters;
- explicit ground-truth files;
- sample-alignment validation;
- documented prevalence/filtering criteria;
- standardized result tables;
- consistent significance thresholds;
- separation of raw, intermediate, and final results where possible; and
- explicit distinction between statistical accuracy and real-data concordance.

---

# Key References

Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J. (2017). *Microbiome datasets are compositional: And this is not optional.* **Frontiers in Microbiology, 8**, 2224.

Hawinkel, S., Mattiello, F., Bijnens, L., & Thas, O. (2019). *A broken promise: Microbiome differential abundance methods do not control the false discovery rate.* **Briefings in Bioinformatics, 20**(1), 210–221.

Lin, H., & Peddada, S. D. (2020). *Analysis of compositions of microbiomes with bias correction.* **Nature Communications, 11**, 3514.

Love, M. I., Huber, W., & Anders, S. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* **Genome Biology, 15**, 550.

Nearing, J. T., Douglas, G. M., Hayes, M. G., et al. (2022). *Microbiome differential abundance methods produce different results across 38 datasets.* **Nature Communications, 13**, 342.

Quinn, T. P., Erb, I., Gloor, G., Notredame, C., Richardson, M. F., & Crowley, T. M. (2019). *A field guide for the compositional analysis of any-omics data.* **GigaScience, 8**(9), giz107.

---

# Team

**Architha Sankar**  
**Erin Maley**  
**Peggy Hall**  
**Raquel De Leon**

University of Maryland Global Campus  
BIOT 670I — Capstone Project

---

## Academic Use

This repository was developed as an academic capstone project examining statistical reliability and agreement in microbiome differential-abundance analysis.

Results should be interpreted within the simulation assumptions, preprocessing decisions, statistical methods, and limitations documented in the repository.
