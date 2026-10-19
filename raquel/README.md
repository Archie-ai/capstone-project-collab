# Microbiome Abundance Consensus Machine — Individual Contribution

**Contributor:** Raquel De Leon  
**Course:** BIOT 670I — Biotechnology Capstone  
**Project:** Microbiome Abundance Consensus Machine

## Overview

This folder contains my individual contributions to the Microbiome Abundance Consensus Machine capstone project. The project evaluates how different statistical methods identify differentially abundant microbial taxa and whether combining their results can improve consistency and reliability.

My work focuses on implementing differential abundance (DA) methods, evaluating their performance using simulated datasets with known ground truth, applying the methods to a real microbiome dataset, and comparing their findings through consensus analysis.

## Statistical Methods

Four differential abundance methods were evaluated:

1. **CLR + Wilcoxon Rank-Sum:** Uses centered log-ratio transformation followed by a nonparametric comparison between groups.
2. **ALDEx2:** Accounts for compositionality and uncertainty in microbiome count data.
3. **ANCOM-BC2:** Applies bias correction to estimate differences in microbial abundance.
4. **DESeq2:** Uses a negative binomial statistical model to analyze count data.

## Simulation Analysis

The known-DA evaluation included **200 simulated datasets**, each containing:

- 50 samples (25 Control and 25 Case)
- 100 microbial taxa
- 10 truly differentially abundant taxa
- A simulated effect size of 2.0

The four methods were evaluated using false discovery proportion (FDP), statistical power, sensitivity, and specificity.

Consensus strategies were also evaluated by requiring agreement among at least two, three, or all four methods. These analyses examine the trade-off between identifying true differences and limiting false discoveries.

## Real Microbiome Dataset

The methods were applied to the **cdi_vincent** dataset from MicrobiomeHD, which contains samples associated with *Clostridioides difficile* infection (CDI).

The analysis included:

- 50 samples (25 Case and 25 Control)
- 618 original microbial features
- 195 taxa retained after prevalence filtering

Results from all four methods were compared to evaluate differences in significant discoveries and identify taxa supported by multiple approaches.

Because the real dataset does not provide known ground truth, agreement among methods was evaluated without interpreting consensus as proof of biological significance.

## Repository Organization

| Folder | Contents |
|---|---|
| `python/` | Python scripts for simulation, CLR + Wilcoxon analysis, evaluation, and consensus |
| `R/` | R scripts for ALDEx2, ANCOM-BC2, and DESeq2 |
| `data/` | Simulated datasets and processed CDI microbiome data |
| `results/da/` | Performance results from known-DA simulations |
| `results/real/cdi_vincent/` | Real-data DA results, consensus comparisons, and validation analyses |
| `results/figures/` | Seven figures provided in PNG and PDF formats |
| `results/null/` | Initial null simulation results |

## Key Findings

The simulation results showed that individual methods differed in their ability to identify true DA taxa while controlling false discoveries.

Consensus approaches generally reduced false discoveries as the agreement requirement increased, but stricter agreement also reduced statistical power.

The real CDI analysis demonstrated substantial differences in the number of significant taxa identified by each method. Consensus analysis provided a way to identify discoveries supported by multiple statistical approaches.

Overall, these findings highlight the importance of considering both statistical performance and agreement when interpreting microbiome differential abundance results.

## Notes

This directory documents my individual analysis contributions and is intended to support integration into the group's final capstone report and presentation.
