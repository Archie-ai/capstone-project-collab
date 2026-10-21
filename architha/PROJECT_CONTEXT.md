# Microbiome Abundance Consensus Machine

## 1. Project Title

**Microbiome Abundance Consensus Machine: Benchmarking Differential-Abundance Methods Across Simulated and Real Microbiome Data**

---

## 2. Central Research Problem

Different microbiome differential-abundance (DA) methods can identify substantially different sets of significant taxa when applied to the same dataset. This creates uncertainty when interpreting biological results because agreement between methods is often incomplete and method performance can depend on characteristics of the underlying microbiome data.

This project will systematically investigate when selected DA methods agree, disagree, or fail and whether these behaviors can be related to measurable properties of the dataset.

---

## 3. Primary Research Question

**How do microbiome data characteristics influence the statistical performance and agreement of differential-abundance methods, and can a consensus strategy provide more robust identification of differentially abundant taxa?**

---

## 4. Primary Hypothesis

The performance of microbiome differential-abundance methods will vary systematically with dataset characteristics such as sample size, sparsity, prevalence, and effect size.

Methods are therefore expected to show different trade-offs between statistical power and false-discovery control under different experimental conditions.

A consensus strategy that incorporates agreement among methods may provide more robust results than relying on a single DA method, although this must be validated against known simulated ground truth.

---

## 5. Core Differential-Abundance Methods

The initial benchmark will include four complementary approaches:

1. **CLR/Wilcoxon**

   * Compositional transformation using the centered log-ratio approach.
   * Statistical comparison using the Wilcoxon rank-sum test.
   * Serves as a relatively simple baseline.

2. **ALDEx2**

   * Designed for compositional data.
   * Uses Monte Carlo sampling from Dirichlet distributions and log-ratio transformations.

3. **ANCOM-BC**

   * Designed specifically for microbiome compositional data.
   * Includes bias correction for compositional effects.

4. **DESeq2**

   * Originally developed for RNA-seq count data.
   * Included intentionally because it is frequently applied to microbiome data despite assumptions that may not always be optimal for highly sparse compositional datasets.

---

## 6. Experimental Design

The project contains two complementary analytical arms.

### Arm A — Simulated Microbiome Data

Simulated microbiome count tables will be generated with known ground truth.

The simulator will define in advance:

* Which taxa are truly differentially abundant.
* Which taxa are non-differential.
* Direction of abundance change.
* Magnitude of effect.
* Sample size.
* Sparsity.
* Taxon prevalence.

Because ground truth is explicitly known, each DA method can be evaluated objectively.

### Arm B — Real Microbiome Data

A selected public microbiome dataset will later be analyzed using the same DA methods.

Potential sources include:

* microbiomeHD
* curatedMetagenomicData

Because true DA status is not known for real biological data, evaluation will focus on agreement and disagreement among methods rather than direct classification accuracy.

---

## 7. Definition of Ground Truth

For simulation experiments, ground truth is defined before the DA methods are run.

Each taxon will have a known status:

* **DA = 1** → deliberately planted as differentially abundant.
* **DA = 0** → not deliberately altered.

For DA taxa, the simulator will also record:

* Expected direction of effect.
* Effect magnitude.

The DA algorithms will not receive this truth information.

Ground truth will only be used during downstream evaluation.

This separation prevents circular evaluation.

---

## 8. Initial Null-Control Experiment

The first computational experiment will contain:

**Zero truly differentially abundant taxa.**

Therefore:

```text
Number of planted DA taxa = 0
```

The purpose of this experiment is to test whether the pipeline produces false discoveries when no biological difference exists.

At a nominal significance threshold of:

```text
alpha = 0.05
```

false-positive behavior should remain reasonably close to the expected nominal level.

If a method reports excessive significant taxa under the null condition, the team must investigate:

* Wrapper implementation.
* Input formatting.
* Statistical testing.
* Multiple-testing correction.
* Filtering.
* Normalization/transformation.
* Simulation assumptions.

Downstream conclusions should not be trusted until the null workflow behaves appropriately.

---

## 9. Initial Toy Dataset

Development will begin with a small dataset before scaling to larger experiments.

Initial design:

```text
Total samples: approximately 50
Groups: 2
Group 1: Control
Group 2: Case
Initial taxa: approximately 100
```

For the first null experiment:

```text
Planted DA taxa: 0
```

The small dataset is intended for debugging and validation rather than final biological conclusions.

---

## 10. Experimental Variables

After null calibration succeeds, the simulator will systematically vary the following properties.

### Sample Size

Potential levels:

```text
20
50
100
200
```

### Sparsity

Potential levels:

```text
Low
Moderate
High
```

Sparsity represents the frequency of zero counts within the abundance table.

### Effect Size

Potential levels:

```text
Weak
Moderate
Strong
```

### Taxon Prevalence

Taxon prevalence represents the proportion of samples in which a given taxon is observed.

Rare and common taxa may behave differently across DA methods.

### Combined Scenarios

After single-variable experiments are validated, combinations will be evaluated, for example:

```text
Small sample size + high sparsity + weak effect

Large sample size + high sparsity + strong effect

Small sample size + low sparsity + strong effect

Large sample size + low sparsity + weak effect
```

---

## 11. Primary Evaluation Metrics

For simulated datasets:

### True Positive

A truly DA taxon correctly identified as significant.

### False Positive

A non-DA taxon incorrectly identified as significant.

### False Negative

A truly DA taxon not identified as significant.

### False Discovery Rate

```text
FDR = FP / (TP + FP)
```

when at least one taxon is called significant.

### Power / Sensitivity

```text
Power = TP / (TP + FN)
```

### Specificity

```text
Specificity = TN / (TN + FP)
```

### Additional Metrics

Where appropriate:

* ROC curve
* ROC-AUC
* Precision
* Recall
* Effect-direction accuracy
* Runtime
* Number of method failures

---

## 12. Real-Data Evaluation

Because the absolute biological truth is unknown for real datasets, evaluation will focus on method concordance.

Potential measures include:

* Number of significant taxa per method.
* Shared significant taxa.
* Pairwise overlap.
* Jaccard similarity.
* Direction-of-effect agreement.
* Correlation between estimated effect sizes.
* Taxa detected uniquely by individual methods.
* Taxa consistently detected across methods.

---

## 13. Consensus Strategy

A consensus method will not be selected in advance.

Candidate strategies may include:

### Intersection

A taxon is retained only when all methods identify it.

### Majority Vote

A taxon is retained when a predefined majority of methods identify it.

Example:

```text
3 of 4 methods
```

### Weighted Consensus

Methods may receive different weights based on their demonstrated reliability under particular simulation conditions.

A conceptual score is:

```text
Consensus score = Σ(weight_method × method_support)
```

The final consensus rule must be evaluated against simulated ground truth.

Agreement alone will not be treated as proof of biological correctness.

---

## 14. Expected Scientific Outcome

The intended outcome is not simply to identify one universally superior DA method.

The stronger objective is to develop a relationship of the form:

```text
DATA CHARACTERISTIC
        ↓
METHOD BEHAVIOR
        ↓
EXPECTED RELIABILITY
```

For example, a method may show:

* Strong FDR control but reduced power under high sparsity.
* High sensitivity but excessive false positives in small samples.
* Improved stability at larger sample sizes.
* Strong performance only when effect sizes are sufficiently large.

The final project should therefore describe method behavior conditionally rather than declaring a universal winner.

---

## 15. Reproducibility Requirements

The project repository should ultimately contain:

* Pinned Python environment.
* Documented R environment.
* Simulation scripts.
* Configuration files.
* Figures.
---

## 16. Immediate Next Experimental Goal
```text
50 samples
100 taxa
2 groups
0 planted DA taxa
```

This null dataset will be used to verify the integrity of the pipeline before introducing true biological effects.

The next objective is to build the first simulated microbiome count table containing:

* Result tables.
* Documentation.

* Final report.

All analysis steps should be reproducible from version-controlled scripts.
* DA method wrappers.
* Evaluation scripts.
* Reproducible workflow driver.

