## Current Progress

- Defined the primary research question and experimental framework.
- Established simulated-data and real-data analysis arms.
- Selected four core differential-abundance methods:
  - CLR/Wilcoxon
  - ALDEx2
  - ANCOM-BC
  - DESeq2
- Created the initial null microbiome simulation framework.
- Generated a 100-taxa × 50-sample null dataset.
- Balanced experimental groups:
  - 25 Control
  - 25 Case
- Confirmed zero truly differential taxa.
- Observed mean sequencing depth: 9,994.22 reads/sample.
- Observed zero fraction: 4.68%.
- Classified the initial simulation as a low-sparsity null baseline.
- Next step: formal simulation QC followed by implementation of the CLR/Wilcoxon baseline.

## Central Workflow
Simulated microbiome data
        ↓
Known ground truth
        ↓
Run multiple DA methods
        ↓
Standardize results
        ↓
Evaluate FDR / Power / ROC
        ↓
Determine when methods disagree
        ↓
Relate data properties to method behavior
        ↓
Evaluate consensus strategies
