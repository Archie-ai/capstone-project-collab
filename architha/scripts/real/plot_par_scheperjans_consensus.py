#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


INPUT_DIR = Path(
    "results/real/par_scheperjans/consolidated"
)

OUTPUT_DIR = Path(
    "results/real/par_scheperjans/figures"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 1. JACCARD HEATMAP
# ============================================================

jaccard = pd.read_csv(
    INPUT_DIR /
    "real_significant_jaccard_matrix.csv",
    index_col=0
)


display_names = {
    "CLR_Wilcoxon":
        "CLR/\nWilcoxon",

    "ALDEx2":
        "ALDEx2",

    "ANCOM_BC2_Robust":
        "ANCOM-BC2\nRobust",

    "DESeq2":
        "DESeq2",
}


jaccard = jaccard.rename(
    index=display_names,
    columns=display_names
)


fig, ax = plt.subplots(
    figsize=(7, 6)
)


image = ax.imshow(
    jaccard.to_numpy(
        dtype=float
    ),
    vmin=0,
    vmax=1
)


ax.set_xticks(
    np.arange(
        len(jaccard.columns)
    )
)

ax.set_yticks(
    np.arange(
        len(jaccard.index)
    )
)


ax.set_xticklabels(
    jaccard.columns
)

ax.set_yticklabels(
    jaccard.index
)


for i in range(
    len(jaccard.index)
):

    for j in range(
        len(jaccard.columns)
    ):

        value = jaccard.iloc[
            i,
            j
        ]

        if pd.isna(value):

            label = "NA"

        else:

            label = f"{value:.2f}"

        ax.text(
            j,
            i,
            label,
            ha="center",
            va="center"
        )


ax.set_title(
    "Pairwise Jaccard Similarity of Significant Taxa\n"
    "par_scheperjans"
)


fig.colorbar(
    image,
    ax=ax,
    label="Jaccard similarity"
)


fig.tight_layout()


fig.savefig(
    OUTPUT_DIR /
    "par_scheperjans_jaccard_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)


plt.close(fig)


# ============================================================
# 2. METHOD DISCOVERY COUNTS
# ============================================================

method_summary = pd.read_csv(
    INPUT_DIR /
    "real_method_summary.csv"
)


method_order = [
    "CLR_Wilcoxon",
    "ALDEx2",
    "ANCOM_BC2_Robust",
    "DESeq2",
]


method_labels = [
    "CLR/Wilcoxon",
    "ALDEx2",
    "ANCOM-BC2\nRobust",
    "DESeq2",
]


method_summary = (
    method_summary
    .set_index("Method")
    .loc[method_order]
)


counts = (
    method_summary[
        "Primary_Significant"
    ]
    .astype(int)
    .to_numpy()
)


fig, ax = plt.subplots(
    figsize=(7, 5)
)


bars = ax.bar(
    method_labels,
    counts
)


ax.set_ylabel(
    "Number of significant taxa"
)


ax.set_title(
    "Primary Differential-Abundance Discoveries\n"
    "par_scheperjans"
)


ax.set_ylim(
    0,
    max(counts) + 1
)


for bar, value in zip(
    bars,
    counts
):

    ax.text(
        bar.get_x()
        +
        bar.get_width() / 2,
        value + 0.05,
        str(value),
        ha="center",
        va="bottom"
    )


fig.tight_layout()


fig.savefig(
    OUTPUT_DIR /
    "par_scheperjans_method_discovery_counts.png",
    dpi=300,
    bbox_inches="tight"
)


plt.close(fig)


# ============================================================
# 3. CONSENSUS DISTRIBUTION
# ============================================================

consensus = pd.read_csv(
    INPUT_DIR /
    "real_consensus_all_taxa.csv"
)


distribution = (
    consensus[
        "Consensus_Count"
    ]
    .value_counts()
    .reindex(
        [0, 1, 2, 3, 4],
        fill_value=0
    )
)


labels = [
    "0/4",
    "1/4",
    "2/4",
    "3/4",
    "4/4",
]


values = distribution.to_numpy()


fig, ax = plt.subplots(
    figsize=(7, 5)
)


bars = ax.bar(
    labels,
    values
)


ax.set_xlabel(
    "Number of methods calling taxon significant"
)

ax.set_ylabel(
    "Number of taxa"
)

ax.set_title(
    "Four-Method Consensus Distribution\n"
    "par_scheperjans"
)


for bar, value in zip(
    bars,
    values
):

    ax.text(
        bar.get_x()
        +
        bar.get_width() / 2,
        value,
        str(value),
        ha="center",
        va="bottom"
    )


fig.tight_layout()


fig.savefig(
    OUTPUT_DIR /
    "par_scheperjans_consensus_distribution.png",
    dpi=300,
    bbox_inches="tight"
)


plt.close(fig)


print(
    "=" * 72
)

print(
    "PAR_SCHEPERJANS FIGURE GENERATION"
)

print(
    "=" * 72
)

print(
    "\nCreated:"
)

print(
    OUTPUT_DIR /
    "par_scheperjans_jaccard_heatmap.png"
)

print(
    OUTPUT_DIR /
    "par_scheperjans_method_discovery_counts.png"
)

print(
    OUTPUT_DIR /
    "par_scheperjans_consensus_distribution.png"
)

print(
    "\nFIGURE GENERATION: PASS"
)
