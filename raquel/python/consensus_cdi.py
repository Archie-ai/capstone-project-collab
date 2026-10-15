# BIOT 670I - Microbiome Abundance Consensus Machine
# Real-data method concordance and consensus analysis
# Dataset: cdi_vincent

from pathlib import Path
from itertools import combinations

import pandas as pd

# --------------------------------------------------
# 1. File locations
# --------------------------------------------------

BASE = Path("results/real/cdi_vincent")
DATA = Path("data/processed/cdi_vincent")
OUTPUT = BASE / "consensus"

OUTPUT.mkdir(parents=True, exist_ok=True)

METHOD_FILES = {
    "CLR_Wilcoxon": BASE / "clr_wilcoxon/clr_wilcoxon_results.csv",
    "ALDEx2": BASE / "aldex2/aldex2_results.csv",
    "ANCOMBC2": BASE / "ancombc/ancombc_results.csv",
    "DESeq2": BASE / "deseq2/deseq2_results.csv",
}

METHODS = [
    "CLR_Wilcoxon",
    "ALDEx2",
    "ANCOMBC2",
    "DESeq2",
]

# --------------------------------------------------
# 2. Load reference OTU list and taxonomy
# --------------------------------------------------

counts = pd.read_csv(
    DATA / "cdi_counts_filtered.csv",
    index_col=0
)

taxonomy = pd.read_csv(
    DATA / "cdi_taxonomy_mapping.csv"
)

if counts.index.duplicated().any():
    raise ValueError("Duplicate OTU IDs in count table")

if taxonomy["OTU_ID"].duplicated().any():
    raise ValueError("Duplicate OTU IDs in taxonomy mapping")

reference_taxa = set(counts.index)

# --------------------------------------------------
# 3. Load method results
# --------------------------------------------------

method_results = {}

for method, path in METHOD_FILES.items():

    df = pd.read_csv(path)

    if df["Taxon"].duplicated().any():
        raise ValueError(f"Duplicate OTUs in {path}")

    unexpected = set(df["Taxon"]) - reference_taxa

    if unexpected:
        raise ValueError(
            f"Unexpected OTUs in {method}: {len(unexpected)}"
        )

    method_results[method] = df.set_index("Taxon")

# --------------------------------------------------
# 4. Analyze both ANCOM-BC2 significance criteria
# --------------------------------------------------

SCENARIOS = {
    "robust": "Robust_Significant",
    "adjusted": "Significant",
}

for scenario, ancombc_column in SCENARIOS.items():

    print(f"\nAnalyzing ANCOM-BC2 scenario: {scenario}")

    # Start with all 195 filtered OTUs
    votes = pd.DataFrame({
        "Taxon": counts.index
    })

    # Add significance calls from each method
    for method in METHODS:

        df = method_results[method]

        sig_column = (
            ancombc_column
            if method == "ANCOMBC2"
            else "Significant"
        )

        if sig_column not in df.columns:
            raise ValueError(
                f"Missing {sig_column} in {method}"
            )

        # Missing/unreported calls count as not detected
        calls = df[sig_column].reindex(votes["Taxon"])

        if calls.notna().any():
            valid_values = set(
                calls.dropna().astype(str).str.lower()
            )

            if not valid_values.issubset({"true", "false"}):
                raise ValueError(
                    f"Unexpected significance values in {method}"
                )

        votes[method] = (
            calls.fillna(False)
            .astype(bool)
            .to_numpy()
        )

    # Count the number of methods detecting each OTU
    votes["Votes"] = votes[METHODS].sum(axis=1)

    # Consensus voting rules
    votes["Consensus_2of4"] = votes["Votes"] >= 2
    votes["Consensus_3of4"] = votes["Votes"] >= 3
    votes["Consensus_4of4"] = votes["Votes"] == 4

    # Attach original taxonomy annotations
    votes = votes.merge(
        taxonomy[["OTU_ID", "Taxonomy"]],
        left_on="Taxon",
        right_on="OTU_ID",
        how="left",
        validate="one_to_one"
    ).drop(columns="OTU_ID")

    if votes["Taxonomy"].isna().any():
        raise ValueError("Missing taxonomy annotations")

    # Save per-OTU method calls and consensus votes
    votes.to_csv(
        OUTPUT / f"otu_votes_{scenario}.csv",
        index=False
    )

    # --------------------------------------------------
    # 5. Method discovery counts
    # --------------------------------------------------

    discovery_counts = pd.DataFrame([
        {
            "Method": method,
            "Significant_OTUs": int(votes[method].sum())
        }
        for method in METHODS
    ])

    discovery_counts.to_csv(
        OUTPUT / f"method_discoveries_{scenario}.csv",
        index=False
    )

    # --------------------------------------------------
    # 6. Pairwise overlap and Jaccard similarity
    # --------------------------------------------------

    pairwise_rows = []

    for method_a, method_b in combinations(METHODS, 2):

        set_a = set(
            votes.loc[votes[method_a], "Taxon"]
        )

        set_b = set(
            votes.loc[votes[method_b], "Taxon"]
        )

        intersection = len(set_a & set_b)
        union = len(set_a | set_b)

        jaccard = (
            intersection / union
            if union > 0
            else 1.0
        )

        pairwise_rows.append({
            "Method_A": method_a,
            "Method_B": method_b,
            "Significant_A": len(set_a),
            "Significant_B": len(set_b),
            "Shared_OTUs": intersection,
            "Union_OTUs": union,
            "Jaccard": jaccard,
        })

    pairwise = pd.DataFrame(pairwise_rows)

    pairwise.to_csv(
        OUTPUT / f"pairwise_agreement_{scenario}.csv",
        index=False
    )

    # --------------------------------------------------
    # 7. Vote distribution
    # --------------------------------------------------

    vote_distribution = pd.DataFrame({
        "Votes": range(5)
    })

    vote_counts = votes["Votes"].value_counts()

    vote_distribution["OTUs"] = [
        int(vote_counts.get(i, 0))
        for i in range(5)
    ]

    vote_distribution.to_csv(
        OUTPUT / f"vote_distribution_{scenario}.csv",
        index=False
    )

    # --------------------------------------------------
    # 8. Consensus summary
    # --------------------------------------------------

    consensus_summary = pd.DataFrame([
        {
            "Rule": "2-of-4",
            "Significant_OTUs": int(
                votes["Consensus_2of4"].sum()
            )
        },
        {
            "Rule": "3-of-4",
            "Significant_OTUs": int(
                votes["Consensus_3of4"].sum()
            )
        },
        {
            "Rule": "4-of-4",
            "Significant_OTUs": int(
                votes["Consensus_4of4"].sum()
            )
        },
    ])

    consensus_summary.to_csv(
        OUTPUT / f"consensus_summary_{scenario}.csv",
        index=False
    )

    # --------------------------------------------------
    # 9. Print results
    # --------------------------------------------------

    print("\nMethod discoveries:")
    print(discovery_counts.to_string(index=False))

    print("\nConsensus results:")
    print(consensus_summary.to_string(index=False))

    print("\nPairwise agreement:")
    print(
        pairwise[
            ["Method_A", "Method_B", "Shared_OTUs", "Jaccard"]
        ].to_string(index=False)
    )

print("\nReal-data consensus analysis completed!")
print(f"Results saved to: {OUTPUT}")
