#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd


EXPECTED_TAXA = 691


FILES = {
    "CLR_Wilcoxon": Path(
        "results/real/par_scheperjans/"
        "clr_wilcoxon/clr_wilcoxon_results.csv"
    ),

    "ALDEx2": Path(
        "results/real/par_scheperjans/"
        "aldex2/aldex2_results.csv"
    ),

    "ANCOM_BC2_Robust": Path(
        "results/real/par_scheperjans/"
        "ancombc2/ancombc2_results.csv"
    ),

    "DESeq2": Path(
        "results/real/par_scheperjans/"
        "deseq2/deseq2_results.csv"
    ),
}


SIG_COL = {
    "CLR_Wilcoxon": "Significant",
    "ALDEx2": "Significant",
    "ANCOM_BC2_Robust": "Robust_Significant",
    "DESeq2": "Significant",
}


PREFIX = {
    "CLR_Wilcoxon": "CLR",
    "ALDEx2": "ALDEx2",
    "ANCOM_BC2_Robust": "ANCOM",
    "DESeq2": "DESeq2",
}


OUT = Path(
    "results/real/par_scheperjans/consolidated"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# HELPERS
# ============================================================

def as_bool(series):

    if pd.api.types.is_bool_dtype(series):

        return (
            series
            .fillna(False)
            .astype(bool)
        )

    text = (
        series
        .astype(str)
        .str.strip()
        .str.lower()
    )

    true_values = {
        "true",
        "t",
        "1",
        "yes",
        "y",
    }

    false_values = {
        "false",
        "f",
        "0",
        "no",
        "n",
        "",
        "nan",
        "none",
    }

    unexpected = (
        set(text.unique())
        -
        true_values
        -
        false_values
    )

    if unexpected:

        raise ValueError(
            "Unexpected boolean values "
            f"in {series.name}: "
            f"{sorted(unexpected)}"
        )

    return text.isin(
        true_values
    )


def require(
    df,
    columns,
    label
):

    missing = [
        column
        for column in columns
        if column not in df.columns
    ]

    if missing:

        raise ValueError(
            f"{label} missing columns: "
            f"{missing}"
        )


def primary(
    df,
    method
):

    sig_col = SIG_COL[
        method
    ]

    require(
        df,
        [
            "Taxon",
            "P_Value",
            "Adjusted_P_Value",
            "Effect_Estimate",
            "Effect_Direction",
            sig_col,
        ],
        method,
    )

    if "Tested" in df.columns:

        tested = as_bool(
            df["Tested"]
        )

    else:

        tested = (
            pd.to_numeric(
                df["P_Value"],
                errors="coerce"
            )
            .notna()
        )


    output = pd.DataFrame(
        {
            "Taxon":
                df["Taxon"]
                .astype(str),

            "Method":
                method,

            "Tested":
                tested,

            "P_Value":
                pd.to_numeric(
                    df["P_Value"],
                    errors="coerce"
                ),

            "Adjusted_P_Value":
                pd.to_numeric(
                    df[
                        "Adjusted_P_Value"
                    ],
                    errors="coerce"
                ),

            "Significant":
                as_bool(
                    df[
                        sig_col
                    ]
                ),

            "Effect_Estimate":
                pd.to_numeric(
                    df[
                        "Effect_Estimate"
                    ],
                    errors="coerce"
                ),

            "Effect_Direction":
                df[
                    "Effect_Direction"
                ]
                .astype(str),
        }
    )


    for column in [
        "Denovo_ID",
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "Source_Taxon",
    ]:

        if column in df.columns:

            output[
                column
            ] = df[
                column
            ]


    return output


# ============================================================
# LOAD
# ============================================================

print(
    "=" * 72
)

print(
    "PAR_SCHEPERJANS FOUR-METHOD "
    "REAL-DATA CONSOLIDATION"
)

print(
    "=" * 72
)


raw = {
    name: pd.read_csv(
        path
    )
    for name, path
    in FILES.items()
}


tables = {
    name: primary(
        raw[name],
        name
    )
    for name
    in FILES
}


print(
    "\nInput dimensions:"
)


for name in FILES:

    print(
        f"{name:18s} "
        f"{raw[name].shape}"
    )


# ============================================================
# TAXON UNIVERSE
# ============================================================

reference = set(
    tables[
        "CLR_Wilcoxon"
    ]["Taxon"]
)


print(
    "\nTaxon-universe checks:"
)


for name, df in tables.items():

    taxa = set(
        df["Taxon"]
    )

    missing = (
        reference -
        taxa
    )

    extra = (
        taxa -
        reference
    )

    print(
        f"{name:18s} "
        f"taxa={len(taxa)} "
        f"missing={len(missing)} "
        f"extra={len(extra)}"
    )

    if taxa != reference:

        raise ValueError(
            "Taxon-universe mismatch "
            f"for {name}"
        )


if len(reference) != EXPECTED_TAXA:

    raise ValueError(
        f"Expected {EXPECTED_TAXA} taxa; "
        f"found {len(reference)}"
    )


# ============================================================
# LONG TABLE
# ============================================================

long_df = pd.concat(
    tables.values(),
    ignore_index=True,
    sort=False
)


long_df.to_csv(
    OUT /
    "real_four_method_long.csv",
    index=False
)


# ============================================================
# METHOD SUMMARY
# ============================================================

summary_rows = []


for name, df in tables.items():

    row = {
        "Method":
            name,

        "Taxa_Available":
            len(df),

        "Taxa_Tested":
            int(
                df[
                    "Tested"
                ].sum()
            ),

        "Primary_Significant":
            int(
                df[
                    "Significant"
                ].sum()
            ),

        "Minimum_P_Value":
            df[
                "P_Value"
            ].min(
                skipna=True
            ),

        "Minimum_Adjusted_P_Value":
            df[
                "Adjusted_P_Value"
            ].min(
                skipna=True
            ),
    }


    if (
        name == "ALDEx2"
        and
        "Welch_Significant"
        in raw[name].columns
    ):

        row[
            "Diagnostic_Welch_Significant"
        ] = int(
            as_bool(
                raw[name][
                    "Welch_Significant"
                ]
            ).sum()
        )


    if name == "ANCOM_BC2_Robust":

        if (
            "Standard_Significant"
            in raw[name].columns
        ):

            row[
                "ANCOM_Standard_Significant"
            ] = int(
                as_bool(
                    raw[name][
                        "Standard_Significant"
                    ]
                ).sum()
            )


        if (
            "Passed_Pseudocount_Sensitivity"
            in raw[name].columns
        ):

            row[
                "ANCOM_Passed_Sensitivity"
            ] = int(
                as_bool(
                    raw[name][
                        "Passed_Pseudocount_Sensitivity"
                    ]
                ).sum()
            )


    if name == "DESeq2":

        if (
            "NoIF_Significant"
            in raw[name].columns
        ):

            row[
                "DESeq2_NoIF_Significant"
            ] = int(
                as_bool(
                    raw[name][
                        "NoIF_Significant"
                    ]
                ).sum()
            )


        if (
            "FullDiagnostic_Significant"
            in raw[name].columns
        ):

            row[
                "DESeq2_NoFilter_NoCooks_Significant"
            ] = int(
                as_bool(
                    raw[name][
                        "FullDiagnostic_Significant"
                    ]
                ).sum()
            )


        if (
            "Dispersion_Fit"
            in raw[name].columns
        ):

            fits = (
                raw[name][
                    "Dispersion_Fit"
                ]
                .dropna()
                .astype(str)
                .unique()
            )

            if len(fits) == 1:

                row[
                    "DESeq2_Dispersion_Fit"
                ] = fits[0]

            else:

                row[
                    "DESeq2_Dispersion_Fit"
                ] = ";".join(
                    fits
                )


    summary_rows.append(
        row
    )


method_summary = pd.DataFrame(
    summary_rows
)


method_summary.to_csv(
    OUT /
    "real_method_summary.csv",
    index=False
)


# ============================================================
# TAXONOMY BASE
# ============================================================

taxonomy_columns = [
    column
    for column in [
        "Taxon",
        "Denovo_ID",
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "Source_Taxon",
    ]
    if column
    in raw[
        "CLR_Wilcoxon"
    ].columns
]


wide = (
    raw[
        "CLR_Wilcoxon"
    ][
        taxonomy_columns
    ]
    .drop_duplicates(
        subset="Taxon"
    )
    .copy()
)


if len(wide) != EXPECTED_TAXA:

    raise ValueError(
        "Taxonomy base does not "
        "contain exactly 691 taxa"
    )


# ============================================================
# MERGE FOUR METHODS
# ============================================================

for name, df in tables.items():

    prefix = PREFIX[
        name
    ]


    add = df[
        [
            "Taxon",
            "Significant",
            "P_Value",
            "Adjusted_P_Value",
            "Effect_Estimate",
            "Effect_Direction",
        ]
    ].copy()


    add.columns = [
        "Taxon",
        f"{prefix}_Significant",
        f"{prefix}_P_Value",
        f"{prefix}_Adjusted_P_Value",
        f"{prefix}_Effect_Estimate",
        f"{prefix}_Effect_Direction",
    ]


    wide = wide.merge(
        add,
        on="Taxon",
        how="left",
        validate="one_to_one"
    )


sig_columns = [
    "CLR_Significant",
    "ALDEx2_Significant",
    "ANCOM_Significant",
    "DESeq2_Significant",
]


for column in sig_columns:

    wide[
        column
    ] = as_bool(
        wide[
            column
        ]
    )


# ============================================================
# CONSENSUS COUNTS
# ============================================================

wide[
    "Consensus_Count"
] = (
    wide[
        sig_columns
    ]
    .astype(int)
    .sum(axis=1)
)


wide[
    "Consensus_2of4"
] = (
    wide[
        "Consensus_Count"
    ] >= 2
)


wide[
    "Consensus_3of4"
] = (
    wide[
        "Consensus_Count"
    ] >= 3
)


wide[
    "Consensus_4of4"
] = (
    wide[
        "Consensus_Count"
    ] == 4
)


wide[
    "Consensus_Class"
] = (
    wide[
        "Consensus_Count"
    ]
    .map(
        {
            0: "0_of_4",
            1: "1_of_4",
            2: "2_of_4",
            3: "3_of_4",
            4: "4_of_4",
        }
    )
)


# ============================================================
# DIRECTION-AWARE CONSENSUS
# ============================================================

method_directions = [
    (
        "CLR_Significant",
        "CLR_Effect_Direction"
    ),

    (
        "ALDEx2_Significant",
        "ALDEx2_Effect_Direction"
    ),

    (
        "ANCOM_Significant",
        "ANCOM_Effect_Direction"
    ),

    (
        "DESeq2_Significant",
        "DESeq2_Effect_Direction"
    ),
]


consensus_direction = []

direction_agreement = []

direction_2of4 = []


for _, row in wide.iterrows():

    directions = []


    for sig_col, direction_col in method_directions:

        if bool(
            row[
                sig_col
            ]
        ):

            direction = str(
                row[
                    direction_col
                ]
            )

            if direction not in {
                "nan",
                "",
                "No_difference",
            }:

                directions.append(
                    direction
                )


    unique_directions = set(
        directions
    )


    if not directions:

        consensus_direction.append(
            "No_primary_discovery"
        )

        direction_agreement.append(
            np.nan
        )

        direction_2of4.append(
            False
        )


    elif len(
        unique_directions
    ) == 1:

        consensus_direction.append(
            directions[0]
        )

        direction_agreement.append(
            True
        )

        direction_2of4.append(
            int(
                row[
                    "Consensus_Count"
                ]
            ) >= 2
        )


    else:

        consensus_direction.append(
            "Discordant"
        )

        direction_agreement.append(
            False
        )

        direction_2of4.append(
            False
        )


wide[
    "Consensus_Direction"
] = consensus_direction


wide[
    "Significant_Direction_Agreement"
] = direction_agreement


wide[
    "Direction_Aware_2of4"
] = direction_2of4


# ============================================================
# SAVE CONSENSUS TABLES
# ============================================================

wide.to_csv(
    OUT /
    "real_consensus_all_taxa.csv",
    index=False
)


wide.loc[
    wide[
        "Consensus_2of4"
    ]
].to_csv(
    OUT /
    "real_consensus_2of4.csv",
    index=False
)


wide.loc[
    wide[
        "Consensus_3of4"
    ]
].to_csv(
    OUT /
    "real_consensus_3of4.csv",
    index=False
)


wide.loc[
    wide[
        "Consensus_4of4"
    ]
].to_csv(
    OUT /
    "real_consensus_4of4.csv",
    index=False
)


wide.loc[
    wide[
        "Direction_Aware_2of4"
    ]
].to_csv(
    OUT /
    "real_direction_aware_consensus_2of4.csv",
    index=False
)


# ============================================================
# COMPARISON TABLE
# ============================================================

comparison_columns = [
    column
    for column in [
        "Taxon",
        "Denovo_ID",
        "Phylum",
        "Family",
        "Genus",

        "CLR_Significant",
        "CLR_Adjusted_P_Value",
        "CLR_Effect_Direction",

        "ALDEx2_Significant",
        "ALDEx2_Adjusted_P_Value",
        "ALDEx2_Effect_Direction",

        "ANCOM_Significant",
        "ANCOM_Adjusted_P_Value",
        "ANCOM_Effect_Direction",

        "DESeq2_Significant",
        "DESeq2_Adjusted_P_Value",
        "DESeq2_Effect_Direction",

        "Consensus_Count",
        "Consensus_Class",
        "Consensus_Direction",
        "Direction_Aware_2of4",
    ]
    if column in wide.columns
]


comparison = (
    wide[
        comparison_columns
    ]
    .sort_values(
        [
            "Consensus_Count",
            "Taxon",
        ],
        ascending=[
            False,
            True,
        ]
    )
)


comparison.to_csv(
    OUT /
    "real_four_method_comparison.csv",
    index=False
)


# ============================================================
# AGREEMENT MATRICES
# ============================================================

method_order = list(
    FILES.keys()
)


indexed = {
    name:
        df.set_index(
            "Taxon"
        )
    for name, df
    in tables.items()
}


significant_sets = {
    name:
        set(
            df.loc[
                df[
                    "Significant"
                ],
                "Taxon"
            ]
        )
    for name, df
    in tables.items()
}


jaccard = pd.DataFrame(
    index=method_order,
    columns=method_order,
    dtype=float
)


call_agreement = pd.DataFrame(
    index=method_order,
    columns=method_order,
    dtype=float
)


direction_matrix = pd.DataFrame(
    index=method_order,
    columns=method_order,
    dtype=float
)


taxon_order = sorted(
    reference
)


for method_a in method_order:

    for method_b in method_order:

        set_a = significant_sets[
            method_a
        ]

        set_b = significant_sets[
            method_b
        ]


        if method_a == method_b:

            jaccard.loc[
                method_a,
                method_b
            ] = 1.0


            direction_matrix.loc[
                method_a,
                method_b
            ] = 1.0


        else:

            union = (
                set_a |
                set_b
            )


            if union:

                jaccard.loc[
                    method_a,
                    method_b
                ] = (
                    len(
                        set_a &
                        set_b
                    )
                    /
                    len(
                        union
                    )
                )

            else:

                jaccard.loc[
                    method_a,
                    method_b
                ] = np.nan


            shared = sorted(
                set_a &
                set_b
            )


            if shared:

                direction_a = (
                    indexed[
                        method_a
                    ]
                    .loc[
                        shared,
                        "Effect_Direction"
                    ]
                    .astype(str)
                    .to_numpy()
                )


                direction_b = (
                    indexed[
                        method_b
                    ]
                    .loc[
                        shared,
                        "Effect_Direction"
                    ]
                    .astype(str)
                    .to_numpy()
                )


                direction_matrix.loc[
                    method_a,
                    method_b
                ] = float(
                    (
                        direction_a
                        ==
                        direction_b
                    ).mean()
                )


            else:

                direction_matrix.loc[
                    method_a,
                    method_b
                ] = np.nan


        calls_a = (
            indexed[
                method_a
            ]
            .loc[
                taxon_order,
                "Significant"
            ]
            .astype(bool)
            .to_numpy()
        )


        calls_b = (
            indexed[
                method_b
            ]
            .loc[
                taxon_order,
                "Significant"
            ]
            .astype(bool)
            .to_numpy()
        )


        call_agreement.loc[
            method_a,
            method_b
        ] = float(
            (
                calls_a
                ==
                calls_b
            ).mean()
        )


jaccard.to_csv(
    OUT /
    "real_significant_jaccard_matrix.csv"
)


call_agreement.to_csv(
    OUT /
    "real_significance_call_agreement_matrix.csv"
)


direction_matrix.to_csv(
    OUT /
    "real_direction_agreement_matrix.csv"
)


# ============================================================
# KEY TAXA
# ============================================================

key_taxa = wide.loc[
    wide[
        "Consensus_Count"
    ] >= 1
].copy()


key_taxa = key_taxa.sort_values(
    [
        "Consensus_Count",
        "Taxon",
    ],
    ascending=[
        False,
        True,
    ]
)


key_taxa.to_csv(
    OUT /
    "real_key_taxa.csv",
    index=False
)


# ============================================================
# OTU_00659 COMPARISON
# ============================================================

otu_659 = long_df.loc[
    long_df[
        "Taxon"
    ].eq(
        "OTU_00659"
    )
].copy()


otu_659.to_csv(
    OUT /
    "real_OTU_00659_comparison.csv",
    index=False
)


# ============================================================
# TERMINAL SUMMARY
# ============================================================

consensus_counts = (
    wide[
        "Consensus_Count"
    ]
    .value_counts()
    .reindex(
        range(5),
        fill_value=0
    )
    .sort_index()
)


print(
    "\nMethod summary:"
)

print(
    method_summary
    .to_string(
        index=False
    )
)


print(
    "\nConsensus distribution:"
)


for number in range(5):

    print(
        f"{number}/4 methods significant: "
        f"{int(consensus_counts.loc[number])}"
    )


print(
    ">=2/4 consensus taxa: "
    f"{int(wide['Consensus_2of4'].sum())}"
)


print(
    ">=3/4 consensus taxa: "
    f"{int(wide['Consensus_3of4'].sum())}"
)


print(
    "4/4 consensus taxa: "
    f"{int(wide['Consensus_4of4'].sum())}"
)


print(
    "Direction-aware >=2/4 consensus taxa: "
    f"{int(wide['Direction_Aware_2of4'].sum())}"
)


print(
    "\nKey taxa with >=1 primary method call:"
)


show_columns = [
    column
    for column in [
        "Taxon",
        "Denovo_ID",
        "Phylum",
        "Family",
        "Genus",
        "CLR_Significant",
        "ALDEx2_Significant",
        "ANCOM_Significant",
        "DESeq2_Significant",
        "Consensus_Count",
        "Consensus_Direction",
    ]
    if column in key_taxa.columns
]


print(
    key_taxa[
        show_columns
    ]
    .to_string(
        index=False
    )
)


print(
    "\nOTU_00659 across all four methods:"
)


otu_columns = [
    column
    for column in [
        "Method",
        "P_Value",
        "Adjusted_P_Value",
        "Significant",
        "Effect_Estimate",
        "Effect_Direction",
    ]
    if column in otu_659.columns
]


print(
    otu_659[
        otu_columns
    ]
    .to_string(
        index=False
    )
)


# ============================================================
# TEXT SUMMARY
# ============================================================

summary_text = [

    "PAR_SCHEPERJANS FOUR-METHOD "
    "REAL-DATA CONSOLIDATION",

    "=" * 72,

    "",

    "PRIMARY SIGNIFICANCE DEFINITIONS",

    "CLR/Wilcoxon: BH q < 0.05",

    "ALDEx2: wi.eBH < 0.05",

    (
        "ANCOM-BC2: q < 0.05 AND "
        "passed pseudo-count sensitivity"
    ),

    (
        "DESeq2: primary BH-adjusted "
        "result < 0.05"
    ),

    "",
]


for _, row in method_summary.iterrows():

    summary_text.append(

        f"{row['Method']}: "
        f"tested="
        f"{int(row['Taxa_Tested'])}, "
        f"significant="
        f"{int(row['Primary_Significant'])}, "
        f"min p="
        f"{row['Minimum_P_Value']:.10g}, "
        f"min q="
        f"{row['Minimum_Adjusted_P_Value']:.10g}"
    )


summary_text.extend(
    [
        "",

        "CONSENSUS",

        (
            f"0/4="
            f"{int(consensus_counts.loc[0])}"
        ),

        (
            f"1/4="
            f"{int(consensus_counts.loc[1])}"
        ),

        (
            f"2/4="
            f"{int(consensus_counts.loc[2])}"
        ),

        (
            f"3/4="
            f"{int(consensus_counts.loc[3])}"
        ),

        (
            f"4/4="
            f"{int(consensus_counts.loc[4])}"
        ),

        (
            f">=2/4="
            f"{int(wide['Consensus_2of4'].sum())}"
        ),

        (
            f">=3/4="
            f"{int(wide['Consensus_3of4'].sum())}"
        ),

        (
            f"4/4="
            f"{int(wide['Consensus_4of4'].sum())}"
        ),

        (
            "Direction-aware >=2/4="
            f"{int(wide['Direction_Aware_2of4'].sum())}"
        ),

        "",

        "INTERPRETATION",

        (
            "These are real-data agreement results. "
            "Ground truth is unknown, so consensus "
            "does not establish biological correctness."
        ),
    ]
)


(
    OUT /
    "real_consolidation_summary.txt"
).write_text(
    "\n".join(
        summary_text
    ),
    encoding="utf-8"
)


print(
    "\n" + "=" * 72
)

print(
    "FOUR-METHOD REAL-DATA "
    "CONSOLIDATION: PASS"
)

print(
    "=" * 72
)
