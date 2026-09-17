"""
Reanalysis of the Fig.4 consensus co-expression network for npjSBA Rev3.
Addresses Reviewer 3 Major Concern 2: treatment-pooled correlations without
multiple-testing correction, and the 28/27/19 consensus-pair-count discrepancy.

Inputs (already exported by the original pipeline, unchanged):
  Figures/Fig_4c_Hes1_normalized_data.csv
  Figures/Fig_4c_Kcnj11_normalized_data.csv
  Figures/Fig_4c_Ins1_normalized_data.csv
  Figures/Fig_4c_Raw_expression_data.csv

Row order (12 rows = 4 reps x 3 groups) inferred from Hes1 expression pattern
(Notch target; DAPT is a Notch inhibitor so Hes1 should be lowest in DAPT):
  rows 0-3  = CTRL
  rows 4-7  = DAPT
  rows 8-11 = DKK-1
"""
import os
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests
from itertools import combinations

# Run from a directory containing the R_scripts_Rev1 pipeline's Fig_4c_*.csv
# outputs (see README.md, Input Data Files Required).
BASE = "."

files = {
    "Hes1-normalized": "Fig_4c_Hes1_normalized_data.csv",
    "Kcnj11-normalized": "Fig_4c_Kcnj11_normalized_data.csv",
    "Ins1-normalized": "Fig_4c_Ins1_normalized_data.csv",
    "Raw": "Fig_4c_Raw_expression_data.csv",
}

group = np.array(["CTRL"] * 4 + ["DAPT"] * 4 + ["DKK1"] * 4)

def sanity_check_group_order(raw_df):
    # Hes1 is not in the normalized files (used as denominator), but IS in Raw.
    if "Hes1" in raw_df.columns:
        means = raw_df.groupby(group)["Hes1"].mean()
        print("Hes1 mean by assumed group (expect DAPT lowest):")
        print(means, "\n")

def pooled_correlations(df):
    genes = df.columns.tolist()
    rows = []
    for g1, g2 in combinations(genes, 2):
        r, p = pearsonr(df[g1], df[g2])
        rows.append((g1, g2, r, p))
    out = pd.DataFrame(rows, columns=["gene1", "gene2", "r", "p_raw"])
    out["p_bh"] = multipletests(out["p_raw"], method="fdr_bh")[1]
    return out

def treatment_adjusted_correlations(df):
    # residualize each gene on group (one-hot, drop first) via OLS, then correlate residuals
    dummies = pd.get_dummies(group, drop_first=True).astype(float)
    X = np.column_stack([np.ones(len(df)), dummies.values])
    resid = {}
    for g in df.columns:
        y = df[g].values.astype(float)
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid[g] = y - X @ beta
    resid_df = pd.DataFrame(resid)
    genes = df.columns.tolist()
    rows = []
    for g1, g2 in combinations(genes, 2):
        r, p = pearsonr(resid_df[g1], resid_df[g2])
        rows.append((g1, g2, r, p))
    out = pd.DataFrame(rows, columns=["gene1", "gene2", "r", "p_raw"])
    out["p_bh"] = multipletests(out["p_raw"], method="fdr_bh")[1]
    return out

# CORRECTION (14-Sep-2026): Fig_4c_*.csv include a Ryr1 column, but the manuscript's
# 20-gene correlation panel deliberately excludes Ryr1 (Methods: "excluding Ryr1,
# undetected across all conditions"; also absent from Fig 3e's gene list and from
# Fig_3g_..._exclude_Ryr1.svg). Dropping Ryr1 here reproduces the manuscript's existing
# per-scheme significant-pair counts (169/118/98/80) and 27-pair consensus exactly --
# confirming Ryr1's presence in the raw CSVs is a pipeline artefact, not intended input.
DROP_GENES = ["Ryr1"]

results = {}
for label, fname in files.items():
    df = pd.read_csv(os.path.join(BASE, fname))
    df = df.drop(columns=[g for g in DROP_GENES if g in df.columns])
    if label == "Raw":
        sanity_check_group_order(df)
    results[label] = {
        "pooled": pooled_correlations(df),
        "adjusted": treatment_adjusted_correlations(df),
    }

def sig_pairs(df, pcol, thresh=0.05):
    sub = df[df[pcol] < thresh]
    return set(tuple(sorted((r.gene1, r.gene2))) for r in sub.itertuples())

print("=" * 70)
print("PER-METHOD SIGNIFICANT PAIR COUNTS")
print("=" * 70)
for scheme in ["pooled", "adjusted"]:
    for pcol, corr_label in [("p_raw", "uncorrected p<0.05"), ("p_bh", "BH-corrected q<0.05")]:
        counts = {label: len(sig_pairs(results[label][scheme], pcol)) for label in files}
        print(f"{scheme:10s} | {corr_label:22s} | " + ", ".join(f"{k}={v}" for k, v in counts.items()))

print()
print("=" * 70)
print("CONSENSUS PAIR COUNTS (intersection across all 4 normalization schemes)")
print("=" * 70)

consensus = {}
for scheme in ["pooled", "adjusted"]:
    for pcol, corr_label in [("p_raw", "uncorrected"), ("p_bh", "BH-corrected")]:
        sets = [sig_pairs(results[label][scheme], pcol) for label in files]
        cons = set.intersection(*sets)
        consensus[(scheme, corr_label)] = cons
        print(f"{scheme:10s} x {corr_label:13s} -> {len(cons)} consensus pairs")

print()
print("=" * 70)
print("REPRODUCE CURRENT MANUSCRIPT METHOD (pooled, uncorrected) — sanity check vs code's '28'")
print("=" * 70)
print(sorted(consensus[("pooled", "uncorrected")]))

print()
print("=" * 70)
print("CORRECTLY-SPECIFIED CONSENSUS (treatment-adjusted, BH-corrected)")
print("=" * 70)
final = consensus[("adjusted", "BH-corrected")]
print(f"N = {len(final)}")
for pair in sorted(final):
    print(" ", pair)

print()
print("=" * 70)
print("OVERLAP: how many of the current 28 survive proper treatment-adjustment + BH?")
print("=" * 70)
old = consensus[("pooled", "uncorrected")]
new = consensus[("adjusted", "BH-corrected")]
print(f"Old (pooled, uncorrected): {len(old)}")
print(f"New (adjusted, BH): {len(new)}")
print(f"Overlap: {len(old & new)}")
print(f"Lost (in old, not new): {len(old - new)}")
print(f"Gained (in new, not old): {len(new - old)}")
