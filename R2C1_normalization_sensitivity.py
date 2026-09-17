"""
R2C1: normalization-gene sensitivity, extended to genes from FUNCTIONAL CATEGORIES
not used as denominators in the original scheme (Hes1=Notch, Kcnj11+Ins1=beta-cell
function). Reviewer 2 explicitly asks about categories not yet tested: Wnt pathway,
Calcium signaling.

Candidates chosen: Tcf7 (Wnt pathway, high-expression, low missingness) and
Cacna1c (Calcium signaling), mirroring the Hey1 (Notch) sensitivity check already
done in Rev2.

Reports pair recovery under BOTH:
  (1) the original pooled+uncorrected method (for continuity with Rev1/Rev2 style)
  (2) the corrected treatment-adjusted+BH method (for consistency with the R3C2 fix)
"""
import os
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests
from itertools import combinations

# Run from a directory containing the R_scripts_Rev1 pipeline's
# Fig_4c_Raw_expression_data.csv output (see README.md, Input Data Files Required).
BASE = "."
raw = pd.read_csv(os.path.join(BASE, "Fig_4c_Raw_expression_data.csv"))
# CORRECTION (14-Sep-2026): Ryr1 is deliberately excluded from the manuscript's 20-gene
# correlation panel (Methods; also absent from Fig 3e/3g elsewhere) but is present in this
# raw CSV as a pipeline artefact. Dropping it reproduces the manuscript's existing 27-pair
# consensus exactly -- see Fig4_treatment_adjusted_reanalysis.py for the full verification.
raw = raw.drop(columns=["Ryr1"])

group = np.array(["CTRL"] * 4 + ["DAPT"] * 4 + ["DKK1"] * 4)

def normalize_by(df, denom_gene):
    out = df.copy()
    for col in out.columns:
        if col != denom_gene:
            out[col] = out[col] / df[denom_gene]
    return out.drop(columns=[denom_gene])

def pooled_sig_pairs(df, thresh=0.05):
    genes = df.columns.tolist()
    pairs = set()
    for g1, g2 in combinations(genes, 2):
        r, p = pearsonr(df[g1], df[g2])
        if p < thresh:
            pairs.add(tuple(sorted((g1, g2))))
    return pairs

def adjusted_bh_sig_pairs(df, thresh=0.05):
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
        rows.append((g1, g2, p))
    out = pd.DataFrame(rows, columns=["g1", "g2", "p"])
    out["p_bh"] = multipletests(out["p"], method="fdr_bh")[1]
    return set(tuple(sorted((r.g1, r.g2))) for r in out.itertuples() if r.p_bh < thresh)

# Original 4 schemes (raw, Hes1, Kcnj11, Ins1) -- pooled/uncorrected
orig_schemes = {
    "raw": raw,
    "Hes1": normalize_by(raw, "Hes1"),
    "Kcnj11": normalize_by(raw, "Kcnj11"),
    "Ins1": normalize_by(raw, "Ins1"),
}
orig_pooled = {k: pooled_sig_pairs(v) for k, v in orig_schemes.items()}
orig_consensus_pooled = set.intersection(*orig_pooled.values())
print("Original 4-scheme consensus (pooled, uncorrected):", len(orig_consensus_pooled))

# New candidate normalizers from unused functional categories
new_candidates = {
    "Tcf7 (Wnt pathway)": "Tcf7",
    "Cacna1c (Calcium signaling)": "Cacna1c",
}

print("\n" + "=" * 70)
print("PART 1: pooled/uncorrected method (continuity with Rev1/Rev2 style)")
print("=" * 70)
for label, gene in new_candidates.items():
    new_norm = normalize_by(raw, gene)
    new_pairs = pooled_sig_pairs(new_norm)
    # A consensus pair involving the substituted gene itself is structurally
    # non-evaluable once that gene is the normalizer denominator (it cannot be
    # correlated against its own normalized ratio) -- exclude it from the
    # evaluable-pair denominator rather than counting it as "lost".
    evaluable = {p for p in orig_consensus_pooled if gene not in p}
    n_excluded = len(orig_consensus_pooled) - len(evaluable)
    retained = evaluable & new_pairs
    lost = evaluable - new_pairs
    print(f"\n{label}: {len(new_pairs)} significant pairs (p<0.05, pooled)")
    print(f"  {n_excluded} of {len(orig_consensus_pooled)} original consensus pairs are structurally "
          f"non-evaluable under this normalizer (involve {gene} itself)")
    print(f"  Retention of {len(evaluable)} evaluable pairs: "
          f"{len(retained)}/{len(evaluable)} ({100*len(retained)/len(evaluable):.1f}%)")
    print(f"  Lost from evaluable consensus under this normalizer: {sorted(lost)}")

# 5-way consensus (original 4 + each new one), pooled/uncorrected
for label, gene in new_candidates.items():
    new_norm = normalize_by(raw, gene)
    new_pairs = pooled_sig_pairs(new_norm)
    five_way = orig_consensus_pooled & new_pairs
    print(f"\n5-way consensus (orig 4 + {label}): {len(five_way)} pairs")

six_way = orig_consensus_pooled.copy()
for label, gene in new_candidates.items():
    six_way &= pooled_sig_pairs(normalize_by(raw, gene))
print(f"\n6-way consensus (orig 4 + BOTH new normalizers): {len(six_way)} pairs")
print(sorted(six_way))

print("\n" + "=" * 70)
print("PART 2: treatment-adjusted + BH-corrected method (consistency w/ R3C2 fix)")
print("=" * 70)
adj_schemes = dict(orig_schemes)
for label, gene in new_candidates.items():
    adj_schemes[label] = normalize_by(raw, gene)

adj_bh = {k: adjusted_bh_sig_pairs(v) for k, v in adj_schemes.items()}
for k, v in adj_bh.items():
    print(f"{k}: {len(v)} pairs survive treatment-adjustment + BH")

all_six_adj = set.intersection(*adj_bh.values())
print(f"\n6-way consensus, treatment-adjusted + BH: {len(all_six_adj)} pairs")
