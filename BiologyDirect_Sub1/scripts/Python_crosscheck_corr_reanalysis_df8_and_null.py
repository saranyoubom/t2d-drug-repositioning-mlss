"""
Treatment-adjusted correlation reanalysis with CORRECT degrees of freedom, plus a
permutation null that shows how many 'consensus' pairs arise from the 3-group
structure alone (no within-group co-regulation).

Inputs (relative path): data/Fig_4c_{Raw_expression,Hes1_normalized,Kcnj11_normalized,Ins1_normalized}_data.csv
Row order (12 rows = 4 technical replicates x 3 groups): rows 0-3 CTRL, 4-7 DAPT, 8-11 DKK-1
(inferred from Hes1 and validated: group-mean Hes1 ratios reproduce the manuscript's 0.14 / 0.52).
Ryr1 is dropped (excluded everywhere else in the study).

Run:  python corr_reanalysis_df8_and_null.py <data_dir> [n_perm]
"""
import sys, numpy as np, pandas as pd
from itertools import combinations
from scipy import stats
from statsmodels.stats.multitest import multipletests

data_dir = sys.argv[1] if len(sys.argv) > 1 else "data"
n_perm = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
rng = np.random.default_rng(20261014)

files = {"Raw": "Fig_4c_Raw_expression_data.csv", "Hes1": "Fig_4c_Hes1_normalized_data.csv",
         "Kcnj11": "Fig_4c_Kcnj11_normalized_data.csv", "Ins1": "Fig_4c_Ins1_normalized_data.csv"}
S = {k: pd.read_csv(f"{data_dir}/{v}").drop(columns=["Ryr1"]) for k, v in files.items()}
grp = np.repeat([0, 1, 2], 4)
X = np.column_stack([np.ones(12), (grp == 1).astype(float), (grp == 2).astype(float)])

def resid(M):                      # M: 12 x g  -> residuals after removing group means
    H = X @ np.linalg.pinv(X)
    return M - H @ M

def corr_p(M, df):                 # pairwise Pearson r and two-sided p with given df
    r = np.corrcoef(M.T)
    iu = np.triu_indices(r.shape[0], 1)
    rv = r[iu]
    t = rv * np.sqrt(df / np.maximum(1 - rv ** 2, 1e-12))
    p = 2 * stats.t.sf(np.abs(t), df)
    return iu, rv, p

def schemes_from_raw(raw_df):       # rebuild the four schemes from a raw (2^-dCt) matrix
    out = {"Raw": raw_df}
    for den in ["Hes1", "Kcnj11", "Ins1"]:
        out[den] = raw_df.drop(columns=[den]).div(raw_df[den], axis=0)
    return out

def analyse(sch, adjusted):
    res = {}
    for k, df in sch.items():
        M = df.values.astype(float)
        if adjusted:
            M = resid(M); dfree = 12 - 2 - 2        # partial-correlation df = 8
        else:
            dfree = 12 - 2                          # pooled df = 10
        iu, r, p = corr_p(M, dfree)
        q = multipletests(p, method="fdr_bh")[1]
        genes = df.columns.values
        res[k] = pd.DataFrame({"g1": genes[iu[0]], "g2": genes[iu[1]], "r": r, "p": p, "q": q})
    return res

def sig_set(df, col):
    d = df[df[col] < 0.05]; return set(zip(d.g1, d.g2))

def summarize(res, label):
    print(f"\n== {label}")
    for col, nm in [("p", "uncorrected p<0.05"), ("q", "BH q<0.05")]:
        counts = {k: len(sig_set(v, col)) for k, v in res.items()}
        sets = [sig_set(v, col) for v in res.values()]
        # canonicalise pair ordering across schemes
        sets = [{tuple(sorted(x)) for x in s} for s in sets]
        cons = set.intersection(*sets)
        print(f"  {nm:20s} per-scheme {counts} | consensus(all 4) = {len(cons)}")
    return res

print("N rows =", 12, "| pooled df = 10 | treatment-adjusted (partial) df = 8")
pooled = summarize(analyse(S, False), "POOLED across 3 treatment groups")
adj8 = summarize(analyse(S, True), "TREATMENT-ADJUSTED (residuals on group), correct df = 8")

# sign consistency of the pooled p<0.05 consensus (the manuscript's stated criterion)
sets = [{tuple(sorted(x)) for x in sig_set(v, "p")} for v in pooled.values()]
cons = sorted(set.intersection(*sets))
sign = {k: {tuple(sorted((a, b))): r for a, b, r in zip(v.g1, v.g2, v.r)} for k, v in pooled.items()}
incons = [pr for pr in cons if len({np.sign(sign[k][pr]) for k in sign}) > 1]
print(f"\nPooled p<0.05 consensus = {len(cons)}; sign-INCONSISTENT across schemes = {len(incons)} -> sign-consistent consensus = {len(cons)-len(incons)}")
print("  inconsistent:", incons)

# ---- permutation null: keep each gene's group means and within-group variance; break within-group co-regulation
raw = S["Raw"].copy()
L = np.log(raw.values.astype(float))
Hm = X @ np.linalg.pinv(X)
fitted = Hm @ L
E = L - fitted
cons_counts, cons_signed, per_scheme = [], [], {k: [] for k in files}
for _ in range(n_perm):
    Ep = E.copy()
    for g in range(E.shape[1]):
        for gi in range(3):
            idx = np.where(grp == gi)[0]
            Ep[idx, g] = rng.permutation(E[idx, g])
    rawp = pd.DataFrame(np.exp(fitted + Ep), columns=raw.columns)
    res = analyse(schemes_from_raw(rawp), False)
    sets = [{tuple(sorted(x)) for x in sig_set(v, "p")} for v in res.values()]
    c = set.intersection(*sets)
    cons_counts.append(len(c))
    sg = {k: {tuple(sorted((a, b))): r for a, b, r in zip(v.g1, v.g2, v.r)} for k, v in res.items()}
    cons_signed.append(sum(1 for pr in c if len({np.sign(sg[k][pr]) for k in sg}) == 1))
    for k, v in res.items(): per_scheme[k].append(len(sig_set(v, "p")))
cc = np.array(cons_counts); cs = np.array(cons_signed)
print(f"\n== PERMUTATION NULL ({n_perm} reps; independent within-group residuals; group means retained)")
for k in files:
    a = np.array(per_scheme[k]); print(f"  {k:7s} pooled p<0.05 pairs: mean {a.mean():.1f} (95% {np.percentile(a,2.5):.0f}-{np.percentile(a,97.5):.0f})")
print(f"  4-scheme consensus pairs: mean {cc.mean():.1f} (95% {np.percentile(cc,2.5):.0f}-{np.percentile(cc,97.5):.0f}); P(null >= 27) = {(cc>=27).mean():.3f}")
print(f"  ... sign-consistent consensus: mean {cs.mean():.1f} (95% {np.percentile(cs,2.5):.0f}-{np.percentile(cs,97.5):.0f}); P(null >= 21) = {(cs>=21).mean():.3f}")
