"""Independent Python port of R_17genes/Fig_6c_MLSS_v4.0.r (MLSS v4.0), plus the upstream network step of
Fig_4c/4d/4e/6b (four normalization schemes -> consensus pairs -> degree) so that the INS-1 matrix can be carried through.

MLSS = [0.45 C + 0.30 B + 0.10 V + 0.15 P] x 9 + Synergy - Antagonism, with
  C = 1 - |A n B| / |A u B| on the sets of hub-gene targets, B = 1 - |nA - nB| / (nA + nB),
  V = min-max normalization (over all combinations) of the mean |r| of the cross-drug target-gene pairs that have a
      correlation in the lookup, P = min(max(mean pIC50 - 4, 0), 5) / 5 (mean of the two drug means; 0 if a drug has none),
  Synergy = sum r_ij x m_ij over cross-drug target pairs with r_ij > 0, Antagonism = sum |r_ij| x m_ij over r_ij < 0,
  m_ij = 0 if either pIC50 is missing, else min(max(pIC50_i, pIC50_j) / 9, 1).
Details of the R code that the header comment does not state (all reproduced here):
  * a drug enters only if a DrugCentral row has GENE equal (after upper-casing) to a hub gene and ORGANISM == 'Homo sapiens';
    rows with GENE 'CACNA1D|CACNA1C' etc. are NOT matched (exact string match), whatever the action type;
  * pIC50 = MAX of ACT_VALUE over rows with ACT_TYPE in IC50/Ki/EC50, RELATION '=', ORGANISM 'Homo sapiens' (types are mixed);
  * the correlation lookup is the Hes1-normalized pooled Pearson matrix, NA where pooled p >= 0.05, independent of consensus status;
  * ties in the final score keep the original row order (expand.grid: Drug_A fastest, Drug_B slowest, Drug_A < Drug_B).
Usage (gate against the existing R output; writes gate_output.txt and gate_scores.csv next to this file):
  python -I mlss_port.py
Imported by mlss_scenarios.py and direction_pic50.py (they add this folder to sys.path)."""
import os, sys, itertools, json, math
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.dirname(os.path.dirname(HERE))                      # .../Analysis_outputs
RDIR = os.environ.get('MLSS_RDIR', os.path.join(ANALYSIS, 'R_17genes'))   # read-only inputs of the author R scripts
SBA_MATRIX = os.environ.get('MLSS_INS1_MATRIX', os.path.join(os.path.dirname(os.path.dirname(ANALYSIS)), 'SBA', 'Figures',
                                                           'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx'))
GAPDH_JSON = os.path.join(ANALYSIS, 'reanalysis_17genes', 'gapdh_exact.json')
W_STD = (0.45, 0.30, 0.10, 0.15)       # C, B, V, P
BASE_MULT = 9.0
GENES17 = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']
GROUP = np.array([0] * 4 + [1] * 4 + [2] * 4)                           # rows 0-3 CTRL, 4-7 DAPT, 8-11 DKK-1
DENOMINATORS = ['Hes1', 'Kcnj11', 'Ins1']
KEY = frozenset(['CACNA1D', 'GLP1R'])                                   # the only correlated cross-drug pair in the reported run


# ------------------------------------------------------------------ inputs of Fig_6c (steps 1-5 of the R script)
def read_drugcentral(path=None):
    """drug.target.interaction.csv as text (no silent NA conversion); ACT_VALUE numeric (NaN where empty)."""
    path = path or os.path.join(RDIR, 'drug.target.interaction.csv')
    dc = pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[''], low_memory=False)
    dc['ACT_VALUE'] = pd.to_numeric(dc['ACT_VALUE'], errors='coerce')
    dc['G'] = dc['GENE'].str.upper()
    return dc


def potency_from(dc):
    """dict (drug, GENE) -> max ACT_VALUE over human IC50/Ki/EC50 rows with RELATION '=' (R: potency_lookup)."""
    p = dc[dc['ACT_TYPE'].isin(['IC50', 'Ki', 'EC50']) & (dc['ORGANISM'] == 'Homo sapiens') & (dc['RELATION'] == '=') & dc['ACT_VALUE'].notna()]
    return {k: float(v) for k, v in p.groupby(['DRUG_NAME', 'G'])['ACT_VALUE'].max().items()}


def targets_from(dc, hubs):
    """dict drug -> list of hub genes (unique, order of first appearance); R: drug_targets %>% get_drug_targets."""
    hubs = {h.upper() for h in hubs}
    t = dc[dc['G'].isin(hubs) & (dc['ORGANISM'] == 'Homo sapiens')][['DRUG_NAME', 'GENE', 'TARGET_NAME', 'ACTION_TYPE', 'G']].drop_duplicates(subset=['DRUG_NAME', 'GENE', 'TARGET_NAME', 'ACTION_TYPE'])
    out = {}
    for d, g in zip(t['DRUG_NAME'], t['G']):
        out.setdefault(d, [])
        if g not in out[d]:
            out[d].append(g)
    return out


def hubs_from_centrality(path=None, min_degree=5):
    path = path or os.path.join(RDIR, 'Fig_6b_Gene_Centrality_Summary.csv')
    h = pd.read_csv(path)
    return [g.upper() for g in h.loc[h['Degree'] >= min_degree, 'Gene']], h


def corr_from_matrix_csv(path=None):
    """Fig_4e_Hes1_significant_correlations_only.csv -> dict frozenset({G1,G2}) -> r (NA cells dropped)."""
    path = path or os.path.join(RDIR, 'Fig_4e_Hes1_significant_correlations_only.csv')
    m = pd.read_csv(path, index_col=0)
    genes = [g.upper() for g in m.index]
    out = {}
    for i, a in enumerate(genes):
        for j, b in enumerate(m.columns):
            b = b.upper()
            if a != b and not pd.isna(m.iloc[i, j]):
                out.setdefault(frozenset((a, b)), float(m.iloc[i, j]))
    return out


# ------------------------------------------------------------------ the scoring (steps 7-14)
def r_sort_key(s):
    """Approximation of R's locale collation (case-insensitive, then code point); identical to code-point order for the
    lower-case names of the reported run ('(' sorts before letters in both)."""
    return (s.casefold(), s)


def combos_in_r_order(drugs):
    """R: expand.grid(Drug_A, Drug_B) %>% filter(Drug_A < Drug_B): Drug_A varies fastest."""
    ds = sorted(drugs, key=r_sort_key)
    out = []
    for j, b in enumerate(ds):
        for i, a in enumerate(ds):
            if r_sort_key(a) < r_sort_key(b):
                out.append((a, b))
    return out


def mean_pic50(drug, tg, pot):
    v = [pot[(drug, g)] for g in dict.fromkeys(tg) if (drug, g) in pot]
    return float(np.mean(v)) if v else 0.0


def mult(da, gi, db, gj, pot):
    if (da, gi) not in pot or (db, gj) not in pot:
        return 0.0
    return min(max(pot[(da, gi)], pot[(db, gj)]) / 9.0, 1.0)


def mlss_table(targets, pot, corr, weights=W_STD, base_mult=BASE_MULT, corr_v=None):
    """targets: drug -> [GENE,...]; pot: (drug, GENE) -> pIC50; corr: frozenset({G1,G2}) -> r (NA pairs absent).
    corr is used for Synergy, Antagonism and Mean_Pos/Neg_r, and for V as well unless corr_v is given (the R script uses one
    lookup for all of them; corr_v only serves the scenario 'change r in the Synergy Bonus but keep V').
    Returns a DataFrame in the final R order (descending MLSS, ties by original row order) with the R column names."""
    wC, wB, wV, wP = weights
    cv = corr if corr_v is None else corr_v
    rows = []
    for a, b in combos_in_r_order(targets):
        ta, tb = list(dict.fromkeys(targets[a])), list(dict.fromkeys(targets[b]))
        inter, union = len(set(ta) & set(tb)), len(set(ta) | set(tb))
        C = 1 - inter / union if union > 0 else 0.0
        B = 1 - abs(len(ta) - len(tb)) / (len(ta) + len(tb)) if (len(ta) + len(tb)) > 0 else 0.0
        rs = [corr[frozenset((gi, gj))] for gi in ta for gj in tb if frozenset((gi, gj)) in corr and gi != gj]
        rsv = [cv[frozenset((gi, gj))] for gi in ta for gj in tb if frozenset((gi, gj)) in cv and gi != gj]
        pos = [r for r in rs if r > 0]; neg = [r for r in rs if r < 0]
        mabs = float(np.mean([abs(r) for r in rsv])) if rsv else 0.0
        mA, mB = mean_pic50(a, ta, pot), mean_pic50(b, tb, pot)
        mp = (mA + mB) / 2.0
        P = min(max(mp - 4, 0), 5) / 5
        syn = ant = 0.0
        for gi in ta:
            for gj in tb:
                k = frozenset((gi, gj))
                if gi == gj or k not in corr:
                    continue
                r = corr[k]
                if r > 0:
                    syn += r * mult(a, gi, b, gj, pot) * 1
                elif r < 0:
                    ant += abs(r) * mult(a, gi, b, gj, pot) * 1
        rows.append(dict(Drug_A=a, Drug_B=b, Targets_A='|'.join(ta), Targets_B='|'.join(tb), n_targets_A=len(ta), n_targets_B=len(tb),
                         n_targets_combined=len(ta) + len(tb), Mean_pIC50_A=mA, Mean_pIC50_B=mB, Mean_pIC50=mp, Complementarity_C=C,
                         Balance_B=B, Mean_Abs_r=mabs, Mean_Pos_r=float(np.mean(pos)) if pos else 0.0,
                         Mean_Neg_r=float(np.mean(neg)) if neg else 0.0, Potency_P=P, Synergy_Bonus=syn, Antagonism_Penalty=ant))
    t = pd.DataFrame(rows)
    if len(t) == 0:
        return t
    lo, hi = t['Mean_Abs_r'].min(), t['Mean_Abs_r'].max()
    t['Network_Coverage_V'] = (t['Mean_Abs_r'] - lo) / (hi - lo) if hi > lo else 0.0
    t['Base_Score'] = (wC * t['Complementarity_C'] + wB * t['Balance_B'] + wV * t['Network_Coverage_V'] + wP * t['Potency_P']) * base_mult
    t['MLSS_v4.0'] = t['Base_Score'] + t['Synergy_Bonus'] - t['Antagonism_Penalty']
    t['Combination'] = t['Drug_A'] + ' + ' + t['Drug_B']
    t = t.iloc[np.argsort(-t['MLSS_v4.0'].values, kind='stable')].reset_index(drop=True)      # arrange(desc()) is stable
    t['Rank'] = np.arange(1, len(t) + 1)
    return t


def run_reported(dc=None, hubs=None, corr=None):
    """The reported run, from the files the R script reads."""
    dc = read_drugcentral() if dc is None else dc
    hubs = hubs if hubs is not None else hubs_from_centrality()[0]
    corr = corr if corr is not None else corr_from_matrix_csv()
    return mlss_table(targets_from(dc, hubs), potency_from(dc), corr)


# ------------------------------------------------------------------ upstream: INS-1 matrix -> four schemes -> consensus -> degree
def load_matrix(path=None):
    """INS-1 relative expression (12 x 17); the plan path (SBA/Figures, 21 genes) is used if present, else the R_17genes copy."""
    for p in ([path] if path else [SBA_MATRIX, os.path.join(RDIR, 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx')]):
        if p and os.path.exists(p):
            return pd.read_excel(p)[GENES17].astype(float).reset_index(drop=True), p
    raise FileNotFoundError('INS-1 matrix not found')


def gapdh_variant(X, which):
    """C0 / C1 recalculation of the 26 August DAPT values of five genes (constants from reanalysis_17genes/gapdh_exact.py)."""
    if which == 'reported':
        return X.copy()
    hk = np.array([22.5685093010354, 22.7982726003435, 22.5685093010354, 22.7982726003435])
    raw = {'A10': 24.9873996418614, 'B10': 24.1318270358876, 'C10': 22.5685093010354, 'D10': 22.7982726003435}
    ref = {'C0': np.mean([raw['C10'], raw['D10']]), 'C1': np.mean(list(raw.values()))}[which]
    d = X.copy()
    for g in ['Tcf7l2', 'Kcnj11', 'Cacna1d', 'Cacna1c', 'Ptbp1']:
        d.loc[4:7, g] = X.loc[4:7, g].values * 2.0 ** (ref - hk)
    return d


def pearson_p(M, df_resid=None):
    M = np.asarray(M, float); n = M.shape[0]
    r = np.corrcoef(M, rowvar=False)
    df = (n - 2) if df_resid is None else df_resid
    with np.errstate(divide='ignore', invalid='ignore'):
        t = r * np.sqrt(df / (1 - r ** 2)); p = 2 * stats.t.sf(np.abs(t), df)
    np.fill_diagonal(p, np.nan)
    return r, p


def four_schemes(X):
    """Fig_4c: raw; each other gene divided by Hes1 / Kcnj11 / Ins1 (the denominator column is dropped)."""
    out = {'Raw': X}
    for d in DENOMINATORS:
        out[d + '-normalized'] = X.drop(columns=[d]).div(X[d], axis=0)
    return out


def network_from_matrix(X):
    """Fig_4c/4d/4e/6b on a 12 x 17 matrix: per-scheme pooled r and p, significant pairs (p < 0.05), consensus (all four), degree."""
    sch = four_schemes(X); res = {}
    for k, m in sch.items():
        r, p = pearson_p(m.values); g = list(m.columns)
        res[k] = {frozenset((g[i], g[j])): (float(r[i, j]), float(p[i, j])) for i, j in itertools.combinations(range(len(g)), 2)}
    sig = {k: {pr for pr, (r, p) in d.items() if p < 0.05} for k, d in res.items()}
    cons = set.intersection(*sig.values())
    deg = {}
    for pr in cons:
        for g in pr:
            deg[g] = deg.get(g, 0) + 1
    return dict(res=res, sig_counts={k: len(v) for k, v in sig.items()}, consensus=sorted(tuple(sorted(c)) for c in cons), degree=deg)


def hes1_lookup(net, alpha=0.05, source='Hes1-normalized', signif_only=True):
    """MLSS correlation lookup from a scheme: frozenset({G1,G2}) -> r, keeping only p < alpha when signif_only (as the R CSV does)."""
    return {frozenset((a.upper(), b.upper())): r for pr, (r, p) in net['res'][source].items() for a, b in [tuple(pr)] if (p < alpha or not signif_only)}


def adjusted_lookup(X, log=False, df_resid=8, alpha=0.05, signif_only=True):
    """Treatment-adjusted r of the raw scheme: residuals after removing group means (linear scale, as pipeline.py / the Methods;
    log=True uses the log scale of protocol 2.2), df = 12 - 2 - 2 = 8."""
    M = np.log(X.values) if log else X.values
    R = M.copy()
    for g in range(3):
        R[GROUP == g] = M[GROUP == g] - M[GROUP == g].mean(0)
    r, p = pearson_p(R, df_resid=df_resid); cols = list(X.columns); out = {}
    for i, j in itertools.combinations(range(len(cols)), 2):
        if p[i, j] < alpha or not signif_only:
            out[frozenset((cols[i].upper(), cols[j].upper()))] = float(r[i, j])
    return out


def centrality_degree_table(deg):
    return pd.DataFrame(sorted(deg.items(), key=lambda kv: (-kv[1], kv[0])), columns=['Gene', 'Degree'])


# ------------------------------------------------------------------ gate
def gate(lines, out_csv=None):
    def log(s=''):
        lines.append(s); print(s)
    R = pd.read_csv(os.path.join(RDIR, 'Fig_6c_MLSS_ALL_Combinations.csv'))
    hubs, hub_df = hubs_from_centrality()
    dc = read_drugcentral(); corr = corr_from_matrix_csv()
    P = run_reported(dc, hubs, corr)
    log('MLSS port gate: Python port vs existing R output Fig_6c_MLSS_ALL_Combinations.csv')
    log('Inputs (read only, R_17genes): Fig_6b_Gene_Centrality_Summary.csv (Degree >= 5 -> hubs %s), Fig_4e_Hes1_significant_correlations_only.csv (%d correlated pairs), drug.target.interaction.csv (%d rows)' % (hubs, len(corr), len(dc)))
    log('Port: %d combinations; R output: %d rows' % (len(P), len(R)))
    R['key'] = [frozenset((a, b)) for a, b in zip(R['Drug_A'], R['Drug_B'])]
    P['key'] = [frozenset((a, b)) for a, b in zip(P['Drug_A'], P['Drug_B'])]
    M = R.merge(P, on='key', suffixes=('_R', '_py'))
    log('Combinations matched by drug pair: %d of %d' % (len(M), len(R)))
    d = (M['MLSS_v4.0_R'] - M['MLSS_v4.0_py']).abs()
    log('MLSS score: max |difference| = %.3e (criterion < 0.005); mean |difference| = %.3e' % (d.max(), d.mean()))
    for c in ['Complementarity_C', 'Balance_B', 'Network_Coverage_V', 'Potency_P', 'Mean_pIC50', 'Mean_Abs_r', 'Mean_Pos_r', 'Mean_Neg_r', 'Base_Score', 'Synergy_Bonus', 'Antagonism_Penalty']:
        log('  %-20s max |difference| = %.3e' % (c, (M[c + '_R'] - M[c + '_py']).abs().max()))
    same_lab = int((M['Combination_R'] == M['Combination_py']).sum())
    log('Combination label (A + B order) identical: %d of %d' % (same_lab, len(M)))
    same_t = int(((M['Targets_A_R'] == M['Targets_A_py']) & (M['Targets_B_R'] == M['Targets_B_py'])).sum())
    log('Target strings identical: %d of %d' % (same_t, len(M)))
    # rank order
    seq_R = list(R.sort_values('Rank')['key']); seq_P = list(P.sort_values('Rank')['key'])
    exact = sum(1 for a, b in zip(seq_R, seq_P) if a == b)
    log('Rank order: identical combination at the same rank in %d of %d positions' % (exact, len(seq_R)))
    # tie-tolerant comparison: group by score rounded to 9 decimals (R scores tie exactly in the CSV at 15 digits)
    R['s9'] = R['MLSS_v4.0'].round(9); P['s9'] = P['MLSS_v4.0'].round(9)
    ok = True
    for (s, g) in R.groupby('s9'):
        setR = set(g['key']); setP = set(P.loc[P['s9'] == s, 'key'])
        if setR != setP:
            ok = False
    ranksR = R.set_index('key')['Rank']; ranksP = P.set_index('key')['Rank']
    shift = (ranksR - ranksP.reindex(ranksR.index)).abs()
    log('Tie-tolerant check (score groups rounded to 9 decimals hold the same combinations): %s' % ('yes' if ok else 'NO'))
    log('Absolute rank difference per combination: max %d, number differing %d' % (shift.max(), int((shift > 0).sum())))
    rho = stats.spearmanr(M['MLSS_v4.0_R'], M['MLSS_v4.0_py'])[0]
    rho_rank = stats.spearmanr(R.set_index('key')['Rank'], P.set_index('key')['Rank'].reindex(R['key']))[0]
    log('Spearman rho of scores = %.6f; of the printed ranks = %.6f' % (rho, rho_rank))
    nties = int(R['s9'].duplicated(keep=False).sum())
    log('Number of combinations whose score is tied with another one in the R output: %d (ties keep the original row order in R, and in the port)' % nties)
    log('Top-15 sets identical: %s; top combination in R: %s (%.6f); in port: %s (%.6f)' % (
        set(seq_R[:15]) == set(seq_P[:15]), R.sort_values('Rank')['Combination'].iloc[0], R['MLSS_v4.0'].max(),
        P['Combination'].iloc[0], P['MLSS_v4.0'].max()))
    passed = (len(M) == 153 == len(R)) and d.max() < 0.005 and ok and shift.max() == 0
    log('GATE (153 scores within 0.005 and rank order reproduced): %s' % ('PASS' if passed else 'FAIL'))
    if out_csv:
        M[['Rank_R', 'Rank_py', 'Combination_R', 'MLSS_v4.0_R', 'MLSS_v4.0_py']].assign(abs_diff=d).sort_values('Rank_R').to_csv(out_csv, index=False)
    return passed, P, R


def gate_upstream(lines):
    """Part 2 (not required by the protocol): the Python network step against the stored R files and gapdh_exact.json."""
    def log(s=''):
        lines.append(s); print(s)
    log(''); log('Part 2: upstream network step (INS-1 matrix -> four schemes -> consensus -> degree) against stored results')
    X, src = load_matrix(); log('INS-1 matrix: %s (%d x %d)' % (src, X.shape[0], X.shape[1]))
    X17 = pd.read_excel(os.path.join(RDIR, 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx'))[GENES17].astype(float)
    log('  max |difference| to the R_17genes copy of the matrix: %.2e' % np.abs(X.values - X17.values).max())
    net = network_from_matrix(X)
    cons_R = {frozenset(s.split('_')) for s in pd.read_csv(os.path.join(RDIR, 'Fig_4d_consensus_pairs.csv'))['Consensus_Pairs']}
    log('  consensus pairs equal to Fig_4d_consensus_pairs.csv: %s (%d vs %d)' % ({frozenset(c) for c in net['consensus']} == cons_R, len(net['consensus']), len(cons_R)))
    hub_df = pd.read_csv(os.path.join(RDIR, 'Fig_6b_Gene_Centrality_Summary.csv'))
    log('  degree table equal to Fig_6b_Gene_Centrality_Summary.csv: %s' % (dict(zip(hub_df['Gene'], hub_df['Degree'])) == net['degree']))
    m = pd.read_csv(os.path.join(RDIR, 'Fig_4e_Hes1_significant_correlations_only.csv'), index_col=0)
    h = hes1_lookup(net); d = max(abs(h[frozenset((a.upper(), b.upper()))] - m.loc[a, b]) for a in m.index for b in m.columns if a != b)
    log('  Hes1-normalized r (120 pairs, all p < 0.05) vs Fig_4e CSV: max |difference| = %.2e; pairs with p >= 0.05 (NA): %d' % (d, 120 - len(h)))
    J = json.load(open(GAPDH_JSON))
    for w in ('reported', 'C0', 'C1'):
        n = network_from_matrix(gapdh_variant(X, w)); S = J[w]
        log('  %-8s consensus %2d (gapdh_exact.json %2d, identical %s); degree identical %s; significant pairs %s' % (
            w, len(n['consensus']), len(S['consensus']), {tuple(c) for c in n['consensus']} == {tuple(sorted(c)) for c in S['consensus']},
            n['degree'] == S['degree'], n['sig_counts']))
    return lines


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    lines = []
    ok, _, _ = gate(lines, os.path.join(HERE, 'gate_scores.csv'))     # part 1 = the protocol gate
    gate_upstream(lines)                                               # part 2 = upstream network step
    with open(os.path.join(HERE, 'gate_output.txt'), 'w', encoding='utf8') as f:
        f.write('\n'.join(lines) + '\n')
    sys.exit(0 if ok else 1)
