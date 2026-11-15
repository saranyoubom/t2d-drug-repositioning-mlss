"""LEAN COPY of ../../human_islet_GSE81608/human_islet_check.py for the WP4 gate (6 Nov 2026). Everything after the loading block is the original text, line for line; only the loading block differs:
the original reads the whole RPKM table (39851 genes x 1600 cells, about 510 MB as float64, three to four copies at peak) into memory, which this 8 GB machine with about 1-2 GB free should not do.
Here the same pandas C parser reads the file in chunks of rows and keeps only the 503 beta-cell columns (and the INS / GCG rows of all cells for the identity check), so the beta-cell values are the same numbers.
Cross-check of the original (Xin et al. 2016; 1,600 cells, 18 donors; RPKM) as in the original docstring: (1) expression and T2D vs non-diabetic (donor level); (2) co-expression of the 21 consensus pairs within beta-cells (donor
effects removed) against the other pairs and expression-matched random pairs; (3) hubs; (4) the Cacna1d-Glp1r pair.
Usage: python -I gate_human_islet_check_lean.py <dir with series_matrix.txt.gz and human_islets_rpkm.txt.gz> <R_17genes dir> <output dir>   (seed 2026)"""
import sys, os, gzip, json
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
sys.stdout.reconfigure(encoding='utf8')
D, R17, OUT = sys.argv[1], sys.argv[2], sys.argv[3]; os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(2026)
ENTREZ = {'Hes1': 3280, 'Wnt5a': 7474, 'Wnt5b': 81029, 'Wnt9a': 7483, 'Lef1': 51176, 'Tcf7': 6932, 'Tcf7l2': 6934, 'Itpr1': 3708, 'Ryr2': 6262, 'Ryr3': 6263,
          'Cacna1c': 775, 'Cacna1d': 776, 'Kcnj11': 3767, 'Ptbp1': 5725, 'Glut2': 6514, 'Ins1': 3630, 'Glp1r': 2740}      # human orthologs (Glut2 = SLC2A2, Ins1 = INS)
HUMAN = {'Glut2': 'SLC2A2', 'Ins1': 'INS'}
sym = lambda g: HUMAN.get(g, g.upper())
# ---------------- metadata and matrix
rows = {}
with gzip.open(os.path.join(D, 'series_matrix.txt.gz'), 'rt', encoding='utf8', errors='replace') as f:
    for l in f:
        if l.startswith(('!Sample_title', '!Sample_characteristics_ch1')): rows.setdefault(l.split('\t')[0], []).append([x.strip('"') for x in l.rstrip('\n').split('\t')[1:]])
ch = rows['!Sample_characteristics_ch1']; title = rows['!Sample_title'][0]
col = lambda i: [x.split(': ', 1)[1] for x in ch[i]]
meta = pd.DataFrame({'sample': ['Sample_' + t.split()[-1] for t in title], 'donor': col(1), 'condition': col(2), 'type': col(6)})
# ---- LEAN LOADING BLOCK (replaces: X = pd.read_csv(...); L = np.log2(X + 1.0); tmeans = ...)
ins, gcg = ENTREZ['Ins1'], 2641
b = meta[meta.type == 'beta'].reset_index(drop=True)
sel = {t: (meta.type == t).values for t in ('alpha', 'beta')}
keep, ident = [], {}
for ch_ in pd.read_csv(os.path.join(D, 'human_islets_rpkm.txt.gz'), sep='\t', index_col=0, chunksize=2500):
    if len(keep) == 0: assert list(ch_.columns) == list(meta['sample']), 'column order differs from the metadata'
    for gid in (ins, gcg):
        if gid in ch_.index: ident[gid] = {t: float(np.log2(ch_.loc[gid].values.astype(float) + 1.0)[sel[t]].mean()) for t in sel}
    keep.append(ch_[list(b['sample'])].astype('float64'))
Xb = pd.concat(keep); del keep                                                      # genes x beta cells, same float64 values as X.loc[:, b['sample']] of the original
Lb = np.log2(Xb + 1.0)
print('sanity: INS beta/alpha mean log2 %.2f / %.2f | GCG alpha/beta %.2f / %.2f' % (ident[ins]['beta'], ident[ins]['alpha'], ident[gcg]['alpha'], ident[gcg]['beta']))
assert ident[ins]['beta'] > ident[ins]['alpha'] + 1 and ident[gcg]['alpha'] > ident[gcg]['beta'] + 1, 'gene id mapping failed the INS/GCG check'
print('beta cells:', len(b), '| per donor:', b.groupby(['condition', 'donor']).size().to_dict())
det = (Xb > 0).mean(axis=1)
# ---- END OF THE LEAN LOADING BLOCK; everything below is the original text
# ---------------- 1. expression and T2D vs non-diabetic (donor level)
don = b.groupby('donor').size(); ok = don[don >= 10].index
res1 = []
for g, e in ENTREZ.items():
    y = Lb.loc[e]; pb = y.groupby(b['donor'].values).mean().loc[ok]; cond = b.drop_duplicates('donor').set_index('donor').loc[ok, 'condition']
    t2d, nd = pb[cond == 'T2D'], pb[cond == 'non-diabetic']
    p = stats.mannwhitneyu(t2d, nd, method='exact').pvalue if len(t2d) > 1 and len(nd) > 1 else np.nan
    pc = stats.mannwhitneyu(y[(b.condition == 'T2D').values], y[(b.condition == 'non-diabetic').values]).pvalue
    res1.append(dict(gene_invitro=g, gene_human=sym(g), entrez=e, beta_cells=len(b), detected_pct=100 * det[e], median_log2_detected=float(np.median(y[y > 0])) if (y > 0).any() else np.nan,
                     donors_T2D=len(t2d), donors_nondiabetic=len(nd), mean_log2_T2D=t2d.mean(), mean_log2_ND=nd.mean(), diff_T2D_minus_ND=t2d.mean() - nd.mean(), p_donor_level=p, p_cell_level_not_independent=pc))
r1 = pd.DataFrame(res1); r1['q_donor_level'] = multipletests(r1.p_donor_level.fillna(1), method='fdr_bh')[1]; r1.to_csv(os.path.join(OUT, 'genes_expression_T2D_vs_ND.csv'), index=False)
# ---------------- 2. co-expression within beta cells (donor effects removed)
elig = det[det >= 0.10].index                                                       # genes detected in >= 10% of beta cells (the matching pool)
G = Lb.loc[elig]; Gc = G.sub(G.T.groupby(b['donor'].values).transform('mean').T)   # donor-centred
Rk = Gc.rank(axis=1).values; Rk = (Rk - Rk.mean(1, keepdims=True)); Rk = Rk / (np.sqrt((Rk ** 2).sum(1, keepdims=True)) + 1e-12)   # rank-standardised; dot product = Spearman rho
idx = {g: i for i, g in enumerate(G.index)}
n = Rk.shape[1]; dfree = n - b.donor.nunique() - 2
rho = lambda a, c: float(Rk[idx[a]] @ Rk[idx[c]])
cons = [tuple(x.split('_')) for x in pd.read_csv(os.path.join(R17, 'Fig_4d_consensus_pairs.csv')).iloc[:, 0]]
edge = pd.read_csv(os.path.join(R17, 'Fig_6b_Edge_Correlation_Statistics.csv')); redge = {frozenset((r['Gene1'], r['Gene2'])): r['Correlation (r)'] for _, r in edge.iterrows()}
consset = {frozenset(c) for c in cons}; genes = list(ENTREZ)
pairs = []
for i, a in enumerate(genes):
    for c in genes[i + 1:]:
        ea, ec = ENTREZ[a], ENTREZ[c]; evaluable = (ea in idx) and (ec in idx)
        r_ = rho(ea, ec) if evaluable else np.nan
        t_ = r_ * np.sqrt(dfree / (1 - r_ ** 2)) if evaluable else np.nan
        pairs.append(dict(gene1=a, gene2=c, consensus_invitro=frozenset((a, c)) in consset, r_invitro_Hes1norm=redge.get(frozenset((a, c)), np.nan), evaluable=evaluable,
                          detected_pct_1=100 * det[ea], detected_pct_2=100 * det[ec], rho_human_beta_donor_adjusted=r_, p=2 * stats.t.sf(abs(t_), dfree) if evaluable else np.nan))
P = pd.DataFrame(pairs); ev = P[P.evaluable].copy(); ev['q'] = multipletests(ev.p, method='fdr_bh')[1]; P = P.merge(ev[['gene1', 'gene2', 'q']], on=['gene1', 'gene2'], how='left')
P.to_csv(os.path.join(OUT, 'pairs_136_human_beta_coexpression.csv'), index=False)
c_ev, o_ev = ev[ev.consensus_invitro], ev[~ev.consensus_invitro]
summ = dict(beta_cells=int(len(b)), donors=int(b.donor.nunique()), df_for_p=int(dfree), genes_evaluable=int(sum(ENTREZ[g] in idx for g in genes)), genes_not_evaluable=[g for g in genes if ENTREZ[g] not in idx],
            consensus_pairs_total=len(cons), consensus_pairs_evaluable=int(len(c_ev)), other_pairs_evaluable=int(len(o_ev)),
            median_rho_consensus=float(c_ev.rho_human_beta_donor_adjusted.median()), median_rho_other=float(o_ev.rho_human_beta_donor_adjusted.median()),
            mannwhitney_p_consensus_vs_other=float(stats.mannwhitneyu(c_ev.rho_human_beta_donor_adjusted, o_ev.rho_human_beta_donor_adjusted).pvalue),
            consensus_positive=int((c_ev.rho_human_beta_donor_adjusted > 0).sum()), consensus_sig_positive_q05=int(((c_ev.rho_human_beta_donor_adjusted > 0) & (c_ev.q < 0.05)).sum()),
            consensus_sig_negative_q05=int(((c_ev.rho_human_beta_donor_adjusted < 0) & (c_ev.q < 0.05)).sum()),
            binomial_p_positive_vs_half=float(stats.binomtest(int((c_ev.rho_human_beta_donor_adjusted > 0).sum()), len(c_ev), 0.5).pvalue))
# expression-matched null: random gene pairs from the beta-cell transcriptome with the same detection-rate deciles as each consensus pair
bins = pd.qcut(det[elig].rank(method='first'), 10, labels=False); own = {ENTREZ[g] for g in genes}
pool = {k: np.array([i for g, i in idx.items() if bins[g] == k and g not in own]) for k in range(10)}
obs_med = float(c_ev.rho_human_beta_donor_adjusted.median()); obs_mean = float(c_ev.rho_human_beta_donor_adjusted.mean())
ev_pairs = [(r.gene1, r.gene2) for r in c_ev.itertuples()]
NN = 2000; meds = np.empty(NN); means = np.empty(NN)
for s in range(NN):
    vals = []
    for a, c in ev_pairs:
        ia = rng.choice(pool[bins[ENTREZ[a]]]); ic = rng.choice(pool[bins[ENTREZ[c]]])
        while ic == ia: ic = rng.choice(pool[bins[ENTREZ[c]]])
        vals.append(float(Rk[ia] @ Rk[ic]))
    meds[s], means[s] = np.median(vals), np.mean(vals)
summ.update(null_median_rho_mean=float(meds.mean()), null_median_rho_95pct=float(np.percentile(meds, 95)), observed_median_rho=obs_med,
            empirical_p_median_ge_observed=float((1 + (meds >= obs_med).sum()) / (NN + 1)), empirical_p_mean_ge_observed=float((1 + (means >= obs_mean).sum()) / (NN + 1)))
# pair-level null (added 4 Nov 2026 with a separate generator, so the results above do not change): percentile of each evaluable consensus pair among 2,000 random pairs with the same detection deciles
rng2 = np.random.default_rng(2027); pn = []
for r in c_ev.itertuples():
    a, c = r.gene1, r.gene2; vals = np.empty(NN)
    for k in range(NN):
        ia = rng2.choice(pool[bins[ENTREZ[a]]]); ic = rng2.choice(pool[bins[ENTREZ[c]]])
        while ic == ia: ic = rng2.choice(pool[bins[ENTREZ[c]]])
        vals[k] = Rk[ia] @ Rk[ic]
    pn.append(dict(gene1=a, gene2=c, rho=r.rho_human_beta_donor_adjusted, q=r.q, null_mean=vals.mean(), null_95pct=np.percentile(vals, 95), empirical_p_ge_observed=(1 + (vals >= r.rho_human_beta_donor_adjusted).sum()) / (NN + 1)))
pd.DataFrame(pn).to_csv(os.path.join(OUT, 'consensus_pairs_vs_matched_null.csv'), index=False)
# ---------------- 3. hubs: degree in the human beta-cell network (q < 0.05 and |rho| >= 0.2) vs in vitro degree
cen = pd.read_csv(os.path.join(R17, 'Fig_6b_Gene_Centrality_Summary.csv')).set_index('Gene')['Degree']
E = ev[(ev.q < 0.05) & (ev.rho_human_beta_donor_adjusted.abs() >= 0.2)]
hd = {g: int(((E.gene1 == g) | (E.gene2 == g)).sum()) for g in genes}
H = pd.DataFrame({'gene': genes, 'degree_invitro_consensus': [int(cen.get(g, 0)) for g in genes], 'degree_human_beta': [hd[g] for g in genes], 'detected_pct_beta': [100 * det[ENTREZ[g]] for g in genes]})
H.to_csv(os.path.join(OUT, 'hubs_invitro_vs_human.csv'), index=False)
sp = stats.spearmanr(H.degree_invitro_consensus, H.degree_human_beta)
summ.update(human_edges_q05_rho02=int(len(E)), spearman_degree_invitro_vs_human=float(sp.statistic), spearman_degree_p=float(sp.pvalue), human_top3_degree=H.sort_values('degree_human_beta', ascending=False).head(3)[['gene', 'degree_human_beta']].values.tolist(),
            invitro_hubs_human_degree={g: hd[g] for g in ('Glp1r', 'Tcf7l2', 'Cacna1d', 'Ptbp1')})
pr = P[(P.gene1 == 'Glp1r') & (P.gene2 == 'Cacna1d') | (P.gene1 == 'Cacna1d') & (P.gene2 == 'Glp1r')].iloc[0]
summ['Cacna1d_Glp1r'] = dict(rho=float(pr.rho_human_beta_donor_adjusted) if pr.evaluable else None, p=float(pr.p) if pr.evaluable else None, evaluable=bool(pr.evaluable), detected_pct_Glp1r=float(100 * det[ENTREZ['Glp1r']]), detected_pct_Cacna1d=float(100 * det[ENTREZ['Cacna1d']]))
json.dump(summ, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)
print(json.dumps(summ, indent=1)); print(r1[['gene_human', 'detected_pct', 'diff_T2D_minus_ND', 'p_donor_level', 'q_donor_level']].round(3).to_string(index=False))
print(H.to_string(index=False))
