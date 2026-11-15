"""perturbation_gse60803.py - protocol 4.6, step 3: log2 fold changes of the 17 genes in GSE60803 (adult human beta-cell-derived cells, beta-catenin shRNA; Affymetrix HuGene-1.1 ST, GPL11532) and their agreement in sign
with the DKK-1 arm of the INS-1 experiment. Design and contrast fixed in PRESPEC_WP4_choices.md addendum G BEFORE any expression value was read:
  per donor: control = the untreated array (S108, S109, S112) or the nontarget array (S114); perturbed = mean of the donor's beta-catenin shRNA arrays (9 and 45; both replicates of S114);
  log2 fold change of a gene = mean over the 4 donors of (perturbed - control) on the log2 scale of the series matrix; the probe set of a gene (Entrez id match in the GPL11532 annotation) with the highest
  mean control expression is used. Secondary, descriptive: donor S114 only (3 beta-catenin arrays against the single nontarget array).
INS-1 side: log2 fold change DKK-1 vs CTRL, 'as reported', RT-qPCR panel, from ../wp5_reporting/effect_sizes.csv. Agreement: genes with both fold changes nonzero, exact two-sided binomial test of the number of
equal signs against 0.5, Spearman correlation of the two vectors. DESCRIPTIVE: the controls are not expanded cells, the model systems differ (rat cell line against human cells), and the sample size is 17 genes at best.
Output: perturbation/gse60803_fold_changes.csv, perturbation/gse60803_concordance.json.   Usage: python -I perturbation_gse60803.py"""
import os, sys, gzip, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd
from scipy import stats
import wp4_core as C

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(os.path.dirname(C.PARSED), 'GSE60803'); OUT = os.path.join(HERE, 'perturbation')
# series matrix
meta, rows = {}, []
with gzip.open(os.path.join(ROOT, 'GSE60803_series_matrix.txt.gz'), 'rt', encoding='utf8', errors='replace') as f:
    inside = False
    for l in f:
        if l.startswith('!series_matrix_table_begin'): inside = True; continue
        if l.startswith('!series_matrix_table_end'): break
        if inside: rows.append(l.rstrip('\n').split('\t'))
        elif l.startswith('!Sample_title'): meta['title'] = [x.strip('"') for x in l.rstrip('\n').split('\t')[1:]]
        elif l.startswith('!Sample_geo_accession'): meta['gsm'] = [x.strip('"') for x in l.rstrip('\n').split('\t')[1:]]
X = pd.DataFrame(rows[1:], columns=[c.strip('"') for c in rows[0]]).set_index('ID_REF').astype(float); X.index = X.index.str.strip('"')
assert list(X.columns) == meta['gsm'], 'sample order differs'
S = pd.DataFrame({'gsm': meta['gsm'], 'title': meta['title']}); S['donor'] = S.title.str.extract(r'donor (S\d+)')[0]
S['arm'] = np.where(S.title.str.contains('untreated'), 'control', np.where(S.title.str.contains('nontarget'), 'control', np.where(S.title.str.contains(r'β-catenin'), 'bcat', 'unknown')))
assert (S.arm != 'unknown').all() and S.donor.nunique() == 4
print('expression range of the series matrix: %.2f to %.2f (median %.2f), probe sets %d' % (X.values.min(), X.values.max(), np.median(X.values), len(X)), flush=True)
assert X.values.max() < 30, 'the matrix is not on a log2 scale'
# annotation: Entrez id -> probe sets
ann = {}
with gzip.open(os.path.join(ROOT, 'GPL11532.annot.gz'), 'rt', encoding='utf8', errors='replace') as f:
    hdr = None
    for l in f:
        if l.startswith(('#', '!', '^')): continue
        p = l.rstrip('\n').split('\t')
        if hdr is None: hdr = p; iG, iS = hdr.index('Gene ID'), hdr.index('Gene symbol'); continue
        for gid in re.split(r'\s*///\s*', p[iG]):
            if gid.strip().isdigit(): ann.setdefault(int(gid), []).append(p[0])
ctrl = S[S.arm == 'control'].set_index('donor'); pert = S[S.arm == 'bcat']
res = []
for g in C.GENES:
    probes = [p for p in ann.get(C.ENTREZ[g], []) if p in X.index]
    if not probes: res.append(dict(gene=g, probe=None)); continue
    cm = X.loc[probes, ctrl.gsm].mean(axis=1); best = cm.idxmax(); d = []
    for donor in sorted(S.donor.unique()):
        d.append(X.loc[best, pert[pert.donor == donor].gsm].mean() - X.loc[best, ctrl.loc[donor, 'gsm']])
    s114 = X.loc[best, pert[pert.donor == 'S114'].gsm].mean() - X.loc[best, ctrl.loc['S114', 'gsm']]
    res.append(dict(gene=g, probe=best, probes_for_gene=len(probes), mean_control_log2=float(cm[best]), **{'delta_' + k: float(v) for k, v in zip(sorted(S.donor.unique()), d)}, log2fc_all_donors=float(np.mean(d)), sd_across_donors=float(np.std(d, ddof=1)), log2fc_S114_only=float(s114)))
R = pd.DataFrame(res)
E = pd.read_csv(os.path.join(HERE, '..', 'wp5_reporting', 'effect_sizes.csv'))
E = E[(E.readout == 'RT-qPCR') & (E.contrast == 'DKK-1 vs CTRL') & (E.variant == 'as reported')].set_index('feature')['log2fc_geo']
R['ins1_log2fc_DKK1'] = R.gene.map(E); R['sign_agree_all'] = np.sign(R.log2fc_all_donors) == np.sign(R.ins1_log2fc_DKK1); R['sign_agree_S114'] = np.sign(R.log2fc_S114_only) == np.sign(R.ins1_log2fc_DKK1)
R.to_csv(os.path.join(OUT, 'gse60803_fold_changes.csv'), index=False)
out = dict(series='GSE60803', platform='GPL11532', donors=sorted(S.donor.unique()), arrays=len(S), genes_with_probe=int(R.probe.notna().sum()), genes_without_probe=R[R.probe.isna()].gene.tolist(), ins1_arm='DKK-1 vs CTRL, as reported', confounded='controls are not expanded cells (untreated and nontarget arrays are described as intact islet cells; nontarget exists for one donor)')
for tag, col, agree in (('all_donors', 'log2fc_all_donors', 'sign_agree_all'), ('S114_only', 'log2fc_S114_only', 'sign_agree_S114')):
    D = R.dropna(subset=[col, 'ins1_log2fc_DKK1']); D = D[(D[col] != 0) & (D.ins1_log2fc_DKK1 != 0)]; k, n = int(D[agree].sum()), len(D)
    rho, p = stats.spearmanr(D[col], D.ins1_log2fc_DKK1) if n > 2 else (np.nan, np.nan)
    out[tag] = dict(n_genes=n, n_same_sign=k, binomial_p=float(stats.binomtest(k, n, 0.5).pvalue) if n else None, spearman_rho=float(rho), spearman_p=float(p), genes=D.gene.tolist(), same_sign_genes=D[D[agree]].gene.tolist(), opposite_sign_genes=D[~D[agree]].gene.tolist())
json.dump(out, open(os.path.join(OUT, 'gse60803_concordance.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
pd.set_option('display.width', 250); print(R.round(3).to_string(index=False)); print(json.dumps({k: v for k, v in out.items() if k in ('all_donors', 'S114_only')}, indent=1))
