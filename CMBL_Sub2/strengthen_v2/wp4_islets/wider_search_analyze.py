"""wider_search_analyze.py - B4 and B5 of PRESPEC_WP4b_positive_control_and_wider_search.md for the kept series of wider_search/frozen_candidates.csv (fixed and fingerprinted before any value of the new series was read).
Per analysable series: log2 fold change of the 17 genes (perturbed arm minus control arm as frozen; arrays on the log2 scale of the series matrix; RNA-seq log2(CPM + 1) of NCBI-generated counts), multiplied by the expected
sign (-1 for MIRROR designs), compared with the INS-1 log2 fold change of the matching arm (../wp5_reporting/effect_sizes.csv, RT-qPCR, 'as reported'; GSE60803 from perturbation/gse60803_fold_changes.csv, contrast of addendum G):
genes with both nonzero, exact binomial test of the equal signs, Spearman rho. Series whose data cannot be obtained by plain download (K4) are listed as not analysable. Pooling (B5): REML + Hartung-Knapp of the
Fisher z of rho (variance 1.06 / (n - 3)); outcome words: CI entirely above 0 "concordant", entirely below 0 "discordant", else "inconclusive"; fewer than 3 series: not applied. Output: wider_search/fold_changes_<acc>.csv,
concordance_results.csv, concordance_meta.csv, concordance_report.md.   Usage: python -I wider_search_analyze.py"""
import os, sys, gzip, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd
from scipy import stats
import wp4_http as H
import wp4_core as C
import meta_analysis as MA
from attempts_log import record

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'wider_search'); ROOT = os.path.join(os.path.dirname(C.PARSED), 'wider')
MOUSE = {'Hes1': 'Hes1', 'Wnt5a': 'Wnt5a', 'Wnt5b': 'Wnt5b', 'Wnt9a': 'Wnt9a', 'Lef1': 'Lef1', 'Tcf7': 'Tcf7', 'Tcf7l2': 'Tcf7l2', 'Itpr1': 'Itpr1', 'Ryr2': 'Ryr2', 'Ryr3': 'Ryr3', 'Cacna1c': 'Cacna1c', 'Cacna1d': 'Cacna1d',
         'Kcnj11': 'Kcnj11', 'Ptbp1': 'Ptbp1', 'Glut2': 'Slc2a2', 'Ins1': 'Ins1', 'Glp1r': 'Glp1r'}
F = pd.read_csv(os.path.join(OUT, 'frozen_candidates.csv')); K = F[F.decision == 'keep'].set_index('accession')


def ins1():
    E = pd.read_csv(os.path.join(HERE, '..', 'wp5_reporting', 'effect_sizes.csv')); E = E[(E.readout == 'RT-qPCR') & (E.variant == 'as reported')]; out = {}
    for arm, con in (('DAPT', 'DAPT vs CTRL'), ('DKK-1', 'DKK-1 vs CTRL')): s = E[E.contrast == con].set_index('feature'); out[arm] = (s['log2fc_geo'].reindex(C.GENES), s['bh_q_17genes'].reindex(C.GENES))
    return out


def fc_gse60803(): d = pd.read_csv(os.path.join(HERE, 'perturbation', 'gse60803_fold_changes.csv')).set_index('gene'); return d['log2fc_all_donors'].reindex(C.GENES)


def titles(acc): return pd.read_csv(os.path.join(OUT, f'sample_records_{acc}.tsv'), sep='\t').set_index('gsm')['title']


def fc_gse99503():
    d = os.path.join(ROOT, 'GSE99503'); p = os.path.join(d, 'GSE99503_counts_GRCh38.p13_NCBI.tsv.gz'); url = 'https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE99503&format=file&file=GSE99503_raw_counts_GRCh38.p13_NCBI.tsv.gz'
    n = H.download(url, p, max_mb=50); record('GSE99503', 'NCBI-generated raw counts GRCh38.p13', url, n, 'NCBI RNA-seq counts', 'downloaded', 'size checked; contrast frozen before reading')
    X = pd.read_csv(p, sep='\t', index_col=0); T = titles('GSE99503'); ctrl = [g for g, t in T.items() if re.fullmatch(r'Neg siRNA #\d', t)]; pert = [g for g, t in T.items() if re.fullmatch(r'Neg siRNA \+ 1-AKP #\d', t)]
    assert len(ctrl) == 3 and len(pert) == 3 and set(ctrl + pert) <= set(X.columns), (ctrl, pert)
    L = np.log2(X[ctrl + pert] / X[ctrl + pert].sum(axis=0) * 1e6 + 1.0); fc = L[pert].mean(axis=1) - L[ctrl].mean(axis=1)
    return pd.Series({g: (float(fc.get(C.ENTREZ[g], np.nan))) for g in C.GENES}), dict(genes_found=int(sum(C.ENTREZ[g] in fc.index for g in C.GENES)), ctrl=len(ctrl), pert=len(pert))


def fc_gse167073():
    d = os.path.join(ROOT, 'GSE167073'); T = titles('GSE167073'); ctrl = [g for g, t in T.items() if 'Fzd9 wt' in t]; pert = [g for g, t in T.items() if 'Fzd9 ko' in t]; assert len(ctrl) == 4 and len(pert) == 4
    tgt = {v: k for k, v in MOUSE.items()}; pm = {}
    with open(os.path.join(d, 'GPL17400_data.txt'), encoding='utf8', errors='replace') as f:
        hdr = None
        for l in f:
            if hdr is None:
                if l.startswith('ID\t'): hdr = l.rstrip('\n').split('\t'); ig = hdr.index('gene_assignment')
                continue
            p = l.rstrip('\n').split('\t')
            if len(p) <= ig or p[ig] in ('---', ''): continue
            for a in p[ig].split(' /// '):
                q = a.split(' // ')
                if len(q) > 1 and q[1].strip() in tgt: pm.setdefault(p[0], set()).add(q[1].strip())
    vals = {}
    with gzip.open(os.path.join(d, 'GSE167073_series_matrix.txt.gz'), 'rt', encoding='utf8', errors='replace') as f:
        inside, cols = False, None
        for l in f:
            if l.startswith('!series_matrix_table_begin'): inside = True; continue
            if l.startswith('!series_matrix_table_end'): break
            if not inside: continue
            p = l.rstrip('\n').split('\t')
            if cols is None: cols = [x.strip('"') for x in p[1:]]; continue
            pid = p[0].strip('"')
            if pid in pm: vals[pid] = np.array([float(x) for x in p[1:]])
    ci = [cols.index(g) for g in ctrl]; pi = [cols.index(g) for g in pert]; out = {}; probes = {}
    for g in C.GENES:
        cand = [pid for pid, s in pm.items() if MOUSE[g] in s and pid in vals]
        if not cand: out[g] = np.nan; continue
        best = max(cand, key=lambda x: vals[x][ci].mean()); v = vals[best]; assert v.max() < 30, 'matrix not on a log2 scale'
        out[g] = float(v[pi].mean() - v[ci].mean()); probes[g] = (best, len(cand))
    return pd.Series(out).reindex(C.GENES), dict(genes_found=int(sum(not np.isnan(x) for x in out.values())), probes=probes)


def concord(fc, sign, vec):
    x = (fc * sign).reindex(C.GENES); D = pd.DataFrame({'fc': x, 'ins1': vec}).dropna(); D = D[(D.fc != 0) & (D.ins1 != 0)]; n = len(D)
    if n < 3: return dict(n_genes=n)
    k = int((np.sign(D.fc) == np.sign(D.ins1)).sum()); rho, p = stats.spearmanr(D.fc, D.ins1)
    return dict(n_genes=n, n_same_sign=k, binomial_p=float(stats.binomtest(k, n, 0.5).pvalue), spearman_rho=float(rho), spearman_p=float(p), same_sign_genes=D.index[np.sign(D.fc) == np.sign(D.ins1)].tolist())


if __name__ == '__main__':
    V = ins1(); res = []; fcs = {}
    for acc in K.index:
        k = K.loc[acc]; base = dict(accession=acc, arm=k.arm, expected_sign=int(k.expected_sign), mirror=int(k.mirror), weak_match_note=k.weak_match_note, title=F[F.accession == acc].title.iloc[0], k4=k.k4)
        if str(k.k4).startswith('NOT'): res.append(dict(base, status='not analysable (K4): ' + k.data_route)); print(acc, 'not analysable'); continue
        try:
            if acc == 'GSE60803': fc, info = fc_gse60803(), {}
            elif acc == 'GSE99503': fc, info = fc_gse99503()
            elif acc == 'GSE167073': fc, info = fc_gse167073()
        except Exception as e: res.append(dict(base, status='not analysable (data step failed: %s)' % repr(e)[:150])); print(acc, 'FAILED', repr(e)[:200]); continue
        pd.DataFrame({'gene': C.GENES, 'log2fc_perturbed_minus_control': fc.values, 'expected_sign': int(k.expected_sign), 'log2fc_signed': (fc * int(k.expected_sign)).values, 'ins1_log2fc_' + k.arm: V[k.arm][0].values}).to_csv(os.path.join(OUT, f'fold_changes_{acc}.csv'), index=False)
        r = concord(fc, int(k.expected_sign), V[k.arm][0]); strong = V[k.arm][1][V[k.arm][1] < 0.05].index.tolist(); sg = (fc * int(k.expected_sign)).reindex(strong).dropna(); sg = sg[sg != 0]
        r.update(base); r['status'] = 'analysed'; r['genes_found'] = info.get('genes_found', int(fc.notna().sum())); r['strong_ins1_genes'] = len(sg); r['strong_same_sign'] = int((np.sign(sg) == np.sign(V[k.arm][0].reindex(sg.index))).sum()); res.append(r)
        print('%-10s %-6s n_genes %2d same sign %2d P %.3f rho %+.3f (P %.3f) | strong INS-1 genes same sign %d of %d' % (acc, k.arm, r.get('n_genes', 0), r.get('n_same_sign', -1), r.get('binomial_p', np.nan), r.get('spearman_rho', np.nan), r.get('spearman_p', np.nan), r['strong_same_sign'], r['strong_ins1_genes']), flush=True)
    R = pd.DataFrame(res); R.to_csv(os.path.join(OUT, 'concordance_results.csv'), index=False)
    A = R[R.status == 'analysed'].copy(); A = A[A.n_genes >= 12] if len(A) else A; metas = []
    def meta(label, S):
        if len(S) < 3: return dict(set=label, k=len(S), applied=False, outcome='rule not applied (fewer than 3 analysable series)')
        z = np.arctanh(np.clip(S.spearman_rho.to_numpy(float), -0.999, 0.999)); se = np.sqrt(1.06 / (S.n_genes.to_numpy(float) - 3)); r = MA.pool(z, se)
        lo, hi = np.tanh(r['lo']), np.tanh(r['hi']); return dict(set=label, k=len(S), applied=True, pooled_rho=float(np.tanh(r['mu'])), lo=float(lo), hi=float(hi), I2=r['I2'], tau2=r['tau2'], outcome='concordant' if lo > 0 else 'discordant' if hi < 0 else 'inconclusive')
    metas.append(meta('PRIMARY: all analysable series', A)); metas.append(meta('Wnt series (DKK-1 arm)', A[A.arm == 'DKK-1'])); metas.append(meta('Notch series (DAPT arm)', A[A.arm == 'DAPT']))
    metas.append(meta('mirror designs excluded', A[A.mirror == 0])); metas.append(meta('GSE60803 excluded', A[A.accession != 'GSE60803']))
    M = pd.DataFrame(metas); M.to_csv(os.path.join(OUT, 'concordance_meta.csv'), index=False)
    L = ['# Wider search: concordance of kept perturbation series with the INS-1 arms (rules fixed in PRESPEC_WP4b... B4/B5 before any value was read)', '', '| series | arm | expected sign | status | genes | same sign | exact binomial P | Spearman rho (P) | strong INS-1 genes same sign |', '|---|---|---|---|---|---|---|---|---|']
    for r in R.itertuples():
        if r.status == 'analysed': L.append('| %s | %s | %+d%s | analysed | %d | %d | %.3f | %+.3f (%.3f) | %d of %d |' % (r.accession, r.arm, r.expected_sign, ' (mirror)' if r.mirror else '', r.n_genes, r.n_same_sign, r.binomial_p, r.spearman_rho, r.spearman_p, r.strong_same_sign, r.strong_ins1_genes))
        else: L.append('| %s | %s | %+d%s | %s | - | - | - | - | - |' % (r.accession, r.arm, r.expected_sign, ' (mirror)' if r.mirror else '', r.status[:110]))
    L += ['', '## Pooled (REML + Hartung-Knapp on Fisher z)', '', '| set | k | pooled rho | 95% CI | I2 | outcome |', '|---|---|---|---|---|---|']
    for m in M.itertuples(): L.append('| %s | %d | %s | %s | %s | %s |' % (m.set, m.k, ('%+.3f' % m.pooled_rho) if m.applied else '-', ('[%+.3f, %+.3f]' % (m.lo, m.hi)) if m.applied else '-', ('%.0f%%' % (100 * m.I2)) if m.applied else '-', m.outcome))
    open(os.path.join(OUT, 'concordance_report.md'), 'w', encoding='utf8').write('\n'.join(L) + '\n'); print('\n'.join(L))
