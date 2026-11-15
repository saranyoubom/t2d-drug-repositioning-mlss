"""gene_graph_null.py - ADD-ON, POST HOC (written after the primary meta-analysis result was known; not in the protocol): does the verdict of 4.4 survive when the null keeps the gene graph of the consensus pairs?
Why: the pre-specified null draws the two genes of every consensus pair independently, but the evaluable consensus pairs share genes (CACNA1D is in 5 of the 11, PTBP1 in 5, CACNA1C and GLP1R and TCF7L2 in 4 each), so the
mean of their correlations varies more than a mean over independent random pairs and the pre-specified SE of E is probably too small.
What: for each dataset (depth-adjusted analysis, all frozen beta-cells, same eligible genes, same detection deciles, same partial-Spearman ranks as wp4_core.analyse) every distinct gene of the evaluable consensus pairs is replaced
by ONE random gene of its own detection decile (not one of the 17, no gene used twice in a draw) and the mean Fisher z of the same pairs is computed on the replaced genes; 2000 draws, seed 2028 (new seed). E_graph = mean z of
the observed pairs minus the mean of the null; SE_graph = SD of the null means. The meta-analysis (same estimators and margin as meta_analysis.py) is repeated with E_graph and SE_graph.
Output: meta/gene_graph_by_dataset.csv, meta/gene_graph_meta.csv, meta/gene_graph_report.md.   Usage: python -I gene_graph_null.py
The result is reported beside the pre-specified one; the rule of 4.4 stays defined on the pre-specified SE."""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd
import wp4_core as C
import meta_analysis as MA
from run_datasets import load

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'meta'); SEED = 2028; NN = 2000
TAGS = ['gse81608', 'gse84133_human', 'gse86469', 'gse81547', 'emtab5061', 'gse84133_mouse']


def graph_null(cells, thr=0.10, nn=NN, seed=SEED):
    cons, consset, redge, cen = C.consensus_tables(); sym = cells.label
    det = np.asarray((cells.B > 0).mean(axis=0)).ravel(); elig = np.nonzero(det >= thr)[0]; labels_e = [cells.genes[j] for j in elig]; idx = {g: r for r, g in enumerate(labels_e)}
    bins, pool = C.make_pool(det[elig], labels_e, {sym(g) for g in C.GENES})
    Rk = C.adjusted_ranks(C.center_by_donor(C.lognorm(cells, elig), cells.donor), C.depth_design(cells.depth, cells.donor))
    pairs = [(a, c) for i, a in enumerate(C.GENES) for c in C.GENES[i + 1:] if frozenset((a, c)) in consset and sym(a) in idx and sym(c) in idx]
    genes = [g for g in C.GENES if any(g in p for p in pairs)]; gi = {g: k for k, g in enumerate(genes)}
    rho = np.array([float(Rk[idx[sym(a)]] @ Rk[idx[sym(c)]]) for a, c in pairs]); z_obs = float(np.arctanh(np.clip(rho, -0.999999, 0.999999)).mean())
    rng = np.random.default_rng(seed); ints = rng.integers; M = np.empty((nn, len(genes)), np.int64)
    for s in range(nn):
        used = set()
        for k, g in enumerate(genes):
            arr = pool[bins[idx[sym(g)]]]
            while True:
                x = int(arr[ints(0, len(arr))])
                if x not in used: break
            used.add(x); M[s, k] = x
    Z = np.zeros(nn)
    for a, c in pairs: Z += np.arctanh(np.clip(np.einsum('ij,ij->i', Rk[M[:, gi[a]]], Rk[M[:, gi[c]]]), -0.999999, 0.999999))
    Z /= len(pairs); share = {g: sum(g in p for p in pairs) for g in genes}
    return dict(pairs=len(pairs), genes_in_pairs=len(genes), max_pairs_per_gene=max(share.values()), z_obs=z_obs, z_null_graph=float(Z.mean()), E_graph=z_obs - float(Z.mean()), SE_graph=float(Z.std(ddof=1)))


if __name__ == '__main__':
    rows = []
    for tag in TAGS:
        t0 = time.time(); cells = load(tag); r = graph_null(cells); s = json.load(open(os.path.join(HERE, 'results', f'{tag}__primary__summary.json')))['adjusted']
        r.update(tag=tag, E_prespecified=s.get('E'), SE_prespecified=s.get('SE_E')); r['SE_ratio'] = r['SE_graph'] / r['SE_prespecified'] if r['SE_prespecified'] else np.nan; rows.append(r)
        print('%-15s pairs %2d genes %2d (max %d pairs per gene) | E %+.3f -> %+.3f | SE %.3f -> %.3f (x%.2f) | %.0f s' % (tag, r['pairs'], r['genes_in_pairs'], r['max_pairs_per_gene'], r['E_prespecified'], r['E_graph'], r['SE_prespecified'], r['SE_graph'], r['SE_ratio'], time.time() - t0), flush=True)
        del cells
    D = pd.DataFrame(rows); D.to_csv(os.path.join(OUT, 'gene_graph_by_dataset.csv'), index=False)
    H = D[D.tag.isin(MA.HUMAN) & (D.pairs >= 3)]; res = []
    L = ['# WP4 add-on (post hoc): meta-analysis with a null that keeps the gene graph of the consensus pairs', '', 'Same estimators and margin as meta_analysis.py; E_graph and SE_graph replace E and SE of the pre-specified analysis. The rule of 4.4 stays defined on the pre-specified SE; this is a sensitivity analysis.', '']
    for name, kw in (('REML+HK', dict(tau='reml', ci='hk', floor=False)), ('REML+HK, factor floored at 1', dict(tau='reml', ci='hk', floor=True)), ('DL+HK', dict(tau='dl', ci='hk')), ('DL, normal CI', dict(tau='dl', ci='normal')), ('fixed effect', dict(ci='fixed'))):
        r = MA.pool(H.E_graph.to_numpy(float), H.SE_graph.to_numpy(float), **kw); r.update(estimator=name, verdict=MA.verdict(r['lo'], r['hi']), analysis='gene-graph null'); res.append(r)
        L.append('* %-30s pooled E %+.3f  95%% CI [%+.3f, %+.3f]  tau2 %.4f  I2 %.0f%%  -> %s' % (name, r['mu'], r['lo'], r['hi'], r['tau2'], 100 * r['I2'], r['verdict']))
    L.append('')
    for t in H.tag:
        R = H[H.tag != t]; r = MA.pool(R.E_graph.to_numpy(float), R.SE_graph.to_numpy(float)); L.append('  leave out %-15s pooled E %+.3f  [%+.3f, %+.3f] -> %s' % (t, r['mu'], r['lo'], r['hi'], MA.verdict(r['lo'], r['hi'])))
    pd.DataFrame(res).to_csv(os.path.join(OUT, 'gene_graph_meta.csv'), index=False); open(os.path.join(OUT, 'gene_graph_report.md'), 'w', encoding='utf8').write('\n'.join(L)); print('\n'.join(L))
