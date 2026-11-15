"""Planning reconnaissance only (not a reported analysis): (1) adjusted counts per scheme under the permutation null (shared-denominator artefact); (2) same-run vs between-run pairs."""
import sys, os, itertools
import numpy as np, pandas as pd
from scipy import stats
sys.stdout.reconfigure(encoding='utf8')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'benchmark_route3'))
import bench_pipeline as B
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'R_17genes', 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx')
G = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']
X = pd.read_excel(SRC)[G].astype(float).values; den = [G.index(g) for g in ('Ins1', 'Kcnj11', 'Hes1')]
# (1) adjusted q<0.05 per scheme: observed vs null with independent within-group residuals
_, obs = B.adjusted(X, den); rng = np.random.default_rng(7); L = np.log(X); gm = np.vstack([L[B.GROUP == g].mean(0) for g in range(3)]); res = L - gm[B.GROUP]; idx = [np.where(B.GROUP == g)[0] for g in range(3)]
null = []
for b in range(1000):
    R = res.copy()
    for ix in idx: R[ix] = rng.permuted(R[ix], axis=0)
    _, c = B.adjusted(np.exp(gm[B.GROUP] + R), den); null.append(c)
null = np.array(null)
print('adjusted q<0.05 [raw, Ins1, Kcnj11, Hes1] observed:', obs, '| null mean:', np.round(null.mean(0), 1).tolist(), '| null 95th pct:', np.percentile(null, 95, axis=0).tolist(), '| P(null>=obs):', [round(float((null[:, k] >= obs[k]).mean()), 3) for k in range(4)])
# (2) run structure (two versions of the June genes)
RUN = {'Glp1r': 'Jun03', 'Glut2': 'Jun03', 'Hes1': 'Jun03', 'Ins1': 'Jun03', 'Lef1': 'Jun24', 'Tcf7': 'Jun24', 'Wnt5a': 'Jun24', 'Wnt5b': 'Jun24', 'Wnt9a': 'Jun24', 'Cacna1c': 'Aug26', 'Cacna1d': 'Aug26', 'Kcnj11': 'Aug26', 'Ptbp1': 'Aug26', 'Tcf7l2': 'Aug26', 'Itpr1': 'Nov09', 'Ryr2': 'Nov09', 'Ryr3': 'Nov09'}
S, SG, cnt = B.pooled(X, den); cons, _ = B.consensus_from(S, SG)
nd = [i for i in range(17) if i not in den]; rows = []
Rres = X - B.DESIGN @ np.linalg.lstsq(B.DESIGN, X, rcond=None)[0]; radj = np.corrcoef(Rres, rowvar=False); rraw = np.corrcoef(X, rowvar=False)
LR = L - gm[B.GROUP]; rlog = np.corrcoef(LR, rowvar=False)
for i, j in itertools.combinations(range(17), 2):
    rows.append(dict(a=G[i], b=G[j], same=RUN[G[i]] == RUN[G[j]], runs=RUN[G[i]] + '/' + RUN[G[j]], eligible=(i in nd and j in nd), cons=bool(cons[min(i, j), max(i, j)]), r_raw=rraw[i, j], r_adj=radj[i, j], r_adj_log=rlog[i, j], sig_raw=bool(np.triu(S[0], 1)[min(i, j), max(i, j)])))
T = pd.DataFrame(rows); E = T[T.eligible]
a, b, c, d = int((E.same & E.cons).sum()), int((E.same & ~E.cons).sum()), int((~E.same & E.cons).sum()), int((~E.same & ~E.cons).sum())
print(f'consensus by run, eligible pairs {len(E)}: same-run {a}/{a + b}, between-run {c}/{c + d}; Fisher OR {stats.fisher_exact([[a, b], [c, d]])[0]:.2f} p {stats.fisher_exact([[a, b], [c, d]])[1]:.3f}')
print('within-group (adjusted, raw scale) |r|: same-run median %.2f (n=%d), between-run median %.2f (n=%d); MWU p %.2g' % (T[T.same].r_adj.abs().median(), T.same.sum(), T[~T.same].r_adj.abs().median(), (~T.same).sum(), stats.mannwhitneyu(T[T.same].r_adj.abs(), T[~T.same].r_adj.abs()).pvalue))
print('within-group signed r: same-run mean %.2f, between-run mean %.2f' % (T[T.same].r_adj.mean(), T[~T.same].r_adj.mean()))
print('per run within-group mean r / mean |r| (same-run pairs):'); print(T[T.same].groupby('runs').agg(n=('r_adj', 'size'), mean_r=('r_adj', 'mean'), mean_abs=('r_adj', lambda s: s.abs().mean())).round(2).to_string())
print('pooled raw significant: same-run %d/%d; between-run %d/%d' % (T[T.same].sig_raw.sum(), T.same.sum(), T[~T.same].sig_raw.sum(), (~T.same).sum()))
print('same-run pairs |r_adj| >= 0.632 (p<0.05, df 8):', T[T.same & (T.r_adj.abs() >= 0.632)][['a', 'b', 'r_adj']].round(2).values.tolist())
print('between-run pairs |r_adj| >= 0.632:', int((T[~T.same].r_adj.abs() >= 0.632).sum()), 'of', int((~T.same).sum()))
