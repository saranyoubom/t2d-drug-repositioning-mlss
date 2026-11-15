"""Tests for coexp_checks.py. (1) Regression on the INS-1 matrix: every number of benchmark_route3/validate_pipeline.py and of the manuscript (pipeline_17genes.json).
(2) The vectorized primitives against reference implementations (statsmodels BH, scipy Spearman, explicit loops, Monte Carlo geometry). (3) The permutation engine (within-group, block sharing,
chunk invariance, uniformity). (4) Synthetic known-answer tests: independent data give approximately uniform P; injected coupling is detected; a shared denominator inflates the
ratio-scheme adjusted counts while the calibrated P stays unremarkable. (5) Command line. (6) Speed.
Usage (runs about 3 minutes, exit code 1 if any check fails; env INS1_MATRIX overrides the matrix path):
  python -I test_coexp_checks.py > wp1_ins1_structure/test_output.txt"""
import sys, os, time, json, tempfile, contextlib, io
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'benchmark_route3'))     # -I removes the script folder from the path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import coexp_checks as cc, bench_pipeline as B

SRC = os.environ.get('INS1_MATRIX', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'R_17genes', 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx'))
REF = json.load(open(os.path.join(HERE, '..', 'reanalysis_17genes', 'pipeline_17genes.json')))
G17 = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']
FAILED, NCHECK = [], [0]


def section(t): print('\n== ' + t)


def check(name, cond, detail=''):
    NCHECK[0] += 1; cond = bool(cond); print(('PASS  ' if cond else 'FAIL  ') + name + ((' | ' + detail) if detail else ''))
    if not cond: FAILED.append(name)


def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n); h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)); return (c - h) / d, (c + h) / d


print('coexp_checks', cc.__version__, '| numpy', np.__version__, '| python', sys.version.split()[0])

# ------------------------------------------------------------------------------------------------------------------ 1. INS-1 regression
section('1. INS-1 regression: 17 genes, 3 groups x 4 wells, denominators Ins1, Kcnj11, Hes1')
X = pd.read_excel(SRC)[G17].astype(float).values; grp = np.repeat(['CTRL', 'DAPT', 'DKK-1'], 4); den = ['Ins1', 'Kcnj11', 'Hes1']; denx = [G17.index(g) for g in den]
des = cc.Design(grp); o = cc.observed_summary(X, grp, den, G17)
check('consensus pairs, any sign = 21', o['consensus_any'] == 21 == len(REF['consensus']), 'got %d' % o['consensus_any'])
check('consensus pairs, same sign in all four schemes = 21', o['consensus_same_sign'] == 21 == len(REF['consensus_signcons']), 'got %d' % o['consensus_same_sign'])
check('same 21 consensus pairs as pipeline_17genes.json', {frozenset(p) for p in o['consensus_pairs']} == {frozenset(p) for p in REF['consensus']})
check('pooled significant pairs raw/Ins1/Kcnj11/Hes1 = 66/80/84/120', o['pooled_sig_per_scheme'] == [66, 80, 84, 120] == [c['sig'] for c in REF['counts'].values()], str(o['pooled_sig_per_scheme']))
check('adjusted q < 0.05 raw/Ins1/Kcnj11/Hes1 = 0/51/11/15', o['adj_q_lt_alpha_per_scheme'] == [0, 51, 11, 15] == [int(v) for v in REF['adj_sig_q05'].values()], str(o['adj_q_lt_alpha_per_scheme']))
check('13 consensus pairs with q < 0.05 in >= 1 scheme', o['cons_adj_any'] == 13 == len(REF['cons_adj_any']), 'got %d' % o['cons_adj_any'])
check('0 consensus pairs with q < 0.05 in all four schemes; 0 pairs significant after adjustment in all four', o['cons_adj_all_schemes'] == 0 and o['adj_pairs_all_schemes'] == 0 == len(REF['adj_all_four']))
check('91 eligible consensus pairs (14 non-denominator genes)', o['n_eligible_pairs'] == 91)
check('hub degrees Glp1r 6; Tcf7l2, Cacna1d, Ptbp1 5; Wnt5a, Cacna1c 4', [o['degree'].get(g) for g in ('Glp1r', 'Tcf7l2', 'Cacna1d', 'Ptbp1', 'Wnt5a', 'Cacna1c')] == [6, 5, 5, 5, 4, 4])
S_old, SG_old, cnt_old = B.pooled(X, denx); Q_old, qcnt_old = B.adjusted(X, denx)
check('pooled significance matrices identical to bench_pipeline.pooled', all(np.array_equal(a, b) for a, b in zip(S_old, o['_S'])))
check('adjusted q matrices equal bench_pipeline.adjusted (statsmodels BH) to 1e-12', all(np.allclose(np.nan_to_num(a, nan=-1), np.nan_to_num(b, nan=-1), atol=1e-12, rtol=0) for a, b in zip(Q_old, o['_Q'])))
sh = cc.group_mean_share(X, grp, np.triu(o['_S'][0], 1)); sh_old = B.between_share(X, np.triu(S_old[0], 1))
check('group-mean covariance share, median over 66 raw significant pairs = 0.97', round(sh['median'], 2) == 0.97 and abs(sh['median'] - sh_old) < 1e-12, 'median %.4f mean %.4f pooled %.4f' % (sh['median'], sh['mean'], sh['pooled']))
cn = cc.consensus_null(X, grp, den, nperm=10000, seed=2026, genes=G17)
check('consensus null mean 22.1 (10000 permutations; manuscript 22.1)', abs(cn['null_mean'] - REF['perm']['mean']) <= 0.4, 'mean %.2f sd %.2f (manuscript %.2f, %.2f)' % (cn['null_mean'], cn['null_sd'], REF['perm']['mean'], REF['perm']['sd']))
check('consensus null SD within 0.3 of the manuscript', abs(cn['null_sd'] - REF['perm']['sd']) <= 0.3)
check('P(null >= 21) about 0.69', abs(cn['p_ge'] - REF['perm']['p_ge_obs']) <= 0.04, 'got %.3f (manuscript %.3f)' % (cn['p_ge'], REF['perm']['p_ge_obs']))
check('same-sign null mean within 0.4 of the manuscript', abs(cn['same_sign']['null_mean'] - REF['perm_signcons']['mean']) <= 0.4, 'got %.2f (manuscript %.2f)' % (cn['same_sign']['null_mean'], REF['perm_signcons']['mean']))
t0 = time.time(); ca_old, cs_old = B.permutation_counts(X, denx, 300, np.random.default_rng(11)); t_old = time.time() - t0; cn2 = cc.consensus_null(X, grp, den, nperm=3000, seed=7, genes=G17)
se = np.sqrt(ca_old.var() / 300 + cn2['_null'].var() / 3000)
check('vectorized null agrees with the loop null (300 loop permutations): means within 3 SE, SDs within 25%', abs(ca_old.mean() - cn2['null_mean']) <= 3 * se and abs(ca_old.std() / cn2['null_sd'] - 1) < 0.25,
      'loop %.2f (sd %.2f), vectorized %.2f (sd %.2f)' % (ca_old.mean(), ca_old.std(), cn2['null_mean'], cn2['null_sd']))
an = cc.adjusted_null(X, grp, den, nperm=2000, seed=2026, genes=G17); a = an['schemes']
check('calibrated adjusted counts: observed 0/51/11/15 from the threshold counter equals the exact BH counts', [a[n]['observed'] for n in an['scheme_names']] == [0, 51, 11, 15])
print('      adjusted q<0.05 [raw, Ins1, Kcnj11, Hes1]: observed %s | null mean %s | P(null >= obs) %s | excess %s' % ([a[n]['observed'] for n in an['scheme_names']], [round(a[n]['null_mean'], 1) for n in an['scheme_names']],
                                                                                                             [round(a[n]['p_ge'], 3) for n in an['scheme_names']], [round(a[n]['excess'], 1) for n in an['scheme_names']]))
check('adjusted nulls close to the planning reconnaissance (1000 permutations: 0.1, 29.0, 11.6, 6.4; P 1.0, 0.087, 0.465, 0.074)', abs(a['raw']['null_mean'] - 0.1) < 0.3 and abs(a['Ins1-normalized']['null_mean'] - 29.0) < 1.5 and
      abs(a['Kcnj11-normalized']['null_mean'] - 11.6) < 1.0 and abs(a['Hes1-normalized']['null_mean'] - 6.4) < 1.0 and abs(a['Hes1-normalized']['p_ge'] - 0.074) < 0.04 and a['raw']['p_ge'] == 1.0)
ct = cc.coupling_test(X, grp, nperm=10000, seed=2026)
print('      INS-1 coupling (10000 permutations): ' + ' | '.join('%s %.4g (null %.4g, P %.4f)' % (k, v['observed'], v['null_mean'], v['p_ge_plus1']) for k, v in ct.items()))
rc = cc.run_checks(X, grp, G17, den, nperm=2000)
check('run_checks warns: consensus count within null, ratio schemes not calibrated, group means carry > 90% of covariance',
      any('consensus count' in w for w in rc['warnings']) and any('Ins1-normalized' in w for w in rc['warnings']) and any('group means carry' in w for w in rc['warnings']), '%d warnings' % len(rc['warnings']))

# ------------------------------------------------------------------------------------------------------------------ 2. primitives
section('2. vectorized primitives against reference implementations')
rng = np.random.default_rng(1)
good = True
for _ in range(200): p = np.round(rng.random(60) ** 2, 3); good &= np.allclose(cc._bh(p), multipletests(p, method='fdr_bh')[1], atol=1e-12)
check('BH adjusted p equals statsmodels fdr_bh (ties and zeros included)', good)
tot, good = 0, True
for _ in range(150):
    Xs, g = cc.simulate([4, 4, 4], 17, rng); Xs = cc.inject_factor(Xs, g, list(range(10)), 0.8, rng); r = cc._corr(Xs.T @ des.Q); iu = np.triu_indices(17, 1)
    exact = int((cc._bh(cc._pval(r[iu], des.df_resid)) < 0.05).sum()); fast = int(cc._bh_count(np.abs(r[iu])[None], des.df_resid)[0]); good &= exact == fast; tot += exact
check('threshold BH counter equals exact BH count on 150 data sets with many significant pairs', good, 'total significant pairs %d' % tot)
Rb = rng.standard_normal((6, 17, 12)); Rb -= Rb.mean(2, keepdims=True); st = cc._coupling_stats(Rb, des); ok = True
for b in range(6):
    r = np.corrcoef(Rb[b]); iu = np.triu_indices(17, 1)
    ok &= abs(st['T_r2'][b] - (r[iu] ** 2).mean()) < 1e-12 and st['T_p'][b] == (np.abs(r[iu]) >= des.rc_resid).sum() and abs(st['T_eig'][b] - np.linalg.eigvalsh(r)[-1]) < 1e-9
Rb2 = rng.standard_normal((4, 6, 12)); st2 = cc._coupling_stats(Rb2, des)
ok &= all(abs(st2['T_eig'][b] - np.linalg.eigvalsh(np.corrcoef(Rb2[b]))[-1]) < 1e-9 for b in range(4))
check('T_r2, T_p, T_eig equal explicit np.corrcoef / eigvalsh computations (G > N via the Gram matrix, G < N directly)', ok)
rho = stats.spearmanr(X).statistic; osp = cc.observed_summary(X, grp, den, G17, kind='spearman'); n_sp = [int((np.abs(np.triu(rho, 1)) >= des.rc_pooled).sum())]
for d in denx:
    keep = [i for i in range(17) if i != d]; rr = stats.spearmanr(X[:, keep] / X[:, [d]]).statistic; n_sp.append(int((np.abs(np.triu(rr, 1)) >= des.rc_pooled).sum()))
check('Spearman pooled counts per scheme equal scipy.stats.spearmanr', osp['pooled_sig_per_scheme'] == n_sp, '%s' % n_sp)
olg = cc.observed_summary(X, grp, den, G17, kind='log'); n_lg = [int((np.abs(np.triu(np.corrcoef(np.log(X).T), 1)) >= des.rc_pooled).sum())]
for d in denx:
    keep = [i for i in range(17) if i != d]; n_lg.append(int((np.abs(np.triu(np.corrcoef(np.log(X[:, keep] / X[:, [d]]).T), 1)) >= des.rc_pooled).sum()))
check('log-scale Pearson pooled counts per scheme equal np.corrcoef of log values', olg['pooled_sig_per_scheme'] == n_lg, '%s' % n_lg)
check('critical values: r_crit(10) = 0.576, r_crit(8) = 0.632 for 3 x 4', abs(des.rc_pooled - 0.5760) < 5e-4 and abs(des.rc_resid - 0.6319) < 5e-4 and des.df_resid == 8)
d2 = cc.Design(['a'] * 3 + ['b'] * 5); d4 = cc.Design([0] * 4 + [1] * 5 + [2] * 6 + [3] * 3)
check('Design handles K = 2 unequal sizes (df_resid 5) and K = 4 (df_resid 13)', d2.df_resid == 5 and d4.df_resid == 13 and d2.Q.shape == (8, 8) and np.allclose(d4.Q.sum(0), 0))
sa = np.random.default_rng(3); a2 = sa.standard_normal((1_000_000, 2)); b2 = sa.standard_normal((1_000_000, 2)); c2 = np.abs((a2 * b2).sum(1) / np.linalg.norm(a2, axis=1) / np.linalg.norm(b2, axis=1))
a3 = sa.standard_normal((1_000_000, 3)); b3 = sa.standard_normal((1_000_000, 3)); c3 = np.abs((a3 * b3).sum(1) / np.linalg.norm(a3, axis=1) / np.linalg.norm(b3, axis=1)); e3, e4 = cc.geometric_expectation(0.576, 3), cc.geometric_expectation(0.576, 4)
check('geometric expectation K = 3: (2/pi) arccos(0.576) = 0.609 and Monte Carlo', abs(e3 - 2 / np.pi * np.arccos(0.576)) < 1e-12 and abs(e3 - (c2 >= 0.576).mean()) < 0.003, 'formula %.4f MC %.4f' % (e3, (c2 >= 0.576).mean()))
check('geometric expectation K = 4: 1 - r_crit = 0.424 and Monte Carlo; K = 2: 1', abs(e4 - 0.424) < 1e-12 and abs(e4 - (c3 >= 0.576).mean()) < 0.003 and cc.geometric_expectation(0.576, 2) == 1.0, 'MC %.4f' % (c3 >= 0.576).mean())

# ------------------------------------------------------------------------------------------------------------------ 3. permutation engine
section('3. permutation engine')
M = np.tile(np.arange(12.0)[:, None], (1, 6)); gr = ['a'] * 4 + ['b'] * 4 + ['c'] * 4; bl = ['x', 'x', 'y', 'y', 'y', None]
bat = np.concatenate(list(cc.permuted_batches(M, gr, 500, seed=5, blocks=bl)))                     # entries are the source row index of each permuted sample
within = all(set(bat[:, g, 4 * k:4 * k + 4].ravel()) <= set(range(4 * k, 4 * k + 4)) and (np.sort(bat[:, g, 4 * k:4 * k + 4], axis=1) == np.arange(4 * k, 4 * k + 4)).all() for g in range(6) for k in range(3))
check('every permutation stays inside its group (each gene a full within-group permutation)', within)
check('genes of one block share the permutation, different blocks differ', (bat[:, 0] == bat[:, 1]).all() and (bat[:, 2] == bat[:, 3]).all() and (bat[:, 3] == bat[:, 4]).all() and not (bat[:, 0] == bat[:, 2]).all() and not (bat[:, 4] == bat[:, 5]).all())
a1 = np.concatenate(list(cc.permuted_batches(M, gr, 700, seed=9, blocks=bl))); a2_ = np.concatenate(list(cc.permuted_batches(M, gr, 700, seed=9, blocks=bl, max_elems=150)))
check('result does not depend on the memory chunking (max_elems)', np.array_equal(a1, a2_))
one = np.concatenate(list(cc.permuted_batches(np.arange(8.0)[:, None], [0] * 4 + [1] * 4, 24000, seed=1))); codes = [tuple(int(v) for v in row[0, :4]) for row in one]
cnts = pd.Series(codes).value_counts(); chi = stats.chisquare(cnts.values) if len(cnts) == 24 else None
check('group of 4: all 24 permutations occur equally often (chi-square P > 0.001)', chi is not None and chi.pvalue > 0.001, 'P %.3f' % (chi.pvalue if chi else float('nan')))
gl = [0] * 10 + [1] * 10; M10 = np.tile(np.arange(20.0)[:, None], (1, 2)); b10 = np.concatenate(list(cc.permuted_batches(M10, gl, 20000, seed=4)))
pos = np.array([[(b10[:, 0, k] == s).mean() for s in range(10)] for k in range(10)])
check('group of 10 (random-key path): each source row lands in each position with probability 0.1', np.abs(pos - 0.1).max() < 0.01 and (np.sort(b10[:, 0, :10], axis=1) == np.arange(10)).all(), 'max deviation %.4f' % np.abs(pos - 0.1).max())
Xr = np.random.default_rng(8).lognormal(0, 1, (12, 6)); one_block = ['r'] * 6
cn1 = cc.consensus_null(Xr, gr, None, nperm=40, seed=3, blocks=one_block); ct1 = cc.coupling_test(Xr, gr, nperm=40, seed=3, blocks=one_block)
check('one block holding all genes leaves every correlation unchanged (null equals observed, P = 1)', cn1['null_sd'] == 0 and cn1['null_mean'] == cn1['observed'] and ct1['T_r2']['null_sd'] < 1e-12 and ct1['T_r2']['p_ge'] == 1.0)
pairs_ok = True
for Mb in cc.permuted_batches(Xr, gr, 50, seed=3, blocks=['p', 'p', 'p', 'q', 'q', 'q']):
    for b in range(50): pairs_ok &= abs(np.corrcoef(Mb[b, 0], Mb[b, 1])[0, 1] - np.corrcoef(Xr[:, 0], Xr[:, 1])[0, 1]) < 1e-12 and abs(np.corrcoef(Mb[b, 3], Mb[b, 5])[0, 1] - np.corrcoef(Xr[:, 3], Xr[:, 5])[0, 1]) < 1e-12
check('run-block permutation keeps same-block gene pairs exactly and breaks between-block pairs', pairs_ok and abs(np.corrcoef(Mb[0, 0], Mb[0, 3])[0, 1] - np.corrcoef(Xr[:, 0], Xr[:, 3])[0, 1]) > 1e-6)

# ------------------------------------------------------------------------------------------------------------------ 4. synthetic known-answer tests
section('4. synthetic known-answer tests (3 groups x 4, 17 genes, group effects SD 1.0, noise SD 0.3 on the log scale)')
rs = np.random.default_rng(2026); ND = 250; P = {k: [] for k in ('T_r2', 'T_p', 'T_eig', 'count')}; t0 = time.time()
for i in range(ND):
    Xs, g = cc.simulate([4, 4, 4], 17, rs); c = cc.coupling_test(Xs, g, nperm=199, seed=1000 + i); k = cc.consensus_null(Xs, g, [14, 15, 16], nperm=199, seed=1000 + i)
    for s in ('T_r2', 'T_p', 'T_eig'): P[s].append(c[s]['p_ge_plus1'])
    P['count'].append(k['p_ge_plus1'])
print('      %d independent data sets, 199 permutations each (%.0fs)' % (ND, time.time() - t0))
for s, lim in (('T_r2', 0.09), ('T_eig', 0.09), ('T_p', 0.09), ('count', 0.09)):
    p = np.array(P[s]); ks = stats.kstest(p, 'uniform'); err = (p <= 0.05).mean()
    lo, hi = wilson(int((p <= 0.05).sum()), ND)
    check('independent data, %s: P approximately uniform (mean 0.5 +- 0.07, type I error at 0.05 within 0.015-%.2f)' % (s, lim), abs(p.mean() - 0.5) < (0.07 if s in ('T_r2', 'T_eig') else 0.12) and 0.015 <= err <= lim,
          'mean P %.3f, KS D %.3f, P(P <= 0.05) = %.3f (Wilson %.3f-%.3f)' % (p.mean(), ks.statistic, err, lo, hi))
pu = np.array(P['T_r2']); ksr = stats.kstest(pu, 'uniform'); check('independent data, T_r2: Kolmogorov-Smirnov against uniform P > 0.01', ksr.pvalue > 0.01, 'KS P %.3f' % ksr.pvalue)
NR = 100; power = {}
for (r, m) in ((0.7, 6), (0.9, 10), (0.5, 3)):
    hit = {'T_r2': 0, 'T_p': 0, 'T_eig': 0, 'count': 0}
    for i in range(NR):
        Xs, g = cc.simulate([4, 4, 4], 17, rs); Xs = cc.inject_factor(cc.decouple(Xs, g, rs), g, list(range(m)), r, rs); c = cc.coupling_test(Xs, g, nperm=199, seed=5000 + i); k = cc.consensus_null(Xs, g, [14, 15, 16], nperm=199, seed=5000 + i)
        for s in ('T_r2', 'T_p', 'T_eig'): hit[s] += c[s]['p_ge_plus1'] < 0.05
        hit['count'] += k['p_ge_plus1'] < 0.05
    power[(r, m)] = hit; print('      injected factor r = %.1f in m = %d genes (decoupled data, %d replicates): detection T_r2 %.2f, T_p %.2f, T_eig %.2f, consensus count %.2f' % (r, m, NR, hit['T_r2'] / NR, hit['T_p'] / NR, hit['T_eig'] / NR, hit['count'] / NR))
check('injected coupling r = 0.7 in 6 genes is detected by T_r2 in >= 80% of data sets (protocol 3.2 requirement, here on synthetic data)', power[(0.7, 6)]['T_r2'] / NR >= 0.80, '%.2f' % (power[(0.7, 6)]['T_r2'] / NR))
check('injected coupling r = 0.9 in 10 genes is detected by T_r2 in >= 95% of data sets', power[(0.9, 10)]['T_r2'] / NR >= 0.95)
check('the consensus count is less sensitive than T_r2 to the same injected coupling (r = 0.7, m = 6), as in benchmark part 1', power[(0.7, 6)]['count'] < power[(0.7, 6)]['T_r2'], 'count %.2f vs T_r2 %.2f' % (power[(0.7, 6)]['count'] / NR, power[(0.7, 6)]['T_r2'] / NR))
Xl, gl_ = cc.simulate([30, 30, 30], 10, rs, noise_sd=0.3); Xi = cc.inject_factor(cc.decouple(Xl, gl_, rs), gl_, list(range(6)), 0.7, rs); rl = cc.residual_corr(Xi, gl_); iu6 = np.triu_indices(6, 1); rn = cc.residual_corr(cc.decouple(Xl, gl_, rs), gl_)
check('inject_factor: mean residual correlation among the 6 injected genes = 0.7 +- 0.08 (3 x 30 samples); decoupled genes about 0', abs(rl[:6, :6][iu6].mean() - 0.7) < 0.08 and abs(rn[np.triu_indices(10, 1)].mean()) < 0.03,
      'injected %.3f, decoupled %.3f' % (rl[:6, :6][iu6].mean(), rn[np.triu_indices(10, 1)].mean()))
dd = cc.decouple(Xl, gl_, rs); check('decouple keeps every group mean of the log values exactly', all(np.allclose(np.log(Xl)[r].mean(0), np.log(dd)[r].mean(0)) for r in cc.Design(gl_).rows))
NS = 60; raw_c, ratio_c, nullm, pr_ = [], [], [], []
for i in range(NS):
    Xs, g = cc.simulate([4, 4, 4], 17, rs, noise_sd=0.3, ref_noise_sd=0.6, ref_gene=[14, 15, 16]); a = cc.adjusted_null(Xs, g, [14, 15, 16], nperm=199, seed=7000 + i)
    s = a['schemes']; raw_c.append(s['raw']['observed']); ratio_c.append(np.mean([s[n]['observed'] for n in a['scheme_names'][1:]])); nullm.append(np.mean([s[n]['null_mean'] for n in a['scheme_names'][1:]])); pr_.append([s[n]['p_ge_plus1'] for n in a['scheme_names'][1:]])
pr_ = np.array(pr_); print('      noisy shared denominators (extra log SD 0.6): adjusted q<0.05 count mean raw %.2f, ratio schemes %.1f (of 120 pairs); null mean of the ratio schemes %.1f; share of ratio-scheme P <= 0.05: %.3f' % (np.mean(raw_c), np.mean(ratio_c), np.mean(nullm), (pr_ <= 0.05).mean()))
check('independent genes with noisy shared denominators: ratio-scheme adjusted counts are inflated (>= 10x the raw count, > 20 pairs)', np.mean(ratio_c) > 20 and np.mean(ratio_c) >= 10 * max(np.mean(raw_c), 0.1))
check('the calibrated null reproduces that inflation: P(null >= obs) of the ratio schemes is not small (share <= 0.05 is at most 0.10)', (pr_ <= 0.05).mean() <= 0.10 and abs(np.mean(ratio_c) - np.mean(nullm)) < 0.15 * np.mean(ratio_c))
# ------------------------------------------------------------------------------------------------------------------ run structure
Dp = []; Ds = []
for i in range(60):
    Xs, g = cc.simulate([4, 4, 4], 12, rs); runs = ['A'] * 5 + ['B'] * 4 + ['C'] * 3; Dp.append(cc.run_structure_test(Xs, g, runs, nperm=999, seed=i)['r']['p_two_sided'])
    Xc = cc.decouple(Xs, g, rs)
    for lab, idx in (('A', range(0, 5)), ('B', range(5, 9))): Xc = cc.inject_factor(Xc, g, list(idx), 0.6, rs)
    Ds.append(cc.run_structure_test(Xc, g, runs, nperm=999, seed=i)['r']['p_two_sided'])
check('run-structure test: P approximately uniform on independent genes (mean 0.5 +- 0.1, P <= 0.05 in at most 12%)', abs(np.mean(Dp) - 0.5) < 0.1 and (np.array(Dp) <= 0.05).mean() <= 0.12, 'mean %.3f, share %.3f' % (np.mean(Dp), (np.array(Dp) <= 0.05).mean()))
check('run-structure test: detects within-run coupling (P <= 0.05 in >= 70% of data sets with r = 0.6 inside runs A and B)', (np.array(Ds) <= 0.05).mean() >= 0.70, 'share %.2f' % (np.array(Ds) <= 0.05).mean())

# ------------------------------------------------------------------------------------------------------------------ 5. other designs and the command line
section('5. other designs and the command line')
Xo, go = cc.simulate([3, 5], 8, rs); ro = cc.run_checks(Xo, ['a'] * 3 + ['b'] * 5, nperm=99, run_nperm=99)
Xk, gk = cc.simulate([4, 5, 6, 3], 9, rs); rk = cc.run_checks(Xk, gk, den=[0, 1], runs=list('aabbbccdd'), nperm=99, run_nperm=99)
check('run_checks runs for K = 2 (unequal sizes) and K = 4 with run labels; geometry expectation 1.0 and 0.42 at the pooled r_crit', ro['meta']['df_resid'] == 5 and rk['meta']['df_resid'] == 13 and ro['geometry']['expected'] == 1.0 and abs(rk['geometry']['expected'] - (1 - rk['geometry']['r_crit'])) < 1e-9 and 'run_structure' in rk)
with tempfile.TemporaryDirectory() as td:
    Xs, g = cc.simulate([4, 4, 4], 10, rs); df = pd.DataFrame(Xs, columns=['g%d' % i for i in range(10)]); df.insert(0, 'sample', ['s%d' % i for i in range(12)]); csv = os.path.join(td, 'd.csv'); df.to_csv(csv, index=False)
    pd.DataFrame({'gene': ['g%d' % i for i in range(10)], 'run': ['A'] * 5 + ['B'] * 4 + ['']}).to_csv(os.path.join(td, 'runs.csv'), index=False)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf): rcode = cc.main([csv, '--groups', '4,4,4', '--denominators', 'g7,g8,g9', '--runs', os.path.join(td, 'runs.csv'), '--nperm', '99', '--run-nperm', '99', '--out', os.path.join(td, 'o.json'), '--report', os.path.join(td, 'o.md')])
    js = json.load(open(os.path.join(td, 'o.json'))); md = open(os.path.join(td, 'o.md'), encoding='utf8').read()
    check('command line: exit code 0, json with warnings/coupling/run_structure, markdown report, unresolved gene excluded from the run test', rcode == 0 and 'warnings' in js and 'T_r2' in js['coupling'] and js['run_structure']['excluded'] == ['g9'] and '## Warnings' in md and 'Within-group coupling' in md)
    with contextlib.redirect_stdout(io.StringIO()): cc.main([csv, '--groups', 'A,A,A,A,B,B,B,B,C,C,C,C', '--nperm', '99', '--out', os.path.join(td, 'o2.json')])
    check('command line: one label per sample accepted; no denominators gives the raw scheme only', json.load(open(os.path.join(td, 'o2.json')))['meta']['denominators'] == [])
msg = ''
try: cc.run_checks(np.array([[1.0, 0.0], [1, 2], [3, 4], [5, 6]]), [0, 0, 1, 1])
except ValueError as e: msg = str(e)
check('non-positive input is refused with a clear message', 'positive' in msg, msg[:60])

# ------------------------------------------------------------------------------------------------------------------ 6. speed
section('6. speed (this machine, G = 17, N = 12, 3 x 4, single process)')
Xs, g = cc.simulate([4, 4, 4], 17, np.random.default_rng(5)); cc.coupling_test(Xs, g, nperm=300)


def rate(f, n, rep=3):
    best = 1e9
    for _ in range(rep): t = time.perf_counter(); f(); best = min(best, time.perf_counter() - t)
    return n / best


sp = dict(T_r2_T_p=rate(lambda: cc.coupling_test(Xs, g, nperm=4000, eig=False), 4000), T_r2_T_p_T_eig=rate(lambda: cc.coupling_test(Xs, g, nperm=4000), 4000),
          consensus=rate(lambda: cc.consensus_null(Xs, g, [14, 15, 16], nperm=4000), 4000), adjusted_4_schemes=rate(lambda: cc.adjusted_null(Xs, g, [14, 15, 16], nperm=4000), 4000),
          permutation_generation=rate(lambda: [b for b in cc.permuted_batches(Xs, g, 20000)], 20000), spearman_consensus=rate(lambda: cc.consensus_null(Xs, g, [14, 15, 16], nperm=1000, kind='spearman'), 1000, 1),
          run_label_permutations=rate(lambda: cc.run_structure_test(Xs, g, list('aaaabbbbbcccccdde'), nperm=20000), 20000))
for k, v in sp.items(): print('      %-26s %8.0f permutations/s' % (k, v))
print('      loop null of bench_pipeline.permutation_counts (consensus): %.0f permutations/s -> vectorized consensus null is %.0fx faster' % (300 / t_old, sp['consensus'] / (300 / t_old)))
check('T_r2 null runs at several thousand permutations per second (>= 3000/s without T_eig)', sp['T_r2_T_p'] >= 3000, '%.0f/s' % sp['T_r2_T_p'])
check('vectorized consensus counting is at least 20x faster than the loop', sp['consensus'] >= 20 * (300 / t_old), '%.0fx' % (sp['consensus'] / (300 / t_old)))

print('\n%d checks, %d failed' % (NCHECK[0], len(FAILED)))
if FAILED:
    print('FAILED: ' + '; '.join(FAILED)); print('TESTS FAILED'); sys.exit(1)
print('ALL TESTS PASSED'); sys.exit(0)
