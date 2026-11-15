"""design_table.py - protocol v2 section 3.7 (WP2 deliverable F): design guidance by simulation.
G = 17 genes, 3 groups, n per group in 4, 6, 8, 12, 20, 30; one module of m = 6 genes with within-group correlation r in 0.3, 0.5, 0.7, 0.9 (coexp_checks.inject_factor), the other 11 genes independent;
r = 0 is added as a type I error check (ADD-ON, same code). Data: independent log-normal values from coexp_checks.simulate (group effects N(0, 1), residual SD 0.3 on the log scale), seed 2026.
Per cell R independent experiments (replicates here are independent experiments: wells from one preparation do not count). Outputs per cell:
  power_T_r2       share of experiments with P < 0.05 for T_r2 (coupling_test, 500 within-group permutations, P = (1 + k)/(500 + 1))
  pair_power       mean over the 15 true module pairs of the share detected at BH q < 0.05 (partial correlation of log values after removing group means, residual df = N - 4, 136 pairs, BH over all pairs)
  any_pair         share of experiments in which at least one true pair is detected
  false_pairs      mean number of detected pairs outside the module (expected near 0)
  r_best_pair      analytic |r| the single best of the 136 pairs needs for BH q < 0.05 (p <= 0.05/136) at the residual df; r_single is |r| for an uncorrected p < 0.05
Usage: python -I design_table.py run [workers] [R=300] | summarize"""
import sys, os, json, time
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))      # strengthen_v2/ holds coexp_checks.py
import numpy as np
COMPACT = os.path.join(HERE, 'compact'); DIR = os.path.join(COMPACT, 'design_cells'); os.makedirs(DIR, exist_ok=True)
NS = [4, 6, 8, 12, 20, 30]; RS = [0.0, 0.3, 0.5, 0.7, 0.9]; G = 17; M = 6; NPERM = 500


def cell(args):
    n, R = args
    import coexp_checks as cc
    rows = []; t0 = time.time(); iu = np.triu_indices(G, 1); true = (iu[0] < M) & (iu[1] < M)
    for ri, r in enumerate(RS):
        for rep in range(R):
            rng = np.random.default_rng(np.random.SeedSequence([2026, 37, n, ri, rep])); X, codes = cc.simulate([n, n, n], G, rng, effect_sd=1.0, noise_sd=0.3)
            if r > 0: X = cc.inject_factor(X, codes, range(M), r, rng)
            des = cc.Design(codes); ct = cc.coupling_test(X, des, nperm=NPERM, seed=int(rng.integers(2 ** 31 - 1)), eig=False)
            Rm = cc.residual_corr(X, des); pv = cc._pval(Rm[iu], des.df_resid); q = cc._bh(pv); det = q < 0.05
            rows.append(dict(n=n, r=r, rep=rep, T_r2_p=ct['T_r2']['p_ge_plus1'], pair_power=float(det[true].mean()), any_pair=bool(det[true].any()), false_pairs=int(det[~true].sum()), n_detected=int(det.sum())))
    json.dump(rows, open(os.path.join(DIR, 'n%d.json' % n), 'w')); return n, len(rows), round(time.time() - t0, 1)


def run(nw, R):
    import multiprocessing as mp
    todo = [(n, R) for n in sorted(NS, reverse=True) if not os.path.exists(os.path.join(DIR, 'n%d.json' % n))]; print('cells to run (n per group):', [t[0] for t in todo], 'R =', R, 'workers', nw, flush=True)
    with mp.get_context('spawn').Pool(nw) as pool:
        for n, k, sec in pool.imap_unordered(cell, todo): print('done n =', n, k, 'experiments', sec, 's', flush=True)
    print('design simulation finished', flush=True)


def summarize():
    import pandas as pd, coexp_checks as cc
    from statsmodels.stats.proportion import proportion_confint
    R = pd.concat([pd.DataFrame(json.load(open(os.path.join(DIR, 'n%d.json' % n)))) for n in NS], ignore_index=True); out = []
    for (n, r), g in R.groupby(['n', 'r']):
        df = 3 * n - 4; k = int((g.T_r2_p < 0.05).sum()); lo, hi = proportion_confint(k, len(g), method='wilson'); k2 = int(g.any_pair.sum()); lo2, hi2 = proportion_confint(k2, len(g), method='wilson')
        out.append(dict(n_per_group=n, N=3 * n, resid_df=df, r=r, replicates=len(g), power_T_r2=k / len(g), T_r2_lo=lo, T_r2_hi=hi, pair_power=float(g.pair_power.mean()), any_pair=k2 / len(g), any_pair_lo=lo2, any_pair_hi=hi2,
                        false_pairs=float(g.false_pairs.mean()), r_best_pair=cc.rcrit(df, 0.05 / 136), r_single=cc.rcrit(df, 0.05)))
    T = pd.DataFrame(out); T.to_csv(os.path.join(HERE, 'design_table.csv'), index=False)
    md = ['# Design table (protocol 3.7): 17 genes, one module of 6 genes, 3 groups, simulated independent experiments', '',
          'Power of T_r2 (P < 0.05, 500 permutations) / share of true module pairs detected at BH q < 0.05 / experiments with at least one true pair detected. r = 0 rows give the type I error of T_r2 and the mean number of false pairs.', '',
          '| n per group (N, residual df) | analytic abs(r) for the best pair | ' + ' | '.join('r = %.1f' % r for r in RS) + ' |', '|---|---|' + '---|' * len(RS)]
    for n in NS:
        g = T[T.n_per_group == n].sort_values('r'); md.append('| %d (%d, %d) | %.2f | ' % (n, 3 * n, 3 * n - 4, g.r_best_pair.iloc[0]) + ' | '.join('%.2f / %.2f / %.2f' % (x.power_T_r2, x.pair_power, x.any_pair) for x in g.itertuples()) + ' |')
    md += ['', 'Replicates per cell: %d independent simulated experiments. Wells or replicates from one preparation are not independent experiments and do not count toward n.' % int(T.replicates.iloc[0]), '']
    open(os.path.join(HERE, 'design_table.md'), 'w', encoding='utf8').write('\n'.join(md)); print('\n'.join(md))


if __name__ == '__main__':
    step = sys.argv[1] if len(sys.argv) > 1 else 'summarize'
    if step == 'run': run(int(sys.argv[2]) if len(sys.argv) > 2 else 2, int(sys.argv[3]) if len(sys.argv) > 3 else 300)
    elif step == 'summarize': summarize()
    else: raise SystemExit('run [workers] [R] | summarize')
