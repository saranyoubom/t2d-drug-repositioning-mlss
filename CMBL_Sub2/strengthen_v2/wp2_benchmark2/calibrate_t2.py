"""calibrate_t2.py - protocol v2 section 3.2: calibration of the within-group coupling tests on decoupled development data (WP2 deliverable B).
For each of the 37 development series, gene set A (the 12 x 17 matrices saved by gate_part1.py extract):
  decouple = permute every gene's log values independently within each group (group means kept, coupling removed), then
  (i)  no injection, 100 replicates per series: type I error of every test at P < 0.05;
  (ii) inject ONE latent factor, centred within group, into m genes (drawn at random among the 14 non-denominator genes, as in the part 1 positive control), loading sqrt(r/(1 - r)) x the gene's pooled
       within-group SD, so that two injected genes have within-group correlation r: r in 0.3, 0.5, 0.7, 0.9; m in 3, 6, 10; 20 replicates per series and condition.
Tests (coexp_checks.coupling_test, 500 within-group permutations each, P = (1 + #{null >= obs})/(500 + 1)): T_r2 (primary), T_p, T_eig. The consensus-COUNT test of part 1 (consensus_null, 500 permutations)
is run at the REDUCED replicate number declared here: m = 6 (all four r) and no injection, the first 20 replicates of each series x condition (the protocol requires at least 5).
Requirement fixed in the protocol: T_r2 type I error <= 0.06 and power >= 80% at r = 0.7, m = 6 (pooled over all series x replicates; the spread over series is reported as well).
Seeds: numpy SeedSequence([2026, GSE number, condition, replicate]). Usage: python -I calibrate_t2.py run [workers] | summarize"""
import sys, os, json, time, itertools
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
COMPACT = os.path.join(HERE, 'compact'); SERDIR = os.path.join(COMPACT, 'calib_series'); os.makedirs(SERDIR, exist_ok=True)
RS = [0.3, 0.5, 0.7, 0.9]; MS = [3, 6, 10]; NREP_INJ = 20; NREP_NULL = 100; NREP_COUNT = 20; NPERM = 500
CONDITIONS = [(0.0, 0, NREP_NULL)] + [(r, m, NREP_INJ) for r in RS for m in MS]


def one_series(gse):
    import bench2_core as C                                  # first: it puts strengthen_v2/ on sys.path, where coexp_checks lives
    import coexp_checks as cc
    X0 = np.load(os.path.join(COMPACT, 'dev_sets.npz'))['%s__A' % gse]; n = int(gse[3:]); rows = []; t0 = time.time()
    for ci, (r, m, nrep) in enumerate(CONDITIONS):
        for rep in range(nrep):
            rng = np.random.default_rng(np.random.SeedSequence([2026, n, ci, rep])); Xd = cc.decouple(X0, C.DES, rng)
            if m > 0: idx = sorted(rng.choice(np.arange(3, 17), m, replace=False).tolist()); Xi = cc.inject_factor(Xd, C.DES, idx, r, rng)
            else: idx = []; Xi = Xd
            sd = int(rng.integers(2 ** 31 - 1)); ct = cc.coupling_test(Xi, C.DES, nperm=NPERM, seed=sd, eig=True)
            row = dict(gse=gse, r=r, m=m, rep=rep, injected=' '.join(map(str, idx)), T_r2_obs=ct['T_r2']['observed'], T_r2_p=ct['T_r2']['p_ge_plus1'], T_p_obs=ct['T_p']['observed'], T_p_p=ct['T_p']['p_ge_plus1'],
                       T_eig_obs=ct['T_eig']['observed'], T_eig_p=ct['T_eig']['p_ge_plus1'], count_obs=None, count_p=None, count_p_plus1=None)
            if (m in (0, 6)) and rep < NREP_COUNT:
                cn = cc.consensus_null(Xi, C.DES, C.DEN, nperm=NPERM, seed=sd); row.update(count_obs=cn['observed'], count_p=cn['p_ge'], count_p_plus1=cn['p_ge_plus1'])
            rows.append(row)
    json.dump(rows, open(os.path.join(SERDIR, gse + '.json'), 'w')); return gse, len(rows), round(time.time() - t0, 1)


def run(nw):
    import multiprocessing as mp
    prep = json.load(open(os.path.join(COMPACT, 'dev_prepared.json'))); todo = [g for g in prep if not os.path.exists(os.path.join(SERDIR, g + '.json'))]
    print('series to run: %d of %d, workers %d' % (len(todo), len(prep), nw), flush=True)
    with mp.get_context('spawn').Pool(nw) as pool:
        for gse, nrow, sec in pool.imap_unordered(one_series, todo): print('done', gse, nrow, 'rows', sec, 's', flush=True)
    print('calibration run finished', flush=True)


def wilson(k, n):
    from statsmodels.stats.proportion import proportion_confint
    lo, hi = proportion_confint(int(k), int(n), method='wilson'); return float(lo), float(hi)


def summarize():
    import pandas as pd
    prep = json.load(open(os.path.join(COMPACT, 'dev_prepared.json'))); frames = []
    for g in prep:
        f = os.path.join(SERDIR, g + '.json')
        if not os.path.exists(f): raise SystemExit('missing series file ' + f)
        frames.append(pd.DataFrame(json.load(open(f))))
    R = pd.concat(frames, ignore_index=True); R['platform'] = R.gse.map(lambda g: prep[g]['platform']); R.to_csv(os.path.join(HERE, 'calibration_replicates.csv'), index=False)
    stats_ = [('T_r2', 'T_r2_p'), ('T_p', 'T_p_p'), ('T_eig', 'T_eig_p'), ('count', 'count_p_plus1')]; per = []; tab = []
    for (r, m), g in R.groupby(['r', 'm']):
        for sname, col in stats_:
            gg = g[g[col].notna()] if sname == 'count' else g
            if len(gg) == 0: continue
            for gse, h in gg.groupby('gse'): per.append(dict(gse=gse, platform=prep[gse]['platform'], r=r, m=m, statistic=sname, n_rep=len(h), n_reject=int((h[col] < 0.05).sum())))
            k = int((gg[col] < 0.05).sum()); n = len(gg); lo, hi = wilson(k, n); ser = gg.groupby('gse')[col].apply(lambda x: (x < 0.05).mean())
            tab.append(dict(r=r, m=m, condition='no injection (type I error)' if m == 0 else 'injected', statistic=sname, n_series=gg.gse.nunique(), n_tests=n, n_reject=k, rate=k / n, ci_lo=lo, ci_hi=hi,
                            series_min=float(ser.min()), series_q1=float(ser.quantile(.25)), series_median=float(ser.median()), series_q3=float(ser.quantile(.75)), series_max=float(ser.max()),
                            n_series_rate_ge_0_80=int((ser >= 0.8).sum()), n_series_rate_gt_0_06=int((ser > 0.06).sum())))
    pd.DataFrame(per).to_csv(os.path.join(HERE, 'calibration_per_condition.csv'), index=False); T = pd.DataFrame(tab).sort_values(['m', 'r', 'statistic']); T.to_csv(os.path.join(HERE, 'calibration_tables.csv'), index=False)
    by_pl = []
    for pl, g in R.groupby('platform'):
        for (r, m), h in g.groupby(['r', 'm']):
            if (r, m) in ((0.0, 0), (0.7, 6)): by_pl.append(dict(platform=pl, r=r, m=m, n_tests=len(h), T_r2_rate=float((h.T_r2_p < 0.05).mean()), T_p_rate=float((h.T_p_p < 0.05).mean()), T_eig_rate=float((h.T_eig_p < 0.05).mean())))
    pd.DataFrame(by_pl).to_csv(os.path.join(HERE, 'calibration_by_platform.csv'), index=False)
    # ---- verdict (rules fixed in the protocol and in the docstring above)
    t1 = T[(T.m == 0) & (T.statistic == 'T_r2')].iloc[0]; pw = T[(T.m == 6) & (T.r == 0.7) & (T.statistic == 'T_r2')].iloc[0]
    ok1 = t1.rate <= 0.06; ok2 = pw.rate >= 0.80; V = ['PROTOCOL 3.2 REQUIREMENT (T_r2 type I error <= 0.06 and power >= 80%% at r = 0.7, m = 6; pooled over %d series)' % t1.n_series, '',
         'T_r2 type I error: %d/%d = %.4f (Wilson 95%% CI %.4f to %.4f) -> %s' % (t1.n_reject, t1.n_tests, t1.rate, t1.ci_lo, t1.ci_hi, 'PASS' if ok1 else 'FAIL'),
         'T_r2 power at r = 0.7, m = 6: %d/%d = %.4f (Wilson 95%% CI %.4f to %.4f) -> %s; per series: min %.2f, median %.2f, %d of %d series at or above 0.80' % (pw.n_reject, pw.n_tests, pw.rate, pw.ci_lo, pw.ci_hi, 'PASS' if ok2 else 'FAIL', pw.series_min, pw.series_median, pw.n_series_rate_ge_0_80, pw.n_series),
         'OVERALL: %s' % ('PASS' if (ok1 and ok2) else 'FAIL'), '']
    cand = []
    for s in ('T_r2', 'T_p', 'T_eig', 'count'):
        a = T[(T.m == 0) & (T.statistic == s)]; b = T[(T.m == 6) & (T.r == 0.7) & (T.statistic == s)]
        if len(a) and len(b): cand.append((s, float(a.iloc[0].rate), float(b.iloc[0].rate)))
    V.append('All pre-listed statistics (type I error, power at r = 0.7, m = 6): ' + '; '.join('%s %.3f, %.3f' % c for c in cand))
    good = [c for c in cand if c[1] <= 0.06 and c[2] >= 0.80]
    if not (ok1 and ok2):
        best = max([c for c in cand if c[1] <= 0.06], key=lambda c: c[2], default=None)
        V.append('Secondary statistic named by the fixed rule (highest power at r = 0.7, m = 6 among statistics with type I error <= 0.06): ' + (best[0] if best else 'none qualifies'))
    else: V.append('Pre-listed statistics meeting both requirements: ' + (', '.join(c[0] for c in good) or 'none'))
    # ---- post hoc sensitivity (protocol section 8, WP2 entry of 5 Nov): without the series whose gene set A has genes with no within-group variation
    deg = sorted({j['gse'] for j in (json.loads(l) for l in open(os.path.join(COMPACT, 'dev_core_results.jsonl'), encoding='utf8')) if j['tag'] == 'A' and j['flat'].get('n_zero_var_genes', 0) > 0})
    Sx = R[~R.gse.isin(deg)]; sens = []
    V += ['', 'POST HOC SENSITIVITY (not the pre-specified analysis): without the %d series whose gene set A has genes with no within-group variation (%s); %d series remain' % (len(deg), ', '.join(deg), Sx.gse.nunique())]
    for sname, col in stats_:
        a = Sx[Sx.m == 0]; b = Sx[(Sx.m == 6) & (Sx.r == 0.7)]
        if sname == 'count': a = a[a[col].notna()]; b = b[b[col].notna()]
        ka, na, kb, nb = int((a[col] < 0.05).sum()), len(a), int((b[col] < 0.05).sum()), len(b); la, ha = wilson(ka, na); lb, hb = wilson(kb, nb)
        sens.append(dict(statistic=sname, type1_k=ka, type1_n=na, type1=ka / na, type1_lo=la, type1_hi=ha, power_k=kb, power_n=nb, power=kb / nb, power_lo=lb, power_hi=hb))
        V.append('  %-6s type I error %d/%d = %.4f (%.4f to %.4f); power at r = 0.7, m = 6: %d/%d = %.4f (%.4f to %.4f)' % (sname, ka, na, ka / na, la, ha, kb, nb, kb / nb, lb, hb))
    pd.DataFrame(sens).to_csv(os.path.join(HERE, 'calibration_sensitivity_no_degenerate.csv'), index=False)
    open(os.path.join(HERE, 'calibration_verdict.txt'), 'w', encoding='utf8').write('\n'.join(V) + '\n')
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
    print('\n'.join(V)); print(T[['r', 'm', 'statistic', 'n_tests', 'n_reject', 'rate', 'ci_lo', 'ci_hi', 'series_min', 'series_median', 'n_series_rate_ge_0_80']].to_string(index=False, float_format=lambda x: '%.3f' % x))


if __name__ == '__main__':
    step = sys.argv[1] if len(sys.argv) > 1 else 'summarize'
    if step == 'run': run(int(sys.argv[2]) if len(sys.argv) > 2 else 2)
    elif step == 'summarize': summarize()
    else: raise SystemExit('run [workers] | summarize')
