"""post_hoc_checks_191126.py - main session, 19 Nov 2026. Two POST HOC checks of the benchmark (method in POSTHOC_METHODS_191126.md, fixed before any result):
  check 1  share of the pooled covariance carried by group means, observed against its value when within-group coupling is removed (independent within-group permutation, selection of
           pooled-significant pairs repeated in every permuted matrix);
  check 2  the within-group coupling test T_r2 after removal of the first principal component of the residuals (the same removal in every permuted matrix; same permutations as the original test).
Read-only on the project: imports coexp_checks.py and bench2_core.py from strengthen_v2 (never modified), reads the stored matrices, writes only into this folder (results_<set>/).
Gates (before any summary is written): the script must reproduce the stored share_between and T_r2_p of dev_core_results.jsonl.
Usage (python -I): python -I post_hoc_checks_191126.py dev            development series (37), stored matrices, no download
                   python -I post_hoc_checks_191126.py validation DATA_DIR   validation sets rebuilt from GEO files in DATA_DIR (only after the user's permission to download)"""
import os, sys, json, csv, hashlib, time
import numpy as np
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.abspath(os.path.join(HERE, '..', '..'))
W2 = os.path.join(V2, 'wp2_benchmark2')
sys.path.insert(0, V2); sys.path.insert(0, W2)
import coexp_checks as cc
import bench2_core as B
DES = B.DES
NPERM = 2000
METHODS = os.path.join(HERE, 'POSTHOC_METHODS_191126.md')
sha256 = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()


def share_expected(X, nperm, seed):
    """observed median share and the distribution of the same statistic without within-group coupling"""
    X = np.asarray(X, float); N, G = X.shape
    obs = cc.group_mean_share(X, DES)
    gm = np.empty_like(X)
    for rows in DES.rows: gm[rows] = X[rows].mean(axis=0)
    bc = gm - X.mean(axis=0); btw = (bc.T @ bc) / N                                  # between-group covariance, unchanged by within-group permutation
    iu = np.triu_indices(G, 1); meds, nsig = [], []
    for Xb in cc.permuted_batches(X, DES, nperm, seed):                              # (c, G, N)
        r = cc._corr(Xb); mask = (np.abs(r) >= DES.rc_pooled)[:, iu[0], iu[1]]
        xc = Xb - Xb.mean(axis=-1, keepdims=True); tot = (xc @ np.swapaxes(xc, -1, -2)) / N
        with np.errstate(divide='ignore', invalid='ignore'): sh = btw[iu][None, :] / tot[:, iu[0], iu[1]]
        for k in range(Xb.shape[0]):
            v = sh[k][mask[k]]; nsig.append(int(v.size)); meds.append(float(np.median(v)) if v.size else float('nan'))
    meds = np.array(meds); nsig = np.array(nsig); ok = ~np.isnan(meds)
    out = dict(obs_median=obs['median'], obs_mean=obs['mean'], obs_n_pairs=obs['n_pairs'], null_n_pairs_mean=float(nsig.mean()), null_no_sig_frac=float(1 - ok.mean()))
    if ok.any():
        out.update(null_mean=float(meds[ok].mean()), null_q025=float(np.percentile(meds[ok], 2.5)), null_q975=float(np.percentile(meds[ok], 97.5)))
        out['excess'] = float(obs['median'] - out['null_mean']) if not np.isnan(obs['median']) else float('nan')
        out['above_q975'] = bool(obs['median'] > out['null_q975']) if not np.isnan(obs['median']) else None
    else:
        out.update(null_mean=float('nan'), null_q025=float('nan'), null_q975=float('nan'), excess=float('nan'), above_q975=None)
    return out


def tr2_after_pc1(Rb):
    """T_r2 of a stack (c, G, N) after subtracting the first singular component of every matrix; also the share of the residual variance in PC1"""
    U, S, Vt = np.linalg.svd(Rb, full_matrices=False)
    R1 = Rb - (U[..., :, :1] * S[..., None, :1]) @ Vt[..., :1, :]
    r = cc._corr(R1); G = r.shape[-1]; r2 = r * r
    return (r2.sum((1, 2)) - np.einsum('bii->b', r2)) / (G * (G - 1)), (S[..., 0] ** 2) / (S ** 2).sum(-1)


def coupling_pc1(X, nperm, seed):
    X = np.asarray(X, float); R = DES.Q @ np.log(X)                                  # N x G residuals, as coupling_test
    obs0 = float(cc._coupling_stats(R.T[None], DES, None, eig=False)['T_r2'][0]); o1, vs = tr2_after_pc1(R.T[None]); obs1 = float(o1[0]); pc1_share = float(vs[0])
    n0, n1 = [], []
    for Rb in cc.permuted_batches(R, DES, nperm, seed):                              # same permutations as the original test (same seed)
        n0.append(cc._coupling_stats(Rb, DES, None, eig=False)['T_r2']); n1.append(tr2_after_pc1(Rb)[0])
    n0 = np.concatenate(n0); n1 = np.concatenate(n1)
    p0 = cc._pvals(n0, obs0)['p_ge_plus1']; p1 = cc._pvals(n1, obs1)['p_ge_plus1']
    return dict(T_r2=obs0, T_r2_p=p0, T_r2_pc1=obs1, T_r2_pc1_p=p1, pc1_share=pc1_share, null_pc1_mean=float(n1.mean()))


def wilson(k, n):
    from statsmodels.stats.proportion import proportion_confint
    lo, hi = proportion_confint(int(k), int(n), method='wilson'); return float(lo), float(hi)


def analyse(matrices, label, stored=None):
    """matrices: list of (gse, platform, tag, X); stored: {(gse, tag): flat} for the gates"""
    rows = []; t0 = time.time(); gate = dict(share_max_abs_diff=0.0, tr2_p_mismatch=0, n_checked=0)
    for gse, plat, tag, X in matrices:
        seed = 2026 + int(gse[3:])
        a = share_expected(X, NPERM, seed); b = coupling_pc1(X, NPERM, seed)
        row = dict(gse=gse, platform=plat, set=tag, **a, **b); rows.append(row)
        if stored is not None and (gse, tag) in stored:
            f = stored[(gse, tag)]; gate['n_checked'] += 1
            if not (np.isnan(a['obs_median']) and (f.get('share_between') is None or np.isnan(f['share_between']))):
                gate['share_max_abs_diff'] = max(gate['share_max_abs_diff'], abs(a['obs_median'] - f['share_between']))
            if abs(b['T_r2_p'] - f['T_r2_p']) > 1e-12: gate['tr2_p_mismatch'] += 1
    print('%s: %d matrices in %.0f s' % (label, len(rows), time.time() - t0), flush=True)
    return rows, gate


def summarize(rows, label, gate):
    S = {'label': label, 'methods_sha256': sha256(METHODS), 'nperm': NPERM, 'gate': gate, 'sets': {}}
    for tag in ('A', 'B', 'A+B'):
        R = [r for r in rows if tag == 'A+B' or r['set'] == tag]; n = len(R)
        defined = [r for r in R if not np.isnan(r['obs_median']) and not np.isnan(r['null_mean'])]
        q = lambda v: [float(x) for x in np.percentile(v, [25, 50, 75])] if len(v) else [float('nan')] * 3
        f0 = sum(r['T_r2_p'] < 0.05 for r in R); f1 = sum(r['T_r2_pc1_p'] < 0.05 for r in R)
        both = sum(r['T_r2_p'] < 0.05 and r['T_r2_pc1_p'] < 0.05 for r in R); only0 = f0 - both; only1 = f1 - both
        S['sets'][tag] = dict(n_matrices=n, n_share_defined=len(defined),
                              obs_share_quartiles=q([r['obs_median'] for r in defined]), null_share_quartiles=q([r['null_mean'] for r in defined]), excess_quartiles=q([r['excess'] for r in defined]),
                              n_above_null_q975=int(sum(1 for r in defined if r['above_q975'])), pc1_share_quartiles=q([r['pc1_share'] for r in R]),
                              flagged_original=f0, flagged_original_ci=wilson(f0, n), flagged_pc1=f1, flagged_pc1_ci=wilson(f1, n),
                              paired=dict(both=both, only_original=only0, only_after_removal=only1, neither=n - f0 - only1))
    return S


def write_md(S, path):
    L = ['# Post hoc checks of the benchmark: results (%s)' % S['label'], '', 'Method: `POSTHOC_METHODS_191126.md` (SHA-256 %s, self-recorded), %d permutations per matrix, seed 2026 + GSE number.' % (S['methods_sha256'], S['nperm']), '',
         'Gates: %d matrices compared with the stored results; largest difference in the observed share %.2e; %d mismatches of the stored T_r2 P.' % (S['gate']['n_checked'], S['gate']['share_max_abs_diff'], S['gate']['tr2_p_mismatch']), '']
    for tag, v in S['sets'].items():
        qf = lambda x: '%.3f / %.3f / %.3f' % tuple(x)
        L += ['## Gene set %s (%d matrices)' % (tag, v['n_matrices']), '',
              '- Observed share of pooled covariance carried by group means, quartiles 25 / 50 / 75 (%d matrices with a defined value): %s' % (v['n_share_defined'], qf(v['obs_share_quartiles'])),
              '- Expected without within-group coupling (mean of the null medians), quartiles: %s' % qf(v['null_share_quartiles']),
              '- Excess (observed minus expected), quartiles: %s; observed above the 97.5th null percentile in %d of %d' % (qf(v['excess_quartiles']), v['n_above_null_q975'], v['n_share_defined']),
              '- Tr2 flagged (P < 0.05): original %d of %d (95%% CI %.2f to %.2f); after removal of PC1 %d of %d (%.2f to %.2f)' % (v['flagged_original'], v['n_matrices'], v['flagged_original_ci'][0], v['flagged_original_ci'][1], v['flagged_pc1'], v['n_matrices'], v['flagged_pc1_ci'][0], v['flagged_pc1_ci'][1]),
              '- Paired: flagged both times %d, only before removal %d, only after removal %d, neither %d' % (v['paired']['both'], v['paired']['only_original'], v['paired']['only_after_removal'], v['paired']['neither']),
              '- Share of the residual variance carried by PC1, quartiles: %s' % qf(v['pc1_share_quartiles']), '']
    open(path, 'w', encoding='utf-8').write('\n'.join(L))
    return '\n'.join(L)


def main(argv):
    mode = argv[1] if len(argv) > 1 else 'dev'
    if mode == 'dev':
        C = os.path.join(W2, 'compact')
        prep = json.load(open(os.path.join(C, 'dev_prepared.json'), encoding='utf8')); npz = np.load(os.path.join(C, 'before_lean', 'dev_sets.npz'))
        stored = {}
        for l in open(os.path.join(C, 'dev_core_results.jsonl'), encoding='utf8'):
            j = json.loads(l); stored[(j['gse'], j['tag'])] = j['flat']
        mats = [(g, prep[g]['platform'], tag, npz['%s__%s' % (g, tag)]) for g in prep for tag in ('A', 'B')]
        label = 'development series (37), stored matrices'
    else:
        raise SystemExit('validation mode is added after the user permits the download')
    out = os.path.join(HERE, 'results_' + mode); os.makedirs(out, exist_ok=True)
    print('methods file sha256', sha256(METHODS), flush=True)
    rows, gate = analyse(mats, label, stored)
    if gate['share_max_abs_diff'] > 1e-9 or gate['tr2_p_mismatch'] != 0 or gate['n_checked'] != len(mats):
        json.dump(dict(gate=gate, n=len(mats)), open(os.path.join(out, 'GATE_FAILED.json'), 'w')); print('GATE FAILED', gate); sys.exit(2)
    keys = list(rows[0].keys())
    with open(os.path.join(out, 'per_matrix.csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    S = summarize(rows, label, gate); json.dump(S, open(os.path.join(out, 'summary.json'), 'w'), indent=1, default=float)
    print(write_md(S, os.path.join(out, 'summary.md')))


if __name__ == '__main__': main(sys.argv)
