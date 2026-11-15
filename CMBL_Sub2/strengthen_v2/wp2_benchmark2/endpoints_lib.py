"""endpoints_lib.py - endpoints of protocol v2 sections 3.3, 3.4 and 3.6 from per-matrix results (flat dicts of bench2_core.analyse_matrix, one row per series and gene set).
Shared by apply_dev.py (development series) and summarize_validation.py (held-out arrays, RNA-seq). Nothing is tuned here: P < 0.05 everywhere, Wilson 95% intervals, medians with quartiles.
P conventions: the consensus-count endpoint E1 uses the plain fraction P(null >= observed) as in part 1 (perm_p); every other test uses (1 + k)/(n + 1) (columns *_p, adjq_*_p)."""
import numpy as np, pandas as pd
from statsmodels.stats.proportion import proportion_confint

GEO_EXPECTED = 0.6086          # (2/pi) arccos(0.576): share of pooled-significant raw pairs if group means carry the correlation (3 groups, df 10)
SCHEMES = ['raw', 'den1', 'den2', 'den3']
SCHEME_LABEL = {'raw': 'raw', 'den1': 'ratio to gene 1', 'den2': 'ratio to gene 2', 'den3': 'ratio to gene 3'}


def wilson(k, n):
    if n == 0: return (float('nan'), float('nan'))
    lo, hi = proportion_confint(int(k), int(n), method='wilson'); return float(lo), float(hi)


def prop(mask):
    mask = np.asarray(mask, bool); k, n = int(mask.sum()), int(mask.size); lo, hi = wilson(k, n)
    return dict(k=k, n=n, prop=(k / n if n else float('nan')), lo=lo, hi=hi)


def med(x):
    x = pd.Series(x, dtype=float).dropna()
    if not len(x): return dict(n=0, median=float('nan'), q1=float('nan'), q3=float('nan'))
    return dict(n=int(len(x)), median=float(x.median()), q1=float(x.quantile(.25)), q3=float(x.quantile(.75)))


def endpoints(T):
    """T: DataFrame, one row per series (ONE gene set). Returns a nested dict of all endpoints."""
    out = dict(n_series=int(len(T)))
    out['E1_count_not_above_null'] = prop(T.perm_p >= 0.05); out['E1_count_flagged'] = prop(T.perm_p < 0.05)
    out['coupling_T_r2'] = prop(T.T_r2_p < 0.05)
    if 'T_p_p' in T: out['coupling_T_p'] = prop(T.T_p_p < 0.05)
    if 'T_eig_p' in T: out['coupling_T_eig'] = prop(T.T_eig_p < 0.05)
    out['T_r2_observed'] = med(T.T_r2_obs); out['T_r2_null_mean'] = med(T.T_r2_null_mean)
    ob = T.observed > 0
    out['E2_retained_any_scheme'] = med(T.retained_any[ob]); out['E2_retained_all_schemes'] = med(T.retained_all4[ob]); out['E2_retained_raw_scheme'] = med(T.retained_raw[ob])
    out['E2_n_series_with_consensus_pairs'] = int(ob.sum()); out['E2_pairs_sig_all_after_adjustment'] = med(T.adj_pairs_all4)
    out['E2_raw_any_consensus_pair_retained'] = prop(T.adj_cons_raw[ob] > 0) if ob.sum() else None
    out['E3_share_between'] = med(T.share_between)
    out['raw_sig_share_group_effect'] = med(T.share_sig_group_effect); out['raw_sig_share_all_pairs'] = med(T.share_sig_all)
    ok = T.share_sig_group_effect.notna(); out['raw_sig_share_above_expectation'] = prop(T.share_sig_group_effect[ok] > GEO_EXPECTED)
    out['observed_consensus'] = med(T.observed); out['null_mean_consensus'] = med(T.null_mean); out['consensus_minus_null'] = med(T.observed - T.null_mean)
    for s in SCHEMES:
        o, nm, p = T['adjq_%s_obs' % s], T['adjq_%s_null_mean' % s], T['adjq_%s_p' % s]
        ratio = (o / nm).where(nm > 0)
        out['adj_' + s] = dict(share_P_lt_0_05=prop(p < 0.05), median_observed=med(o), median_null_mean=med(nm), median_ratio_obs_over_null=med(ratio), n_ratio_defined=int(ratio.notna().sum()),
                               median_excess=med(o - nm), n_series_observed_gt_0=int((o > 0).sum()))
    return out


def _row(arm, gset, key, scheme, d, kind='proportion'):
    if kind == 'proportion': return dict(arm=arm, gene_set=gset, endpoint=key, scheme=scheme, k=d['k'], n=d['n'], estimate=d['prop'], lo=d['lo'], hi=d['hi'])
    return dict(arm=arm, gene_set=gset, endpoint=key, scheme=scheme, k=np.nan, n=d['n'], estimate=d['median'], lo=d['q1'], hi=d['q3'])


def long_table(arm, gset, e):
    """flatten endpoints() into rows (proportions with Wilson limits, medians with quartiles as lo/hi)"""
    rows = []
    for key, d in e.items():
        if d is None or not isinstance(d, dict): continue
        if 'prop' in d: rows.append(_row(arm, gset, key, '', d))
        elif 'median' in d: rows.append(_row(arm, gset, key, '', d, 'median'))
        elif key.startswith('adj_'):
            s = key[4:]
            rows.append(_row(arm, gset, 'adj_share_P_lt_0_05', s, d['share_P_lt_0_05']))
            for sub in ('median_observed', 'median_null_mean', 'median_ratio_obs_over_null', 'median_excess'): rows.append(_row(arm, gset, 'adj_' + sub, s, d[sub], 'median'))
    return rows


def fmt_prop(d): return '%d/%d = %.0f%% (%.0f-%.0f%%)' % (d['k'], d['n'], 100 * d['prop'], 100 * d['lo'], 100 * d['hi']) if d and d['n'] else 'n/a'
def fmt_med(d, pct=False, nd=2):
    if not d or d['n'] == 0: return 'n/a'
    f = (lambda v: '%.0f%%' % (100 * v)) if pct else (lambda v: ('%.' + str(nd) + 'f') % v)
    return '%s [%s, %s] (n=%d)' % (f(d['median']), f(d['q1']), f(d['q3']), d['n'])
