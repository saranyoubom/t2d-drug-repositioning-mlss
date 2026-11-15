"""bench2_core.py - per-matrix analysis of protocol v2 section 3 for one 12 x 17 matrix (3 groups x 4 replicates, linear scale, the first three genes are the denominators),
shared by the development series, the held-out arrays and the RNA-seq series (same code path). Every statistic comes from coexp_checks (WP1); this module only organizes the calls.
Permutation counts: within-group coupling tests 2000 (3.1), consensus count 2000 (as in part 1), calibrated adjusted counts 500 (3.4). The random stream is seeded with 2026 + GSE number
(coexp_checks streams differ from bench_pipeline; observed values are identical, null means agree within Monte Carlo error, see gate_part1.py)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE)
if V2 not in sys.path: sys.path.insert(0, V2)
import numpy as np
import coexp_checks as cc

GROUPS = [0] * 4 + [1] * 4 + [2] * 4
DES = cc.Design(GROUPS)
DEN = [0, 1, 2]
NPAIR = 17 * 16 // 2
SCHEMES = ['raw', 'den1', 'den2', 'den3']


def analyse_matrix(X, genes, n, nperm_cp=2000, nperm_cons=2000, nperm_adj=500, eig=True):
    """all endpoints of 3.1, 3.3, 3.4 and 3.6 for one matrix; returns (flat dict of scalars, detail dict)"""
    X = np.asarray(X, float); seed = 2026 + int(n); names = [str(g) for g in genes]
    obs = cc.observed_summary(X, DES, DEN, names); S0 = np.triu(obs['_S'][0], 1)
    cons = obs['_cons_any']; Q = np.stack(obs['_Q']); n_cons = int(cons.sum())
    q_raw = Q[0]; cons_raw = int((cons & (q_raw < 0.05)).sum())                                     # E2 on the raw scheme alone
    qall = np.nanmax(Q, axis=0); pairs_all4 = int(np.triu(np.nan_to_num(qall, nan=1.0) < 0.05, 1).sum())  # part 1 definition: significant in every scheme that contains the pair
    sh = cc.group_mean_share(X, DES, S0); geo = cc.geometry_check(X, DES)
    cp = cc.coupling_test(X, DES, nperm=nperm_cp, seed=seed, eig=eig)
    cn = cc.consensus_null(X, DES, DEN, nperm=nperm_cons, seed=seed, genes=names)
    an = cc.adjusted_null(X, DES, DEN, nperm=nperm_adj, seed=seed, genes=names)
    flat = dict(observed=n_cons, observed_same_sign=obs['consensus_same_sign'], null_mean=cn['null_mean'], null_sd=cn['null_sd'], null_q95=cn['null_q95'], z=cn['z'], perm_p=cn['p_ge'], perm_p_plus1=cn['p_ge_plus1'],
                perm_p_same_sign=cn['same_sign']['p_ge'], pooled_raw=obs['pooled_sig_per_scheme'][0], pooled_den1=obs['pooled_sig_per_scheme'][1], pooled_den2=obs['pooled_sig_per_scheme'][2], pooled_den3=obs['pooled_sig_per_scheme'][3],
                adj_any=obs['cons_adj_any'], adj_all4=obs['cons_adj_all_schemes'], adj_pairs_all4=pairs_all4, adj_cons_raw=cons_raw,
                retained_any=(obs['cons_adj_any'] / n_cons if n_cons else float('nan')), retained_all4=(obs['cons_adj_all_schemes'] / n_cons if n_cons else float('nan')), retained_raw=(cons_raw / n_cons if n_cons else float('nan')),
                share_between=sh['median'], share_sig_all=geo['share_sig_all'], share_sig_group_effect=geo['share_sig_group_effect'], n_group_effect_genes=geo['n_group_effect_genes'], geo_expected=geo['expected'])
    sdw = np.sqrt(((DES.Q @ np.log(X)) ** 2).sum(0) / (DES.N - DES.K))                                 # pooled within-group SD (natural log) per gene
    flat['n_zero_var_genes'] = int((sdw < 1e-6).sum())                                                  # degenerate genes (ADD-ON, 5 Nov 01:35): no within-group variation, residual correlation undefined
    for k in ('T_r2', 'T_p', 'T_eig'):
        if k in cp: flat.update({k + '_obs': cp[k]['observed'], k + '_null_mean': cp[k]['null_mean'], k + '_p': cp[k]['p_ge_plus1']})
    for s, nm in zip(SCHEMES, an['scheme_names']):
        v = an['schemes'][nm]; flat.update({'adjq_%s_obs' % s: v['observed'], 'adjq_%s_null_mean' % s: v['null_mean'], 'adjq_%s_null_sd' % s: v['null_sd'], 'adjq_%s_null_q95' % s: v['null_q95'],
                                            'adjq_%s_p' % s: v['p_ge_plus1'], 'adjq_%s_p_plain' % s: v['p_ge']})
    detail = dict(genes=names, consensus_pairs=obs['consensus_pairs'], degree=obs['degree'], adj_q05_per_scheme=obs['adj_q_lt_alpha_per_scheme'])
    return flat, detail
