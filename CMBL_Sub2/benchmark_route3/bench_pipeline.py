"""Generalized pooled co-expression pipeline of the manuscript (PROTOCOL_prespecified.md, section 2), for any 12 x G matrix of linear expression values (3 groups x 4 replicates in file order)
and any three denominator genes. Same definitions as ../reanalysis_17genes/pipeline.py (consensus = p <= 0.05 in all four schemes; group-preserving permutation null on the log scale;
treatment adjustment by residualizing on group, df = 8, Benjamini-Hochberg within scheme); the permutation null is vectorized, so its random stream differs from the original loop (validated in validate_pipeline.py)."""
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

N, K = 12, 3
GROUP = np.repeat(np.arange(K), 4)
DESIGN = np.column_stack([np.ones(N), (GROUP == 1).astype(float), (GROUP == 2).astype(float)])
_T = stats.t.ppf(0.975, N - 2); RC10 = _T / np.sqrt(N - 2 + _T ** 2)            # |r| that gives p = 0.05 at df 10 (pooled analysis)
DF_ADJ = N - (K - 1) - 2                                                          # 8: treatment-adjusted partial correlation


def scheme_list(G, den):
    """(kept gene indices, denominator index or None) for the four schemes: raw, then ratio to each denominator gene (dropped from its own scheme)"""
    return [(list(range(G)), None)] + [([i for i in range(G) if i != d], d) for d in den]


def scheme_matrix(X, keep, d): return X[:, keep] if d is None else X[:, keep] / X[:, [d]]


def pooled(X, den, rc=RC10):
    """per scheme: significance (|r| >= rc) and sign matrices embedded in the full G x G gene space (False where a gene is absent from the scheme)"""
    G = X.shape[1]; S, SG, cnt = [], [], []
    for keep, d in scheme_list(G, den):
        M = scheme_matrix(X, keep, d)
        with np.errstate(invalid='ignore', divide='ignore'): r = np.corrcoef(M, rowvar=False)
        r = np.nan_to_num(r, nan=0.0); np.fill_diagonal(r, 0.0)
        sig = np.zeros((G, G), bool); sgn = np.zeros((G, G)); ix = np.ix_(keep, keep)
        sig[ix] = np.abs(r) >= rc; sgn[ix] = np.sign(r); S.append(sig); SG.append(sgn); cnt.append(int(np.triu(sig, 1).sum()))
    return S, SG, cnt


def consensus_from(S, SG):
    """pairs significant in all four schemes (any sign) and the subset with the same sign in all four"""
    allsig = np.logical_and.reduce(S); same = allsig & (np.abs(np.sum(SG, axis=0)) == len(SG))
    return np.triu(allsig, 1), np.triu(same, 1)


def consensus_counts(X, den):
    S, SG, _ = pooled(X, den); a, s = consensus_from(S, SG); return int(a.sum()), int(s.sum())


def permutation_counts(X, den, nperm, rng):
    """null distribution of the consensus count: keep each gene's group means, permute its within-group log-scale residuals independently across genes"""
    L = np.log(X); gm = np.vstack([L[GROUP == g].mean(0) for g in range(K)]); resid = L - gm[GROUP]
    ca, cs = np.empty(nperm, int), np.empty(nperm, int); idx = [np.where(GROUP == g)[0] for g in range(K)]
    for b in range(nperm):
        R = resid.copy()
        for ix in idx: R[ix] = rng.permuted(R[ix], axis=0)
        ca[b], cs[b] = consensus_counts(np.exp(gm[GROUP] + R), den)
    return ca, cs


def adjusted(X, den):
    """treatment-adjusted partial correlations per scheme: q-value matrices in the full gene space (nan where absent) and q < 0.05 counts"""
    G = X.shape[1]; Q, cnt = [], []
    for keep, d in scheme_list(G, den):
        M = scheme_matrix(X, keep, d); R = M - DESIGN @ np.linalg.lstsq(DESIGN, M, rcond=None)[0]
        with np.errstate(invalid='ignore', divide='ignore'):
            r = np.corrcoef(R, rowvar=False); t = r * np.sqrt(DF_ADJ / (1 - r ** 2))
        p = 2 * stats.t.sf(np.abs(t), DF_ADJ); m = len(keep); iu = np.triu_indices(m, 1); pv = np.nan_to_num(p[iu], nan=1.0)
        q = np.full((m, m), np.nan); q[iu] = multipletests(pv, method='fdr_bh')[1]
        full = np.full((G, G), np.nan); full[np.ix_(keep, keep)] = q; Q.append(full); cnt.append(int((q[iu] < 0.05).sum()))
    return Q, cnt


def between_share(X, sigmask):
    """median share of the pooled covariance carried by the group means, over the pairs in sigmask (raw linear scale)"""
    xm = np.vstack([X[GROUP == g].mean(0) for g in range(K)])[GROUP]; xc = X - X.mean(0); bc = xm - X.mean(0)
    cov_tot = (xc.T @ xc) / N; cov_btw = (bc.T @ bc) / N; ii, jj = np.where(sigmask)
    if len(ii) == 0: return np.nan
    return float(np.median(cov_btw[ii, jj] / cov_tot[ii, jj]))


def analyse(X, den, nperm=2000, seed=2026):
    """all endpoints of the protocol for one 12 x G matrix"""
    X = np.asarray(X, float); S, SG, cnt = pooled(X, den); cons_any, cons_same = consensus_from(S, SG)
    rng = np.random.default_rng(seed); ca, cs = permutation_counts(X, den, nperm, rng)
    Q, qcnt = adjusted(X, den); obs = int(cons_any.sum())
    qmin = np.nanmin(np.stack(Q), axis=0); qall = np.nanmax(np.stack(Q), axis=0)           # q in the best scheme / in the worst scheme (nan if the pair is absent from a scheme)
    adj_any = cons_any & (qmin < 0.05); adj_all4 = cons_any & (np.nan_to_num(qall, nan=1.0) < 0.05)
    allq4 = np.triu(np.nan_to_num(qall, nan=1.0) < 0.05, 1)                                   # any pair significant after adjustment in all four schemes
    return dict(observed=obs, observed_same_sign=int(cons_same.sum()), pooled_sig_per_scheme=cnt, null_mean=float(ca.mean()), null_sd=float(ca.std()), null_q05=float(np.percentile(ca, 5)), null_q95=float(np.percentile(ca, 95)),
                perm_p=float((1 + (ca >= obs).sum()) / (nperm + 1)), perm_p_plain=float((ca >= obs).mean()), z=float((obs - ca.mean()) / ca.std()) if ca.std() > 0 else float('nan'),
                null_mean_same_sign=float(cs.mean()), perm_p_same_sign=float((1 + (cs >= int(cons_same.sum())).sum()) / (nperm + 1)),
                adj_q05_per_scheme=qcnt, cons_adj_any=int(adj_any.sum()), cons_adj_all4=int(adj_all4.sum()), adj_pairs_all4=int(allq4.sum()),
                share_between_raw=between_share(X, np.triu(S[0], 1)))
