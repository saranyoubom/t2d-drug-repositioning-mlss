"""coexp_checks.py - calibrated checks for pooled co-expression analyses of small designed experiments
(protocol v2 sections 2 and 3.1; WP1 of PLAN_v2_strengthen_JTM.md). Permutations are stacked and correlated in batched matrix products, so a null of thousands of permutations takes
fractions of a second (the loop in benchmark_route3/bench_pipeline.py does about 50-64 consensus permutations/s). Measured at G = 17, N = 12 on the 8-core laptop (test_coexp_checks.py, section 6): T_r2 + T_p about
5000 permutations/s (3000/s with T_eig), consensus counting 2500/s (50x the loop), adjusted counts over 4 schemes 1800/s, Spearman consensus 400/s. Random streams differ from bench_pipeline:
observed values are identical, null means and P agree within Monte Carlo error.

DATA CONVENTIONS (every function)
  X       samples x genes array on the LINEAR scale, strictly positive (2^-dCq, CPM, intensity ...); rows in the order of `groups`.
  groups  one label per sample (any number of groups >= 2, any sizes >= 2) or a Design object. Residual df = N - K - 1 (partial correlation after removing group means).
  genes   optional list of gene names (default G1..Gn).   den  optional denominator genes (names or column indices): ratio schemes divide every other gene by each of them.
  runs    optional list of run labels, one per gene (None / '' = unresolved: excluded from run tests, singleton block in block nulls).
  Seeds   every random function takes `seed` (default 2026) and draws from numpy default_rng; results do not depend on the memory chunking (MAX_ELEMS).
  Schemes 'raw' plus one ratio scheme per denominator gene (that gene is dropped from its own scheme). Consensus = pooled |r| >= r_crit(N-2) in all schemes, so only pairs of non-denominator
          genes qualify. P(null >= obs) is returned plain (p_ge) and as (1 + k)/(n + 1) (p_ge_plus1); use the second as a test P.

PUBLIC API (import coexp_checks as cc)
  everything in one call      cc.run_checks(X, groups, genes=None, den=None, runs=None, nperm=2000, run_nperm=20000, seed=2026, which=None) -> JSON-ready dict with a 'warnings' list
                              cc.to_markdown(result) -> short report;  which = subset of cc.ALL = observed, share, coupling, consensus_null, adjusted_null, geometry, runs
  design                      cc.Design(groups) (.N .K .sizes .rows .df_pooled .df_resid .rc_pooled .rc_resid .Q); cc.rcrit(df, alpha=0.05)
  observed pipeline           cc.observed_summary(X, groups, den, genes, kind='pearson') -> dict: pooled significant counts per scheme, consensus pairs (any sign / same sign), degrees, BH counts per scheme
                              after removing group means, consensus pairs with adjusted q < 0.05 in >= 1 / all schemes; keys starting with '_' hold arrays (_S, _SG, _Q q-matrices, _cons_any ...)
                              kind: 'pearson' (linear scale, the manuscript), 'log' (Pearson on log values, ratio = difference of logs), 'spearman' (ranks of each scheme matrix)
  group-mean share            cc.group_mean_share(X, groups, mask=None, kind='pearson') -> dict(median, mean, pooled, n_pairs): share of the pooled covariance carried by group means
                              cc.group_mean_share_expected(X, groups, nperm=2000, seed=2026, kind='pearson') -> the same statistic against its value when within-group coupling is removed (version 1.1):
                              observed_median, null_mean, null_q025, null_q975, excess; a share near 1 is expected for genes selected for a group effect, so read the excess
  within-group coupling       cc.coupling_test(X, groups, nperm=2000, seed=2026, blocks=None, pair_masks=None, eig=True) -> {'T_r2','T_p','T_eig',...}: observed, null_mean/sd/q95, p_ge_plus1 (one-sided).
                              T_r2 = mean squared off-diagonal correlation of log-scale residuals (primary), T_p = pairs with p < 0.05 at the residual df, T_eig = largest eigenvalue (eig=False skips it, 3x faster);
                              pair_masks {name: GxG bool} add T_r2_<name> (e.g. same-run pairs); cc.residual_corr(X, groups) is the GxG matrix behind them
  consensus null             cc.consensus_null(X, groups, den, nperm=2000, seed=2026, blocks=None, kind='pearson') -> observed count, null mean/sd/q05/q95, z, p_ge, p_ge_plus1, same-sign version;
                              blocks=run labels gives the run-block null (replicate labels permuted jointly for all genes of a run, within group)
  calibrated adjusted counts  cc.adjusted_null(X, groups, den, nperm=2000, seed=2026, blocks=None) -> per scheme: observed count of BH q < 0.05, null mean/sd/q95, excess, p_ge, p_ge_plus1
  run structure               cc.run_structure_test(X, groups, runs, nperm=20000, seed=2026) -> D = mean residual r of same-run minus between-run pairs (r and |r|), label-permutation P, per run;
                              cc.same_run_mask(runs) -> (same, between) pair masks
  geometry                    cc.geometric_expectation(r_crit, K) = P(|r| >= r_crit) if group means carry the correlation with uniform orientations ((2/pi) arccos r_crit for K = 3);
                              cc.geometry_check(X, groups) compares it with the share of pooled-significant raw pairs among genes with an ANOVA group effect; cc.anova_group_effect(X, groups)
  simulation (protocol 3.2)   cc.decouple(X, groups, rng); cc.inject_factor(X, groups, idx, r, rng); cc.simulate(n_per_group, G, rng, ...) -> (X, group codes)
  null engine                 cc.permuted_batches(M, groups, nperm, seed, blocks=None) yields (c, G, N) stacks (genes x samples) of within-group sample permutations of M (N x G);
                              genes of one block share the permutation (blocks=None: every gene its own block)
CLI  python coexp_checks.py DATA.csv --groups 4,4,4 [--runs RUNS] [--denominators G1,G2,G3] [--nperm 2000] [--run-nperm 20000] [--seed 2026] [--out result.json] [--report report.md]
     DATA.csv: header = gene names, one row per sample, linear scale; a non-numeric first column is read as sample ids. --groups: sizes of consecutive blocks ("4,4,4"), one label per sample,
     or col:NAME for a column of the csv. --runs: one label per gene (comma separated, in column order) or a csv with columns gene,run (missing = unresolved).
     The exit code is always 0; failed checks are listed under 'warnings' in the json and the markdown report."""
import argparse, itertools, json, os, sys
from functools import lru_cache
import numpy as np
from scipy import stats, special

__version__ = '1.1 (module of 4 Nov 2026 plus group_mean_share_expected and two warnings, 19 Nov 2026)'
MAX_ELEMS = 1_000_000          # elements per working array in a permutation chunk (memory guard; does not change results)
ALPHA = 0.05


# ------------------------------------------------------------------------------------------------------------------ design
def rcrit(df, alpha=ALPHA):
    """smallest |r| with two-sided p <= alpha at df degrees of freedom"""
    t = stats.t.ppf(1 - alpha / 2, df); return float(t / np.sqrt(df + t * t))


class Design:
    """group structure: codes, sizes, row indices, df and critical |r| values; Q is the N x N matrix that removes group means (M @ Q for genes x samples arrays)"""
    def __init__(self, groups):
        g = np.asarray(groups)
        if g.ndim != 1: raise ValueError('groups must be one label per sample')
        self.labels = list(dict.fromkeys(g.tolist())); self.codes = np.array([self.labels.index(x) for x in g.tolist()])
        self.N, self.K = len(g), len(self.labels)
        if self.K < 2: raise ValueError('need at least two groups')
        self.sizes = np.bincount(self.codes, minlength=self.K)
        if (self.sizes < 2).any(): raise ValueError('every group needs at least two samples')
        self.rows = [np.flatnonzero(self.codes == k) for k in range(self.K)]
        self.df_pooled, self.df_resid = self.N - 2, self.N - self.K - 1
        if self.df_resid < 3: raise ValueError('residual df = N - K - 1 = %d is below 3' % self.df_resid)
        self.rc_pooled, self.rc_resid = rcrit(self.df_pooled), rcrit(self.df_resid)
        A = np.zeros((self.N, self.N))
        for r in self.rows: A[np.ix_(r, r)] = 1.0 / len(r)
        self.Q = np.eye(self.N) - A

    def __repr__(self): return 'Design(N=%d, K=%d, sizes=%s, df_resid=%d)' % (self.N, self.K, self.sizes.tolist(), self.df_resid)


def as_design(groups): return groups if isinstance(groups, Design) else Design(groups)


def _check_X(X):
    X = np.asarray(X, float)
    if X.ndim != 2: raise ValueError('X must be samples x genes')
    if not np.isfinite(X).all() or (X <= 0).any(): raise ValueError('X must be finite and strictly positive (linear scale); add a pseudocount for zeros')
    return X


def _names(G, genes): return [str(g) for g in genes] if genes is not None else ['G%d' % (i + 1) for i in range(G)]


def _den_idx(den, genes):
    if den is None: return []
    out = []
    for d in den:
        if isinstance(d, (int, np.integer)): out.append(int(d))
        else:
            if genes is None or d not in list(genes): raise ValueError('denominator %r not among the genes' % (d,))
            out.append(list(genes).index(d))
    return out


def _unresolved(b): return b is None or (isinstance(b, str) and b == '') or (isinstance(b, float) and np.isnan(b))


# ------------------------------------------------------------------------------------------------------------------ building blocks (arrays are genes x samples, ..., G, N)
def _corr(M, return_z=False):
    """Pearson correlation of the rows of M (..., m, N) over samples; constant rows give r = 0 (never significant)"""
    mu = M.mean(axis=-1, keepdims=True); Z = M - mu; ss = np.einsum('...i,...i->...', Z, Z)[..., None]; ok = ss > 1e-24 * (ss + M.shape[-1] * mu * mu)
    Z *= np.where(ok, 1.0 / np.sqrt(np.where(ok, ss, 1.0)), 0.0); r = Z @ np.swapaxes(Z, -1, -2)
    return (r, Z) if return_z else r


def _pval(r, df):
    with np.errstate(divide='ignore', invalid='ignore'): t = r * np.sqrt(df / (1.0 - r * r))
    return 2.0 * special.stdtr(df, -np.abs(t))


def _bh(p):
    """Benjamini-Hochberg adjusted p along the last axis (same values as statsmodels fdr_bh)"""
    m = p.shape[-1]; order = np.argsort(p, axis=-1); ps = np.take_along_axis(p, order, -1)
    q = np.minimum.accumulate((ps * m / np.arange(1, m + 1))[..., ::-1], axis=-1)[..., ::-1]
    out = np.empty_like(q); np.put_along_axis(out, order, np.minimum(q, 1.0), -1); return out


@lru_cache(maxsize=64)
def _bh_thresholds(df, m, alpha):
    """critical |r| for rank k (k = 1..m): p(|r|) = alpha * k / m"""
    t = stats.t.isf(alpha * np.arange(1, m + 1) / m / 2.0, df); return t / np.sqrt(df + t * t)


def _bh_count(absr, df, alpha=ALPHA):
    """number of pairs with BH q < alpha from |r| (..., m pairs) without computing p-values: largest k whose k-th largest |r| exceeds the rank-k threshold"""
    m = absr.shape[-1]; ok = -np.sort(-absr, axis=-1) > _bh_thresholds(df, m, alpha)
    return np.where(ok.any(axis=-1), m - np.argmax(ok[..., ::-1], axis=-1), 0)


def _transform(M, kind):
    if kind == 'pearson': return M
    if kind == 'log': return np.log(M)
    if kind == 'spearman': return stats.rankdata(M, axis=-1)
    raise ValueError("kind must be 'pearson', 'log' or 'spearman'")


def _schemes(G, den):
    """(kept gene indices, denominator index or None): raw, then ratio to each denominator (dropped from its own scheme)"""
    return [(list(range(G)), None)] + [([i for i in range(G) if i != d], d) for d in den]


def _scheme_matrix(Xg, keep, d):
    """Xg (..., G, N) -> scheme matrix (..., m, N)"""
    M = Xg[..., keep, :]; return M if d is None else M / Xg[..., [d], :]


def _scheme_names(genes, den): return ['raw'] + ['%s-normalized' % genes[d] for d in den]


def _group_means(M, des):
    """group means broadcast back to the sample axis for M (N, m)"""
    out = np.empty_like(M)
    for rows in des.rows: out[rows] = M[rows].mean(axis=0)
    return out


def _blocks(G, blocks):
    """block index per gene: None -> every gene its own block; labels -> genes with the same label share a block; unresolved labels are singletons"""
    if blocks is None: return np.arange(G), G
    out, key, nxt = np.empty(G, int), {}, 0
    for i, b in enumerate(list(blocks)):
        if _unresolved(b): out[i] = nxt; nxt += 1
        else:
            if b not in key: key[b] = nxt; nxt += 1
            out[i] = key[b]
    return out, nxt


@lru_cache(maxsize=16)
def _perm_table(n): return np.array(list(itertools.permutations(range(n))), dtype=np.intp)


def permuted_batches(M, groups, nperm, seed=2026, blocks=None, max_elems=None):
    """yield (c, G, N) stacks (genes x samples) of within-group sample permutations of M (N x G). Genes in one block share the permutation: blocks=None permutes every gene independently
    (the null of 'no within-group coupling'); blocks=run labels keeps the coupling inside a run and destroys it between runs. Groups of up to 8 samples draw from the table of all
    permutations, larger ones sort random keys. The random stream does not depend on the chunk size."""
    des = as_design(groups); M = np.asarray(M, float); N, G = M.shape
    if N != des.N: raise ValueError('groups and M disagree on the number of samples')
    bog, nb = _blocks(G, blocks); cols = np.arange(G)[None, :, None]; W = sum(nb if n <= 8 else n * nb for n in des.sizes)
    rng = np.random.default_rng(seed); per = max(1, (max_elems or MAX_ELEMS) // (N * max(G, nb))); done = 0
    while done < nperm:
        c = min(per, nperm - done); u = rng.random((c, W)); off = 0; out = np.empty((c, G, N), dtype=M.dtype)
        for rows, n in zip(des.rows, des.sizes):
            if n <= 8: tab = _perm_table(int(n)); perm = tab[(u[:, off:off + nb] * len(tab)).astype(np.intp)][:, bog, :]; off += nb
            else: perm = np.argsort(u[:, off:off + n * nb].reshape(c, nb, n), axis=2)[:, bog, :]; off += n * nb
            out[:, :, rows] = M[rows][perm, cols]
        yield out; done += c


def _pairs(G): return np.triu_indices(G, 1)


def _pvals(null, obs):
    null = np.asarray(null); n = len(null); k = int((null >= obs - 1e-12 * max(1.0, abs(obs))).sum()); return dict(p_ge=float(k / n), p_ge_plus1=float((1 + k) / (n + 1)))


def _desc(null): return dict(null_mean=float(np.mean(null)), null_sd=float(np.std(null)), null_q05=float(np.percentile(null, 5)), null_q95=float(np.percentile(null, 95)))


def residual_corr(X, groups, log=True):
    """G x G Pearson correlation of residuals after removing group means (log scale by default): the within-group correlation matrix"""
    X = _check_X(X); des = as_design(groups); return _corr((np.log(X) if log else X).T @ des.Q)


# ------------------------------------------------------------------------------------------------------------------ observed pipeline
def observed_summary(X, groups, den=None, genes=None, kind='pearson', alpha=ALPHA):
    """pooled significant pairs per scheme, consensus pairs, treatment-adjusted BH counts per scheme and their overlap (the manuscript pipeline for any group structure).
    Pooled significance: |r| >= r_crit(N-2). Adjusted: group means removed from each (transformed) scheme matrix, partial r at df N-K-1, BH within scheme, count of q < alpha."""
    X = _check_X(X); des = as_design(groups); N, G = X.shape
    if N != des.N: raise ValueError('groups and X disagree on the number of samples')
    names = _names(G, genes); den = _den_idx(den, genes); sch = _schemes(G, den); Xg = X.T; S, SG, cnt, Q, qcnt = [], [], [], [], []
    for keep, d in sch:
        T = _transform(_scheme_matrix(Xg, keep, d), kind); r = _corr(T); np.fill_diagonal(r, 0.0); m = len(keep)
        sig = np.zeros((G, G), bool); sgn = np.zeros((G, G)); ix = np.ix_(keep, keep); sig[ix] = np.abs(r) >= des.rc_pooled; sgn[ix] = np.sign(r)
        S.append(sig); SG.append(sgn); cnt.append(int(np.triu(sig, 1).sum()))
        iu = _pairs(m); q = np.full((m, m), np.nan); q[iu] = _bh(_pval(_corr(T @ des.Q)[iu], des.df_resid)); full = np.full((G, G), np.nan); full[ix] = q; Q.append(full); qcnt.append(int((q[iu] < alpha).sum()))
    allsig = np.logical_and.reduce(S); same = allsig & (np.abs(np.sum(SG, axis=0)) == len(SG)); cons, cons_same = np.triu(allsig, 1), np.triu(same, 1)
    Qs = np.stack(Q); present = ~np.isnan(Qs); qmin = np.where(present, Qs, np.inf).min(axis=0); qmax = np.where(present, Qs, -np.inf).max(axis=0); in_all = present.all(axis=0)
    adj_any = cons & (qmin < alpha); adj_all = cons & in_all & (qmax < alpha); adj_pairs_all = in_all & (qmax < alpha)
    deg = np.zeros(G, int)
    for i, j in zip(*np.nonzero(cons)): deg[i] += 1; deg[j] += 1
    pr = lambda mask: [[names[i], names[j]] for i, j in zip(*np.nonzero(mask))]
    return dict(scheme_names=_scheme_names(names, den), kind=kind, pooled_sig_per_scheme=cnt, consensus_any=int(cons.sum()), consensus_same_sign=int(cons_same.sum()),
                consensus_pairs=pr(cons), consensus_pairs_same_sign=pr(cons_same), degree={names[i]: int(deg[i]) for i in np.argsort(-deg, kind='stable') if deg[i] > 0},
                adj_q_lt_alpha_per_scheme=qcnt, cons_adj_any=int(adj_any.sum()), cons_adj_all_schemes=int(adj_all.sum()), adj_pairs_all_schemes=int(adj_pairs_all.sum()),
                cons_adj_any_pairs=pr(adj_any), n_eligible_pairs=int(in_all.sum()), _S=S, _SG=SG, _Q=Q, _cons_any=cons, _cons_same=cons_same, _den=den, _genes=names)


def group_mean_share(X, groups, mask=None, kind='pearson'):
    """share of the pooled covariance carried by the group means, per gene pair: cov(group means) / cov(all samples), raw scheme. mask: boolean G x G (upper triangle used); default = pairs with
    pooled |r| >= r_crit(N-2) under the same `kind`. Returns median (manuscript: 0.97 for INS-1), mean, pooled (ratio of sums) and the number of pairs."""
    X = _check_X(X); des = as_design(groups); Y = _transform(X.T, kind).T; N = Y.shape[0]
    if mask is None: mask = np.abs(_corr(Y.T)) >= des.rc_pooled
    ii, jj = np.nonzero(np.triu(np.asarray(mask, bool), 1))
    if len(ii) == 0: return dict(median=float('nan'), mean=float('nan'), pooled=float('nan'), n_pairs=0)
    xc = Y - Y.mean(0); bc = _group_means(Y, des) - Y.mean(0); tot = (xc.T @ xc) / N; btw = (bc.T @ bc) / N; sh = btw[ii, jj] / tot[ii, jj]
    return dict(median=float(np.median(sh)), mean=float(np.mean(sh)), pooled=float(btw[ii, jj].sum() / tot[ii, jj].sum()), n_pairs=int(len(ii)))


def group_mean_share_expected(X, groups, nperm=2000, seed=2026, kind='pearson'):
    """the share of group_mean_share against what independence would give (added in version 1.1, 19 Nov 2026). Null: every gene's values are permuted independently within its group (group means kept,
    within-group coupling removed), the pooled-significant pairs are recomputed in every permuted matrix and the median share is taken again. Returns the observed median, the mean of the null medians
    with its 2.5th and 97.5th percentiles, excess = observed - null mean, above_q975 (observed above the 97.5th null percentile), the mean number of pooled-significant pairs in the null and the share of
    permutations without any (left out of the summary). A share near 1 is expected without any within-group coupling when the genes were selected for a group effect: read the excess, not the share."""
    X = _check_X(X); des = as_design(groups); Y = _transform(X.T, kind).T; N, G = Y.shape
    obs = group_mean_share(X, des, kind=kind); bc = _group_means(Y, des) - Y.mean(0); btw = (bc.T @ bc) / N; iu = _pairs(G); meds, nsig = [], []
    for Yb in permuted_batches(Y, des, nperm, seed):                                    # (c, G, N)
        r = _corr(Yb); mask = (np.abs(r) >= des.rc_pooled)[:, iu[0], iu[1]]
        xc = Yb - Yb.mean(axis=-1, keepdims=True); tot = (xc @ np.swapaxes(xc, -1, -2)) / N
        with np.errstate(divide='ignore', invalid='ignore'): sh = btw[iu][None, :] / tot[:, iu[0], iu[1]]
        for k in range(Yb.shape[0]):
            v = sh[k][mask[k]]; nsig.append(int(v.size)); meds.append(float(np.median(v)) if v.size else float('nan'))
    meds = np.array(meds); nsig = np.array(nsig); ok = ~np.isnan(meds); nan = float('nan')
    out = dict(observed_median=obs['median'], observed_n_pairs=obs['n_pairs'], null_n_pairs_mean=float(nsig.mean()), null_no_sig_frac=float(1 - ok.mean()), nperm=int(nperm), kind=kind)
    if ok.any():
        out.update(null_mean=float(meds[ok].mean()), null_q025=float(np.percentile(meds[ok], 2.5)), null_q975=float(np.percentile(meds[ok], 97.5)))
        out.update(excess=float(obs['median'] - out['null_mean']) if obs['n_pairs'] else nan, above_q975=bool(obs['median'] > out['null_q975']) if obs['n_pairs'] else None)
    else: out.update(null_mean=nan, null_q025=nan, null_q975=nan, excess=nan, above_q975=None)
    return out


# ------------------------------------------------------------------------------------------------------------------ within-group coupling
def _coupling_stats(Rb, des, masks=None, eig=True):
    """T_r2, T_p, T_eig (and T_r2 over pair masks) for a stack of log-scale residual matrices Rb (c, G, N)"""
    r, Z = _corr(Rb, return_z=True); G = r.shape[-1]; r2 = r * r; up = np.triu(np.ones((G, G), bool), 1); out = {}
    out['T_r2'] = (r2.sum((1, 2)) - np.einsum('bii->b', r2)) / (G * (G - 1)); out['T_p'] = ((np.abs(r) >= des.rc_resid) & up).sum((1, 2)).astype(float)
    if eig: out['T_eig'] = np.linalg.eigvalsh(np.swapaxes(Z, 1, 2) @ Z if Z.shape[-1] < G else r)[:, -1]       # same non-zero eigenvalues from the smaller Gram matrix
    for name, mk in (masks or {}).items():
        mk = np.triu(np.asarray(mk, bool), 1); out['T_r2_' + name] = (r2 * mk).sum((1, 2)) / max(int(mk.sum()), 1)
    return out


def coupling_test(X, groups, nperm=2000, seed=2026, blocks=None, pair_masks=None, eig=True):
    """within-group coupling on log-scale residuals (group means removed), raw scheme. Null: residuals permuted within group independently per gene (blocks=None) or jointly within blocks.
    Statistics: T_r2 (primary), T_p (pairs with p < 0.05 at the residual df), T_eig. P is one-sided: p_ge_plus1 = (1 + #{null >= obs}) / (nperm + 1).
    pair_masks: {name: G x G boolean} adds T_r2_<name> (mean squared residual r over those pairs, e.g. same-run / between-run)."""
    X = _check_X(X); des = as_design(groups); R = des.Q @ np.log(X); obs = {k: float(v[0]) for k, v in _coupling_stats(R.T[None], des, pair_masks, eig).items()}; acc = {k: [] for k in obs}
    for Rb in permuted_batches(R, des, nperm, seed, blocks):
        for k, v in _coupling_stats(Rb, des, pair_masks, eig).items(): acc[k].append(v)
    null = {k: np.concatenate(v) for k, v in acc.items()}
    return {k: dict(observed=obs[k], **_desc(null[k]), **_pvals(null[k], obs[k])) for k in obs}


# ------------------------------------------------------------------------------------------------------------------ consensus and adjusted counts
def _consensus_batch(Xb, elig, den, rc, kind='pearson'):
    """consensus counts (any sign, same sign) for a stack Xb (c, G, N); only pairs among eligible (non-denominator) genes can be consensus pairs"""
    allsig = agree = pos0 = None
    for d in [None] + list(den):
        M = Xb[:, elig, :] if d is None else Xb[:, elig, :] / Xb[:, [d], :]; r = _corr(_transform(M, kind)); s, pos = np.abs(r) >= rc, r > 0
        allsig = s if allsig is None else allsig & s
        if pos0 is None: pos0, agree = pos, np.ones_like(s)
        else: agree &= pos == pos0                                          # significant pairs have r != 0, so the sign is r > 0
    up = np.triu(np.ones(allsig.shape[-2:], bool), 1)
    return (allsig & up).sum((1, 2)), (allsig & agree & up).sum((1, 2))


def consensus_null(X, groups, den=None, nperm=2000, seed=2026, blocks=None, kind='pearson', genes=None):
    """consensus count (pooled significant in all schemes) against the permutation null. Null: group means kept, within-group log residuals permuted independently per gene (blocks=None) or
    jointly within blocks (e.g. runs). P(null >= obs) is reported plain (p_ge) and as (1 + k)/(n + 1) (p_ge_plus1)."""
    X = _check_X(X); des = as_design(groups); G = X.shape[1]; den = _den_idx(den, genes); elig = [i for i in range(G) if i not in den]
    oa, os_ = _consensus_batch(X.T[None], elig, den, des.rc_pooled, kind); oa, os_ = int(oa[0]), int(os_[0]); ca, cs = [], []
    for Xb in permuted_batches(X, des, nperm, seed, blocks):
        a, s = _consensus_batch(Xb, elig, den, des.rc_pooled, kind); ca.append(a); cs.append(s)
    ca, cs = np.concatenate(ca), np.concatenate(cs); sd = float(np.std(ca))
    return dict(observed=oa, **_desc(ca), z=float((oa - ca.mean()) / sd) if sd > 0 else float('nan'), **_pvals(ca, oa), n_eligible_pairs=len(elig) * (len(elig) - 1) // 2,
                same_sign=dict(observed=os_, **_desc(cs), **_pvals(cs, os_)), nperm=int(nperm), kind=kind, _null=ca, _null_same=cs)


def _adjusted_batch(Xb, des, sch, alpha=ALPHA):
    """BH q < alpha counts per scheme (c, S) for a stack Xb (c, G, N): group means removed from each scheme matrix, partial r at df N-K-1, BH within scheme"""
    out = np.empty((Xb.shape[0], len(sch)), int)
    for s, (keep, d) in enumerate(sch):
        r = _corr(_scheme_matrix(Xb, keep, d) @ des.Q); iu = _pairs(len(keep)); out[:, s] = _bh_count(np.abs(r[:, iu[0], iu[1]]), des.df_resid, alpha)
    return out


def adjusted_null(X, groups, den=None, nperm=2000, seed=2026, blocks=None, genes=None, alpha=ALPHA):
    """calibrated adjusted counts: per scheme, the observed number of pairs with BH q < alpha after removing group means, against its distribution when the within-group residuals are
    independent across genes (blocks=None). Ratio schemes share a denominator, so their null count is not zero and their q values are not calibrated; the raw-scheme count is."""
    X = _check_X(X); des = as_design(groups); G = X.shape[1]; den = _den_idx(den, genes); sch = _schemes(G, den); names = _scheme_names(_names(G, genes), den)
    obs = _adjusted_batch(X.T[None], des, sch, alpha)[0]; nul = np.concatenate([_adjusted_batch(Xb, des, sch, alpha) for Xb in permuted_batches(X, des, nperm, seed, blocks)])
    return dict(scheme_names=names, nperm=int(nperm), alpha=alpha, schemes={n: dict(observed=int(obs[s]), **_desc(nul[:, s]), excess=float(obs[s] - nul[:, s].mean()), **_pvals(nul[:, s], obs[s]))
                                                                           for s, n in enumerate(names)}, _null=nul)


# ------------------------------------------------------------------------------------------------------------------ run structure
def same_run_mask(runs):
    """G x G boolean masks (upper triangle) of same-run and between-run pairs; pairs with an unresolved gene are in neither"""
    G = len(runs); ok = np.array([not _unresolved(b) for b in runs]); lab = np.array([str(b) for b in runs]); eq = lab[:, None] == lab[None, :]; both = ok[:, None] & ok[None, :]
    up = np.triu(np.ones((G, G), bool), 1); return eq & both & up, ~eq & both & up


def run_structure_test(X, groups, runs, nperm=20000, seed=2026, genes=None, R=None):
    """same-run versus between-run within-group correlation. r = Pearson correlation of log-scale residuals (group means removed), raw scheme (or pass a G x G matrix as R).
    D = mean r over same-run pairs minus mean r over between-run pairs (tag 'abs_r': the same for |r|). Null: run labels permuted across genes preserving run sizes; two-sided
    P = (1 + #{|D_null| >= |D_obs|}) / (nperm + 1). Per run k: D_k = mean r over pairs inside run k minus mean r over pairs between run k and the other runs (same permutations).
    Genes without a run label are excluded."""
    X = _check_X(X); des = as_design(groups); names = _names(X.shape[1], genes); r = residual_corr(X, des) if R is None else np.asarray(R)
    keep = [i for i, b in enumerate(runs) if not _unresolved(b)]; lab = [runs[i] for i in keep]; uniq = list(dict.fromkeys(lab)); code = np.array([uniq.index(x) for x in lab]); g = len(keep)
    rr = r[np.ix_(keep, keep)]; iu = _pairs(g); P = len(iu[0]); nk = np.bincount(code, minlength=len(uniq)); ni = nk * (nk - 1) // 2; nc = nk * (g - nk); nS = int(ni.sum()); nB = P - nS
    vals = dict(r=rr[iu], abs_r=np.abs(rr[iu])); rng = np.random.default_rng(seed); per = max(1, 2_000_000 // max(P, 1))

    def stat(cp):
        si, sj = cp[:, iu[0]], cp[:, iu[1]]; same = (si == sj).astype(float); out = {}
        for tag, v in vals.items():
            Dk = np.full((cp.shape[0], len(uniq)), np.nan)
            for k in range(len(uniq)):
                if ni[k] > 0 and nc[k] > 0: Dk[:, k] = (((si == k) & (sj == k)).astype(float) @ v) / ni[k] - (((si == k) ^ (sj == k)).astype(float) @ v) / nc[k]
            out[tag] = ((same @ v) / max(nS, 1) - ((1 - same) @ v) / max(nB, 1), Dk)
        return out
    obs = stat(code[None]); nul = {t: ([], []) for t in vals}; done = 0
    while done < nperm:
        c = min(per, nperm - done); cp = code[np.argsort(rng.random((c, g)), axis=1)]
        for tag, (d, dk) in stat(cp).items(): nul[tag][0].append(d); nul[tag][1].append(dk)
        done += c
    res = dict(n_genes=g, runs={u: int(nk[k]) for k, u in enumerate(uniq)}, n_same_pairs=nS, n_between_pairs=nB, nperm=int(nperm), excluded=[names[i] for i in range(len(names)) if i not in keep])
    same0 = code[iu[0]] == code[iu[1]]
    for tag, v in vals.items():
        d0 = float(obs[tag][0][0]); dn = np.concatenate(nul[tag][0]); dkn = np.concatenate(nul[tag][1])
        res[tag] = dict(mean_same=float(v[same0].mean()) if nS else float('nan'), mean_between=float(v[~same0].mean()) if nB else float('nan'), D=d0, null_sd=float(dn.std()),
                        p_two_sided=float((1 + (np.abs(dn) >= abs(d0) - 1e-12).sum()) / (nperm + 1)), per_run={})
        for k, u in enumerate(uniq):
            dk0 = float(obs[tag][1][0, k]); ok = not np.isnan(dk0); inside = (code[iu[0]] == k) & (code[iu[1]] == k); cross = (code[iu[0]] == k) ^ (code[iu[1]] == k)
            res[tag]['per_run'][u] = dict(n_genes=int(nk[k]), n_pairs_inside=int(ni[k]), mean_inside=float(v[inside].mean()) if ni[k] else float('nan'), mean_cross=float(v[cross].mean()) if nc[k] else float('nan'),
                                          D=dk0, p_two_sided=float((1 + (np.abs(dkn[:, k]) >= abs(dk0) - 1e-12).sum()) / (nperm + 1)) if ok else float('nan'))
    return res


# ------------------------------------------------------------------------------------------------------------------ geometry
def geometric_expectation(r_crit, K):
    """P(|r| >= r_crit) for two genes whose pooled correlation is carried entirely by the group means, if the orientations of the centred group-mean vectors are uniform. The centred
    group-mean vector lives in K-1 dimensions; P = I_{1 - r_crit^2}((K-2)/2, 1/2) (regularized incomplete beta) = (2/pi) arccos(r_crit) for K = 3 (0.61 at r_crit = 0.576) and 1 for K = 2."""
    return 1.0 if K < 3 else float(special.betainc((K - 2) / 2.0, 0.5, 1.0 - r_crit ** 2))


def anova_group_effect(X, groups, log=True):
    """one-way ANOVA per gene (log scale by default): F, p and Benjamini-Hochberg q across genes"""
    X = _check_X(X); des = as_design(groups); Y = np.log(X) if log else X; gm = _group_means(Y, des); N, K = des.N, des.K
    ssw = ((Y - gm) ** 2).sum(0); ssb = ((gm - Y.mean(0)) ** 2).sum(0); F = (ssb / (K - 1)) / (ssw / (N - K)); p = stats.f.sf(F, K - 1, N - K); return dict(F=F, p=p, q=_bh(p))


def geometry_check(X, groups, fdr=ALPHA, log=True):
    """geometric expectation versus the observed share of pooled-significant raw pairs (all pairs, and pairs of genes with an ANOVA group effect at BH q < fdr). An anchor, not a test."""
    X = _check_X(X); des = as_design(groups); G = X.shape[1]; up = np.triu(np.ones((G, G), bool), 1); sig = (np.abs(_corr(X.T)) >= des.rc_pooled) & up
    de = anova_group_effect(X, des, log)['q'] < fdr; both = np.outer(de, de) & up
    return dict(r_crit=des.rc_pooled, expected=geometric_expectation(des.rc_pooled, des.K), n_pairs_all=int(up.sum()), n_sig_all=int(sig.sum()), share_sig_all=float(sig.sum() / up.sum()),
                n_group_effect_genes=int(de.sum()), n_pairs_group_effect=int(both.sum()), n_sig_group_effect=int((sig & both).sum()),
                share_sig_group_effect=float((sig & both).sum() / both.sum()) if both.sum() else float('nan'))


# ------------------------------------------------------------------------------------------------------------------ simulation helpers (protocol 3.2 and the tests)
def decouple(X, groups, rng):
    """remove within-group coupling: permute each gene's log values independently within every group (group means kept); returns linear scale"""
    X = _check_X(X); des = as_design(groups); L = np.log(X)
    for rows in des.rows: L[rows] = np.take_along_axis(L[rows], np.argsort(rng.random((len(rows), L.shape[1])), axis=0), axis=0)
    return np.exp(L)


def inject_factor(X, groups, idx, r, rng):
    """add one latent factor (standard normal, centred within group, scaled to unit within-group variance) to the log values of genes idx, loading sqrt(r/(1-r)) x the gene's pooled
    within-group SD; two injected genes then have within-group correlation r. Returns linear scale."""
    X = _check_X(X); des = as_design(groups); L = np.log(X); f = rng.standard_normal(des.N)
    for rows in des.rows: f[rows] -= f[rows].mean()
    f /= np.sqrt((f ** 2).sum() / (des.N - des.K)); sd = np.sqrt(((des.Q @ L) ** 2).sum(0) / (des.N - des.K))
    for j in idx: L[:, j] += np.sqrt(r / (1.0 - r)) * sd[j] * f
    return np.exp(L)


def simulate(n_per_group, G, rng, effect_sd=1.0, noise_sd=0.3, ref_noise_sd=0.0, ref_gene=None):
    """independent log-normal data: groups of the given sizes (a list), G genes with group effects ~ N(0, effect_sd) and independent residuals ~ N(0, noise_sd) on the log scale.
    With ref_gene (an index or a list of indices) and ref_noise_sd > 0 those genes get extra independent noise of that SD, so ratios to them share noise (a stand-in for a shared
    reference or denominator). Returns (X linear, group codes)."""
    codes = np.repeat(np.arange(len(n_per_group)), n_per_group); mu = rng.normal(0, effect_sd, (len(n_per_group), G)); L = mu[codes] + rng.normal(0, noise_sd, (len(codes), G))
    if ref_gene is not None and ref_noise_sd > 0:
        for j in np.atleast_1d(ref_gene): L[:, j] += rng.normal(0, ref_noise_sd, len(codes))
    return np.exp(L), codes


# ------------------------------------------------------------------------------------------------------------------ all-in-one
ALL = ('observed', 'share', 'coupling', 'consensus_null', 'adjusted_null', 'geometry', 'runs')


def run_checks(X, groups, genes=None, den=None, runs=None, nperm=2000, run_nperm=20000, seed=2026, which=None):
    """everything for one data set; returns a JSON-ready dict (arrays dropped) with a 'warnings' list. which: subset of ALL (default: all; 'runs' needs run labels)."""
    X = _check_X(X); des = as_design(groups); N, G = X.shape; names = _names(G, genes); denx = _den_idx(den, genes); which = set(which or ALL); res = dict(version=__version__)
    res['meta'] = dict(N=N, K=des.K, group_sizes=des.sizes.tolist(), G=G, genes=names, denominators=[names[d] for d in denx], df_pooled=des.df_pooled, df_resid=des.df_resid, r_crit_pooled=des.rc_pooled,
                       r_crit_resid=des.rc_resid, nperm=nperm, run_nperm=run_nperm, seed=seed)
    W = []
    if (des.sizes < 3).any(): W.append('a group has fewer than 3 samples: within-group correlations rest on very few degrees of freedom')
    if des.df_resid < 6: W.append('residual df = %d: the within-group tests have little power (r_crit = %.2f)' % (des.df_resid, des.rc_resid))
    if G < 5: W.append('only %d genes: T_r2 and the network statistics are unstable' % G)
    sdw = np.sqrt(((des.Q @ np.log(X)) ** 2).sum(0) / (des.N - des.K)); zero = [names[i] for i in np.flatnonzero(sdw < 1e-6)]        # added in version 1.1 (same rule as bench2_core: pooled within-group SD < 1e-6)
    if zero: W.append('%d gene(s) have no within-group variation (%s%s): their residual correlations are numerical noise and inflate the within-group tests; remove them before interpreting T_r2' % (len(zero), ', '.join(zero[:5]), ', ...' if len(zero) > 5 else ''))
    if runs is not None and len(runs) != G: raise ValueError('runs needs one label per gene')
    use_runs = runs is not None and 'runs' in which
    if 'observed' in which or 'share' in which:
        obs = observed_summary(X, des, denx, names)
        if 'observed' in which: res['observed'] = {k: v for k, v in obs.items() if not k.startswith('_')}
        if 'share' in which:
            res['group_mean_share'] = gm = group_mean_share(X, des, np.triu(obs['_S'][0], 1)); res['group_mean_share_expected'] = ge = group_mean_share_expected(X, des, nperm, seed)
            if gm['median'] > 0.9:
                if np.isnan(ge['null_mean']): W.append('group means carry %.0f%% of the pooled covariance of significant raw pairs (no pooled-significant pair in the permuted matrices, so no expected share): a share near 1 reflects the group effects of the genes and is not by itself evidence about within-group coupling' % (100 * gm['median']))
                else: W.append('group means carry %.0f%% of the pooled covariance of significant raw pairs; without any within-group coupling %.0f%% would be expected (95%% null interval %.0f%% to %.0f%%), so a share near 1 reflects the group effects of the genes and is not by itself evidence about within-group coupling' % (100 * gm['median'], 100 * ge['null_mean'], 100 * ge['null_q025'], 100 * ge['null_q975']))
    if 'geometry' in which: res['geometry'] = geometry_check(X, des)
    if 'coupling' in which:
        pm = None
        if use_runs: sm, bm = same_run_mask(runs); pm = dict(same_run=sm, between_run=bm)
        res['coupling'] = ct = coupling_test(X, des, nperm, seed, pair_masks=pm)
        if ct['T_r2']['p_ge_plus1'] < ALPHA: W.append('T_r2 = %.3f is above the within-group permutation null (P = %.4f): within-group coupling exists; it does not by itself identify biology' % (ct['T_r2']['observed'], ct['T_r2']['p_ge_plus1']))
    if 'consensus_null' in which:
        cn = consensus_null(X, des, denx, nperm, seed, genes=names); res['consensus_null'] = {k: v for k, v in cn.items() if not k.startswith('_')}
        if cn['p_ge_plus1'] >= ALPHA: W.append('consensus count %d is within the permutation null (mean %.1f, P = %.2f): weak evidence of no coupling, because the count test is insensitive to injected coupling in public benchmark data; read T_r2' % (cn['observed'], cn['null_mean'], cn['p_ge_plus1']))
        if use_runs: cb = consensus_null(X, des, denx, nperm, seed, blocks=runs, genes=names); res['consensus_null_run_block'] = {k: v for k, v in cb.items() if not k.startswith('_')}
    if 'adjusted_null' in which:
        an = adjusted_null(X, des, denx, nperm, seed, genes=names); res['adjusted_calibrated'] = {k: v for k, v in an.items() if not k.startswith('_')}
        for n, v in an['schemes'].items():
            if n != 'raw' and v['observed'] > 0 and v['p_ge_plus1'] >= ALPHA: W.append('%s scheme: %d pairs with BH q < 0.05 is not above its null (mean %.1f, P = %.2f): BH is not calibrated for a ratio scheme (shared denominator)' % (n, v['observed'], v['null_mean'], v['p_ge_plus1']))
    if use_runs:
        rs = run_structure_test(X, des, runs, run_nperm, seed, names); res['run_structure'] = rs
        if rs['r']['p_two_sided'] < ALPHA: W.append('within-group r differs between same-run and between-run pairs (D = %.2f, P = %.4f): compatible with shared reference, cDNA or plate effects, not co-regulation' % (rs['r']['D'], rs['r']['p_two_sided']))
    res['warnings'] = W; return res


def to_markdown(res, title='Co-expression checks'):
    """short readable report of a run_checks result"""
    m = res['meta']; L = ['# ' + title, '', 'N = %d samples, %d groups (sizes %s), %d genes, denominators: %s; residual df = %d; %d permutations (run test %d), seed %d.' % (m['N'], m['K'], m['group_sizes'], m['G'], ', '.join(m['denominators']) or 'none', m['df_resid'], m['nperm'], m['run_nperm'], m['seed']), '']
    L += ['## Warnings', ''] + (['- ' + w for w in res['warnings']] or ['- none']) + ['']
    if 'group_mean_share' in res:
        g = res['group_mean_share']; L += ['## Group-mean share', '', 'Median share of the pooled covariance carried by group means (raw, pooled-significant pairs, n = %d): %.2f (mean %.2f).' % (g['n_pairs'], g['median'], g['mean'])]
        if 'group_mean_share_expected' in res:
            e = res['group_mean_share_expected']; L.append('Expected without within-group coupling (%d permutations, mean of the null medians): %.2f (2.5th to 97.5th percentile %.2f to %.2f); observed minus expected %+.3f.' % (e['nperm'], e['null_mean'], e['null_q025'], e['null_q975'], e['excess']))
        L.append('')
    if 'coupling' in res:
        L += ['## Within-group coupling (log-scale residuals, raw scheme)', '', '| statistic | observed | null mean | null 95th pct | P (one-sided) |', '|---|---|---|---|---|']
        for k, v in res['coupling'].items(): L.append('| %s | %.4g | %.4g | %.4g | %.4f |' % (k, v['observed'], v['null_mean'], v['null_q95'], v['p_ge_plus1']))
        L.append('')
    if 'consensus_null' in res:
        c = res['consensus_null']; L += ['## Consensus count', '', 'Observed %d (same sign %d) of %d eligible pairs; permutation null mean %.1f (SD %.1f), P(null >= obs) = %.3f.' % (c['observed'], c['same_sign']['observed'], c['n_eligible_pairs'], c['null_mean'], c['null_sd'], c['p_ge'])]
        if 'consensus_null_run_block' in res: b = res['consensus_null_run_block']; L.append('Run-block null (replicate labels permuted jointly within each run): mean %.1f (SD %.1f), P(null >= obs) = %.3f.' % (b['null_mean'], b['null_sd'], b['p_ge']))
        L.append('')
    if 'adjusted_calibrated' in res:
        L += ['## Calibrated adjusted counts (BH q < 0.05 after removing group means)', '', '| scheme | observed | null mean | null 95th pct | P(null >= obs) |', '|---|---|---|---|---|']
        for n, v in res['adjusted_calibrated']['schemes'].items(): L.append('| %s | %d | %.1f | %.0f | %.3f |' % (n, v['observed'], v['null_mean'], v['null_q95'], v['p_ge']))
        L.append('')
    if 'run_structure' in res:
        r = res['run_structure']['r']; L += ['## Run structure', '', 'Mean residual r: same-run %.3f, between-run %.3f; D = %.3f, permutation P (two-sided) = %.4f.' % (r['mean_same'], r['mean_between'], r['D'], r['p_two_sided']), '']
    if 'geometry' in res: g = res['geometry']; L += ['## Geometry', '', 'r_crit %.3f; expected share if group means carry the correlation: %.2f; observed share of significant raw pairs: all %.2f, genes with a group effect %.2f (n genes %d).' % (g['r_crit'], g['expected'], g['share_sig_all'], g['share_sig_group_effect'], g['n_group_effect_genes']), '']
    return '\n'.join(L)


# ------------------------------------------------------------------------------------------------------------------ command line
def _jsonable(o):
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.ndarray): return o.tolist()
    raise TypeError(type(o))


def _parse_groups(spec, N):
    parts = [s.strip() for s in spec.split(',')]
    if len(parts) == N: return parts
    try: sizes = [int(s) for s in parts]
    except ValueError: raise ValueError('--groups: %d labels for %d samples' % (len(parts), N))
    if sum(sizes) != N: raise ValueError('--groups sizes sum to %d, not %d samples' % (sum(sizes), N))
    return [i for i, n in enumerate(sizes) for _ in range(n)]


def main(argv=None):
    import pandas as pd
    ap = argparse.ArgumentParser(description='Calibrated checks for pooled co-expression analyses (see the docstring of this file).')
    ap.add_argument('data'); ap.add_argument('--groups', required=True); ap.add_argument('--runs', default=None); ap.add_argument('--denominators', default=None)
    ap.add_argument('--nperm', type=int, default=2000); ap.add_argument('--run-nperm', type=int, default=20000); ap.add_argument('--seed', type=int, default=2026)
    ap.add_argument('--out', default=None); ap.add_argument('--report', default=None); a = ap.parse_args(argv); df = pd.read_csv(a.data)
    if a.groups.startswith('col:'): groups = df.pop(a.groups[4:]).tolist()
    if not pd.api.types.is_numeric_dtype(df.iloc[:, 0]): df = df.iloc[:, 1:]
    if not a.groups.startswith('col:'): groups = _parse_groups(a.groups, len(df))
    genes = [str(c) for c in df.columns]; X = df.to_numpy(float); runs = None
    if a.runs:
        if os.path.exists(a.runs): m = pd.read_csv(a.runs); mp = dict(zip(m.iloc[:, 0].astype(str), m.iloc[:, 1])); runs = [None if (g not in mp or pd.isna(mp[g])) else str(mp[g]) for g in genes]
        else: runs = [s.strip() or None for s in a.runs.split(',')]
    den = [s.strip() for s in a.denominators.split(',')] if a.denominators else None
    res = run_checks(X, groups, genes, den, runs, a.nperm, a.run_nperm, a.seed); md = to_markdown(res, 'Co-expression checks: ' + os.path.basename(a.data))
    if a.out: json.dump(res, open(a.out, 'w'), indent=1, default=_jsonable)
    if a.report: open(a.report, 'w', encoding='utf8').write(md)
    if hasattr(sys.stdout, 'reconfigure'): sys.stdout.reconfigure(encoding='utf8')
    print(md); return 0


if __name__ == '__main__':
    sys.exit(main())
