"""metacell.py - protocol 4.5 (metacell sensitivity), details fixed in PRESPEC_WP4_choices.md section C.
Per dataset and metacell size (20 primary; 10 and 50 sensitivity): 20 random partitions (partition i uses numpy default_rng(5000 + i)); donors with >= 40 beta-cells; the cells of a donor are shuffled and cut into
consecutive blocks of exactly `size` cells, the remainder is dropped, a donor needs >= 2 blocks. UMI and counts are summed (then log2(CP10K + 1) or log2(CPM + 1) with the summed library size), RPKM is averaged
(log2(mean + 1)). Donor-centred, ranked, Spearman, no depth term. Evaluable gene = nonzero in >= 50% of the metacells. Matched null rebuilt on the metacell data (deciles of metacell detection among evaluable genes
that are not one of the 17, 2000 draws, seed 2026) with the same excess statistic E as the cell-level analysis (mean Fisher z of the evaluable consensus pairs minus the mean over the null pair sets).
Output: metacell/metacell_runs.csv (one row per dataset, size, partition), metacell/metacell_summary.csv (median and range over the partitions), compared with the cell-level 10% analysis of results/index.csv.
Usage: python -I metacell.py [tag ...]      tags as in run_datasets.py (default: all human datasets that exist, then mouse)"""
import os, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd
from scipy import sparse
import wp4_core as C
from run_datasets import load, TAGS

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'metacell'); os.makedirs(OUT, exist_ok=True)
SIZES, NPART, MIN_CELLS, DET = (20, 10, 50), 20, 40, 0.50


def fast_draws(pool, bins, pair_idx, nn, seed):
    """same draws as wp4_core.draws_main (rng.choice(a) is a[rng.integers(0, len(a))], tested below: identical values and identical generator state) at about one fifth of the cost"""
    rng = np.random.default_rng(seed); IA = np.empty((nn, len(pair_idx)), np.int64); IC = np.empty_like(IA); arrs = [(pool[bins[a]], pool[bins[c]]) for a, c in pair_idx]; ints = rng.integers
    for s in range(nn):
        for j, (pa, pc) in enumerate(arrs):
            ia = pa[ints(0, len(pa))]; ic = pc[ints(0, len(pc))]
            while ic == ia: ic = pc[ints(0, len(pc))]
            IA[s, j], IC[s, j] = ia, ic
    return IA, IC


def _selftest():
    r = np.random.default_rng(1); bins = r.integers(0, 10, 400); pool = {k: np.nonzero(bins == k)[0] for k in range(10)}; pidx = [(int(a), int(b)) for a, b in r.integers(0, 400, (6, 2))]
    for seed in (2026, 7):
        A = C.draws_main(pool, bins, pidx, 40, seed); B = fast_draws(pool, bins, pidx, 40, seed); assert (A[0] == B[0]).all() and (A[1] == B[1]).all(), 'fast_draws differs from wp4_core.draws_main'


_selftest()


def one_run(cells, size, i):
    rng = np.random.default_rng(5000 + i); blocks, bdon = [], []
    for d in dict.fromkeys(cells.donor.tolist()):
        idx = np.nonzero(cells.donor == d)[0]
        if len(idx) < MIN_CELLS: continue
        nb = len(idx) // size
        if nb < 2: continue
        perm = rng.permutation(idx)
        for b in range(nb): blocks.append(perm[b * size:(b + 1) * size]); bdon.append(d)
    nm = len(blocks)
    if nm < 4: return dict(metacells=nm, donors_used=len(set(bdon)))
    M = sparse.csr_matrix((np.ones(nm * size), (np.repeat(np.arange(nm), size), np.concatenate(blocks))), shape=(nm, cells.n))
    S = (M @ cells.B).tocsr(); tot = M @ cells.tot; det = np.asarray((S > 0).mean(axis=0)).ravel(); elig = np.nonzero(det >= DET)[0]
    D = S[:, elig].toarray()
    X = np.log2(D / tot[:, None] * 1e4 + 1.0) if cells.kind == 'umi' else np.log2(D / tot[:, None] * 1e6 + 1.0) if cells.kind == 'count' else np.log2(D / size + 1.0)
    Xc = C.center_by_donor(X, np.array(bdon)); Rk = C.ranks_of(Xc)
    sym = cells.label; labels_e = [cells.genes[j] for j in elig]; idx_e = {g: r for r, g in enumerate(labels_e)}; det_e = det[elig]
    cons, consset, redge, cen = C.consensus_tables(); own = {sym(g) for g in C.GENES}; bins, pool = C.make_pool(det_e, labels_e, own)
    ev_genes = [g for g in C.GENES if sym(g) in idx_e]
    pairs = [(a, c) for i1, a in enumerate(C.GENES) for c in C.GENES[i1 + 1:] if frozenset((a, c)) in consset and sym(a) in idx_e and sym(c) in idx_e]
    out = dict(metacells=nm, donors_used=len(set(bdon)), eligible_genes=len(elig), genes_evaluable=len(ev_genes), consensus_pairs_evaluable=len(pairs))
    if not pairs: return out
    pidx = [(idx_e[sym(a)], idx_e[sym(c)]) for a, c in pairs]; rho = np.array([float(Rk[a] @ Rk[c]) for a, c in pidx]); z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    IA, IC = fast_draws(pool, bins, pidx, C.NN, C.SEED_NULL); V = C.null_values(Rk, IA, IC); zm = np.arctanh(np.clip(V, -0.999999, 0.999999)).mean(axis=1)
    out.update(median_rho=float(np.median(rho)), null_median_rho=float(np.median(V, axis=1).mean()), mean_z=float(z.mean()), mean_z_null=float(zm.mean()), E=float(z.mean() - zm.mean()), SE_E=float(zm.std(ddof=1)))
    return out


def main(tags):
    rows = []
    for tag in tags:
        cells = load(tag); t0 = time.time()
        for size in SIZES:
            for i in range(NPART):
                r = one_run(cells, size, i); r.update(tag=tag, size=size, partition=i); rows.append(r)
            sub = pd.DataFrame([x for x in rows if x['tag'] == tag and x['size'] == size])
            print('%-15s size %2d | metacells %s donors %s | evaluable genes %s pairs %s | median rho %s | E %s' % (tag, size, sub.metacells.median(), sub.donors_used.median(), sub.get('genes_evaluable', pd.Series([np.nan])).median(),
                  sub.get('consensus_pairs_evaluable', pd.Series([np.nan])).median(), '%.3f' % sub.median_rho.median() if 'median_rho' in sub else '-', '%.3f' % sub.E.median() if 'E' in sub else '-'), flush=True)
        print('%s done in %.0f s' % (tag, time.time() - t0), flush=True); del cells
        pd.DataFrame(rows).to_csv(os.path.join(OUT, 'metacell_runs.csv'), index=False)
    R = pd.DataFrame(rows); S = []
    cell = pd.read_csv(os.path.join(HERE, 'results', 'index.csv')); cell = cell[(cell.variant == 'primary') & (cell['mode'] == 'unadjusted')].set_index('tag')
    for (tag, size), g in R.groupby(['tag', 'size'], sort=False):
        f = lambda c: (g[c].median(), g[c].min(), g[c].max()) if c in g and g[c].notna().any() else (np.nan, np.nan, np.nan)
        S.append(dict(tag=tag, size=size, partitions=len(g), partitions_with_pairs=int(g.get('consensus_pairs_evaluable', pd.Series(dtype=float)).fillna(0).gt(0).sum()),
                      donors_used_median=g.donors_used.median(), metacells_median=g.metacells.median(), **{f'{c}_{s}': v for c in ('genes_evaluable', 'consensus_pairs_evaluable', 'median_rho', 'null_median_rho', 'E', 'SE_E') for s, v in zip(('med', 'min', 'max'), f(c))},
                      cell_level_pairs=int(cell.loc[tag, 'consensus_pairs_evaluable']) if tag in cell.index else np.nan, cell_level_median_rho=cell.loc[tag, 'median_rho_consensus'] if tag in cell.index else np.nan))
    S = pd.DataFrame(S); S['gain_pairs_median'] = S.consensus_pairs_evaluable_med - S.cell_level_pairs; S.to_csv(os.path.join(OUT, 'metacell_summary.csv'), index=False); print(S.round(3).to_string(index=False))


if __name__ == '__main__':
    main(sys.argv[1:] or [t for t in TAGS if os.path.exists(os.path.join(HERE, 'results', f'index_{t}.csv'))])
