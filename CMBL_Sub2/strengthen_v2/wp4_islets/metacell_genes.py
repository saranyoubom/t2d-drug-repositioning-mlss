"""metacell_genes.py - which of the 17 genes (and which consensus pairs) are evaluable (nonzero in >= 50% of the metacells) in the metacell partitions of metacell.py. Same partition code as metacell.one_run (same
seeds, donors with >= 40 beta-cells, exact blocks, remainder dropped, >= 2 blocks per donor), no ranks and no null; the number of evaluable genes per run is checked against metacell/metacell_runs.csv (genes_evaluable).
Output: metacell/metacell_gene_evaluability.csv (dataset, size, gene, partitions in which the gene is evaluable of 20, cell-level evaluable) and metacell/metacell_pair_evaluability.csv.   Usage: python -I metacell_genes.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd
from scipy import sparse
import wp4_core as C
from run_datasets import load, TAGS
import metacell as MC

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'metacell')


def evaluable_genes(cells, size, i):
    rng = np.random.default_rng(5000 + i); blocks, bdon = [], []
    for d in dict.fromkeys(cells.donor.tolist()):
        idx = np.nonzero(cells.donor == d)[0]
        if len(idx) < MC.MIN_CELLS: continue
        nb = len(idx) // size
        if nb < 2: continue
        perm = rng.permutation(idx)
        for b in range(nb): blocks.append(perm[b * size:(b + 1) * size]); bdon.append(d)
    nm = len(blocks)
    if nm < 4: return None, nm
    M = sparse.csr_matrix((np.ones(nm * size), (np.repeat(np.arange(nm), size), np.concatenate(blocks))), shape=(nm, cells.n)); S = (M @ cells.B).tocsr(); det = np.asarray((S > 0).mean(axis=0)).ravel()
    return {g: bool(cells.label(g) in cells.gi and det[cells.gi[cells.label(g)]] >= MC.DET) for g in C.GENES}, nm


if __name__ == '__main__':
    runs = pd.read_csv(os.path.join(OUT, 'metacell_runs.csv')); cons, consset, redge, cen = C.consensus_tables(); allpairs = [(a, c) for i, a in enumerate(C.GENES) for c in C.GENES[i + 1:] if frozenset((a, c)) in consset]
    rows, prow, bad = [], [], 0
    for tag in TAGS:
        cells = load(tag); det0 = np.asarray((cells.B > 0).mean(axis=0)).ravel(); cell_ev = {g: bool(cells.label(g) in cells.gi and det0[cells.gi[cells.label(g)]] >= 0.10) for g in C.GENES}
        for size in MC.SIZES:
            cnt = {g: 0 for g in C.GENES}; pcnt = {p: 0 for p in allpairs}; n_ok = 0
            for i in range(MC.NPART):
                ev, nm = evaluable_genes(cells, size, i)
                if ev is None: continue
                n_ok += 1
                for g, v in ev.items(): cnt[g] += v
                for a, c in allpairs: pcnt[(a, c)] += ev[a] and ev[c]
                r = runs[(runs.tag == tag) & (runs['size'] == size) & (runs.partition == i)].iloc[0]
                if int(r.genes_evaluable) != sum(ev.values()) or int(r.metacells) != nm: bad += 1
            for g in C.GENES: rows.append(dict(tag=tag, size=size, gene=g, partitions_evaluable=cnt[g], partitions_analysed=n_ok, cell_level_evaluable=cell_ev[g]))
            for a, c in allpairs: prow.append(dict(tag=tag, size=size, gene1=a, gene2=c, partitions_evaluable=pcnt[(a, c)], partitions_analysed=n_ok))
        print(tag, 'done', flush=True); del cells
    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'metacell_gene_evaluability.csv'), index=False); pd.DataFrame(prow).to_csv(os.path.join(OUT, 'metacell_pair_evaluability.csv'), index=False)
    print('runs where the recomputed counts differ from metacell_runs.csv:', bad)
