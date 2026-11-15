"""Code validation (PROTOCOL section 6): the generalized pipeline must reproduce the manuscript's INS-1 numbers from the 17-gene matrix before it is used on any benchmark dataset.
Usage: python -I validate_pipeline.py"""
import sys, os, json
import numpy as np, pandas as pd
sys.stdout.reconfigure(encoding='utf8')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))          # -I removes the script folder from the path
import bench_pipeline as B
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'R_17genes', 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx')
REF = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'reanalysis_17genes', 'pipeline_17genes.json')))
GENES17 = [g for g in ['Wnt2', 'Wnt2b', 'Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Hey1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r'] if g not in ('Hey1', 'Wnt2', 'Wnt2b', 'Ryr1')]
df = pd.read_excel(SRC)[GENES17].astype(float); X = df.values
den = [GENES17.index(g) for g in ('Ins1', 'Kcnj11', 'Hes1')]
o = B.analyse(X, den, nperm=10000, seed=2026)
ok = True
def check(name, got, want, tol=0):
    global ok; good = abs(got - want) <= tol; ok &= good; print(f'{"OK " if good else "FAIL"} {name}: got {got} | manuscript {want}')
check('consensus pairs (any sign)', o['observed'], len(REF['consensus'])); check('consensus pairs (same sign)', o['observed_same_sign'], len(REF['consensus_signcons']))
for (k, c), got in zip(REF['counts'].items(), o['pooled_sig_per_scheme']): check(f'pooled significant pairs, {k}', got, c['sig'])
for (k, v), got in zip(REF['adj_sig_q05'].items(), o['adj_q05_per_scheme']): check(f'adjusted q<0.05, {k}', got, int(v))
check('consensus pairs with q<0.05 in >= 1 scheme after adjustment', o['cons_adj_any'], len(REF['cons_adj_any'])); check('pairs significant in all four schemes after adjustment', o['adj_pairs_all4'], len(REF['adj_all_four']))
check('permutation null mean (10000 permutations; manuscript 2000)', round(o['null_mean'], 2), round(REF['perm']['mean'], 2), tol=0.4); check('permutation null SD', round(o['null_sd'], 2), round(REF['perm']['sd'], 2), tol=0.3)
check('P(null >= observed)', round(o['perm_p_plain'], 3), round(REF['perm']['p_ge_obs'], 3), tol=0.04)
check('same-sign null mean', round(o['null_mean_same_sign'], 2), round(REF['perm_signcons']['mean'], 2), tol=0.4)
print('share of pooled covariance carried by group means (raw, significant pairs): %.2f (manuscript: median 0.97)' % o['share_between_raw'])
print('VALIDATION', 'PASSED' if ok else 'FAILED')
sys.exit(0 if ok else 1)
