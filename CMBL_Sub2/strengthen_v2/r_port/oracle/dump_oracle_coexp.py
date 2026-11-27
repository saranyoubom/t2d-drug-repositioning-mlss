"""dump_oracle_coexp.py - writes the raw null distributions of coexp_checks.py on the INS-1 matrix (seed 2026) as float64 files, so the R port can be compared
count by count (R_PORT_PLAN.md, gate 1.1). Run: python -I dump_oracle_coexp.py   (reads Analysis_outputs/R_17genes/Fig_4b_...xlsx, writes ./ *.f64)"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(HERE, '..', '..'))
import coexp_checks as cc
G17 = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']
SRC = os.environ.get('INS1_MATRIX', os.path.join(ROOT, 'R_17genes', 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx'))
X = pd.read_excel(SRC)[G17].astype(float).values; grp = np.repeat(['CTRL', 'DAPT', 'DKK-1'], 4); den = ['Ins1', 'Kcnj11', 'Hes1']
def save(name, a): np.ascontiguousarray(a, dtype='<f8').tofile(os.path.join(HERE, name)); print(name, np.shape(a))
an = cc.adjusted_null(X, grp, den, nperm=2000, seed=2026, genes=G17); save('adjusted_null_2000.f64', an['_null'])          # (2000, 4) C order
cn = cc.consensus_null(X, grp, den, nperm=10000, seed=2026, genes=G17); save('consensus_null_10000.f64', cn['_null']); save('consensus_null_same_10000.f64', cn['_null_same'])
des = cc.Design(grp); R = des.Q @ np.log(X); acc = {'T_r2': [], 'T_p': [], 'T_eig': []}
for Rb in cc.permuted_batches(R, des, 10000, 2026, None):
    for k, v in cc._coupling_stats(Rb, des).items(): acc[k].append(v)
for k, v in acc.items(): save('coupling_%s_10000.f64' % k, np.concatenate(v))
ct = cc.coupling_test(X, grp, nperm=10000, seed=2026)
print({k: (v['observed'], v['null_mean'], v['p_ge_plus1']) for k, v in ct.items()})
print('adjusted null p_ge', {n: v['p_ge'] for n, v in an['schemes'].items()})
