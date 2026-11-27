"""dump_oracle_spearman.py - Spearman consensus null (10000 permutations, seed 2026) of coexp_checks.py on the INS-1 matrix, for tracing a difference of the R port.
Writes spearman_null.f64 (any-sign counts, then same-sign counts, 10000 each) and, for the first permutation where a count is recomputed per scheme, the matrices of that permutation."""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..')); sys.path.insert(0, os.path.join(HERE, '..', '..'))
import coexp_checks as cc
G17 = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']
X = pd.read_excel(os.path.join(ROOT, 'R_17genes', 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx'))[G17].astype(float).values; grp = np.repeat(['CTRL', 'DAPT', 'DKK-1'], 4); den = ['Ins1', 'Kcnj11', 'Hes1']
cn = cc.consensus_null(X, grp, den, nperm=10000, seed=2026, kind='spearman', genes=G17)
np.concatenate([cn['_null'], cn['_null_same']]).astype('<f8').tofile(os.path.join(HERE, 'spearman_null.f64'))
print('null mean %.4f sd %.4f p %.4f' % (cn['null_mean'], cn['null_sd'], cn['p_ge']))
