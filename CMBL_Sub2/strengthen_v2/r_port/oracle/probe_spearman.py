"""probe_spearman.py - why does the Spearman consensus null differ between an isolated call and the full WP1 script? Prints null mean / SD / P for
(a) consensus_null alone, (b) after observed_summary of the same kind (as in the WP1 loop), (c) with the denominators given in the order of the MLSS port."""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..')); sys.path.insert(0, os.path.join(HERE, '..', '..'))
import coexp_checks as cc
G17 = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']
X = pd.read_excel(os.path.join(ROOT, 'R_17genes', 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx'))[G17].astype(float).values; grp = np.repeat(['CTRL', 'DAPT', 'DKK-1'], 4); DEN = ['Ins1', 'Kcnj11', 'Hes1']
f = lambda c: '%.4f sd %.5f P %.4f' % (c['null_mean'], c['null_sd'], c['p_ge'])
for kind in ('pearson', 'log', 'spearman'):
    a = cc.consensus_null(X, grp, DEN, 10000, 2026, kind=kind, genes=G17)
    ov = cc.observed_summary(X, grp, DEN, G17, kind=kind); b = cc.consensus_null(X, grp, DEN, 10000, 2026, kind=kind, genes=G17)
    c = cc.consensus_null(X, grp, ['Hes1', 'Kcnj11', 'Ins1'], 10000, 2026, kind=kind, genes=G17)
    print(kind, '| alone:', f(a), '| after observed_summary:', f(b), '| other denominator order:', f(c))
print('archived: log 22.111 sd 2.897 P 0.297 | spearman 16.424 sd 2.984 P 0.341')
