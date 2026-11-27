"""probe_matrices.py - how do the two copies of the INS-1 matrix differ (SBA/Figures, 21 genes, used by the working scripts and the archived results; R_17genes, 17 genes, shipped in the repository)?"""
import os
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); AO = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
SBA = os.path.normpath(os.path.join(AO, '..', '..', 'SBA', 'Figures', 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx')); R17 = os.path.join(AO, 'R_17genes', 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx')
G17 = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']
a, b = pd.read_excel(SBA), pd.read_excel(R17)
print('SBA columns', len(a.columns), list(a.columns)); print('R_17genes columns', len(b.columns), list(b.columns)[:3], '...')
A, B = a[G17].astype(float).values, b[G17].astype(float).values
d = A != B
print('cells that differ (17 genes x 12 rows):', int(d.sum()), '| max abs difference %.3e' % np.abs(A - B).max(), '| max relative %.3e' % (np.abs(A - B) / np.abs(B)).max())
for j in np.flatnonzero(d.any(0)): print('  gene', G17[j], 'rows', np.flatnonzero(d[:, j]).tolist(), 'SBA', A[d[:, j], j].tolist(), 'R_17genes', B[d[:, j], j].tolist())
tiesA = {G17[j]: len(A[:, j]) - len(set(A[:, j].tolist())) for j in range(17)}; tiesB = {G17[j]: len(B[:, j]) - len(set(B[:, j].tolist())) for j in range(17)}
print('duplicate values per gene: SBA', {k: v for k, v in tiesA.items() if v}, '| R_17genes', {k: v for k, v in tiesB.items() if v})
