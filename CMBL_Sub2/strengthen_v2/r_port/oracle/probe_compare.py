"""probe_compare.py ROOT - run the Spearman consensus null on the INS-1 matrix found under ROOT (an Analysis_outputs folder) with the coexp_checks.py found under ROOT/strengthen_v2.
Prints checksums of the matrix, the module path, and the first null counts, so two roots can be compared."""
import hashlib, os, sys
import numpy as np, pandas as pd
ROOT = os.path.abspath(sys.argv[1]); sys.path.insert(0, os.path.join(ROOT, 'strengthen_v2'))
import coexp_checks as cc
G17 = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']
path = os.path.join(ROOT, 'R_17genes', 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx'); xl = pd.read_excel(path)
X = xl[G17].astype(float).values; grp = np.repeat(['CTRL', 'DAPT', 'DKK-1'], 4); DEN = ['Ins1', 'Kcnj11', 'Hes1']
print('module:', cc.__file__); print('xlsx sha256', hashlib.sha256(open(path, 'rb').read()).hexdigest()[:16], '| columns', list(xl.columns)[:5], '... n =', len(xl.columns), '| rows', len(xl))
print('X sha256', hashlib.sha256(np.ascontiguousarray(X).tobytes()).hexdigest()[:16], '| sum %.17g' % X.sum(), '| flags C', X.flags.c_contiguous, 'F', X.flags.f_contiguous)
cn = cc.consensus_null(X, grp, DEN, 10000, 2026, kind='spearman', genes=G17)
print('spearman null mean %.5f sd %.5f P %.4f | first 12 counts %s' % (cn['null_mean'], cn['null_sd'], cn['p_ge'], cn['_null'][:12].tolist()))
Xc = np.ascontiguousarray(X); cn2 = cc.consensus_null(Xc, grp, DEN, 10000, 2026, kind='spearman', genes=G17)
print('same data made C-contiguous: mean %.5f' % cn2['null_mean'])
Xf = np.asfortranarray(X); cn3 = cc.consensus_null(Xf, grp, DEN, 10000, 2026, kind='spearman', genes=G17)
print('same data made F-contiguous: mean %.5f' % cn3['null_mean'])
