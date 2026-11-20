"""post_hoc_ins1_191126.py - main session, 19 Nov 2026. The two post hoc checks of post_hoc_checks_191126.py on the INS-1 matrix of the manuscript (addendum POSTHOC_METHODS_ADDENDUM_191126.md, written after the
development results were known). Imports the functions of post_hoc_checks_191126.py unchanged; reads the matrix file only (never writes next to it); writes into results_ins1/.
Usage (python -I): python -I post_hoc_ins1_191126.py   (env INS1_MATRIX overrides the matrix path)"""
import os, sys, json, hashlib
import numpy as np, pandas as pd
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import post_hoc_checks_191126 as P
ADD = os.path.join(HERE, 'POSTHOC_METHODS_ADDENDUM_191126.md')
MATRIX = os.environ.get('INS1_MATRIX', os.path.join(HERE, '..', '..', '..', 'R_17genes', 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx'))
G17 = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']
NPERM, SEED = 10000, 2026
print('main methods sha256', P.sha256(P.METHODS), '| addendum sha256', P.sha256(ADD), flush=True)
X = pd.read_excel(MATRIX)[G17].astype(float).values
assert X.shape == (12, 17)
a = P.share_expected(X, NPERM, SEED); b = P.coupling_pc1(X, NPERM, SEED)
gate = dict(share_ok=round(a['obs_median'], 3) == 0.971, tr2_ok=round(b['T_r2'], 3) == 0.133, p_ok=abs(b['T_r2_p'] - 0.061) <= 0.004)
print('gates', gate, '| observed share %.4f, T_r2 %.4f, P %.4f' % (a['obs_median'], b['T_r2'], b['T_r2_p']), flush=True)
out = os.path.join(HERE, 'results_ins1'); os.makedirs(out, exist_ok=True)
res = dict(main_methods_sha256=P.sha256(P.METHODS), addendum_sha256=P.sha256(ADD), nperm=NPERM, seed=SEED, gate=gate, share=a, coupling=b)
if not all(gate.values()):
    json.dump(res, open(os.path.join(out, 'GATE_FAILED.json'), 'w'), indent=1, default=float); print('GATE FAILED'); sys.exit(2)
json.dump(res, open(os.path.join(out, 'summary.json'), 'w'), indent=1, default=float)
L = ['# Post hoc checks on the INS-1 matrix (addendum; %d permutations, seed %d)' % (NPERM, SEED), '',
     'Methods: `POSTHOC_METHODS_191126.md` (SHA-256 %s) and the addendum (SHA-256 %s), both self-recorded.' % (res['main_methods_sha256'], res['addendum_sha256']), '',
     'Gates: observed share 0.971, T_r2 0.133 and P within 0.004 of 0.061 (manuscript): %s.' % gate, '',
     '## Check 1: share of the pooled covariance carried by group means',
     '- Observed (median over %d pooled-significant pairs): %.3f' % (a['obs_n_pairs'], a['obs_median']),
     '- Expected without within-group coupling (mean of the null medians, 10000 permutations): %.3f (2.5th to 97.5th percentile %.3f to %.3f); observed minus expected %+.3f; above the 97.5th percentile: %s' % (a['null_mean'], a['null_q025'], a['null_q975'], a['excess'], a['above_q975']),
     '- Pooled-significant pairs in the permuted matrices: mean %.1f (observed %d); permutations without a significant pair: %.3f' % (a['null_n_pairs_mean'], a['obs_n_pairs'], a['null_no_sig_frac']), '',
     '## Check 2: T_r2 before and after removal of the first principal component',
     '- Original: T_r2 %.4f, P = %.4f' % (b['T_r2'], b['T_r2_p']),
     '- After PC1 removal: T_r2 %.4f, P = %.4f (null mean %.4f)' % (b['T_r2_pc1'], b['T_r2_pc1_p'], b['null_pc1_mean']),
     '- Share of the residual variance carried by PC1: %.3f' % b['pc1_share'], '']
open(os.path.join(out, 'summary.md'), 'w', encoding='utf-8').write('\n'.join(L)); print('\n'.join(L))
