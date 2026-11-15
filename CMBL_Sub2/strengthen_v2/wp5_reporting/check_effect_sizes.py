"""WP5a: independent check of effect_sizes.csv.
Recomputes, with other code paths, (1) the Welch intervals of the 34 gene contrasts with scipy.stats.ttest_ind(...).confidence_interval, (2) the pooled-variance
intervals with statsmodels OLS over the three groups, (3) the Benjamini-Hochberg q with a hand-written step-up routine. Prints the maximum absolute differences
(expected below 1e-5, the rounding of the CSV to 6 significant digits).
Usage: python -I check_effect_sizes.py     (CMBL_ROOT environment variable overrides the project root = three folders above this script)"""
import os, sys, csv
from pathlib import Path
import numpy as np
from scipy import stats
import statsmodels.api as sm
import openpyxl

sys.stdout.reconfigure(encoding='utf8')
HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get('CMBL_ROOT', str(HERE.parents[2])))
AF4 = ROOT / 'Submission_package_JTM' / 'Additional file 4.xlsx'

wb = openpyxl.load_workbook(AF4, read_only=True, data_only=True)
rows = [[c for c in r if c is not None] for r in wb['Fig_4c_Raw_expression_data'].iter_rows(values_only=True)]
wb.close()
h = next(i for i, r in enumerate(rows) if r and r[0] == 'Wnt5a')
genes, X = rows[h], np.array(rows[h + 1:h + 13], float)
main = [r for r in csv.DictReader(open(HERE / 'effect_sizes.csv', encoding='utf8')) if r['readout'] == 'RT-qPCR' and r['variant'] == 'as reported']
assert len(main) == 34
design = np.column_stack([np.ones(12), [0] * 4 + [1] * 4 + [0] * 4, [0] * 8 + [1] * 4]).astype(float)
d_welch = d_pool = 0.0
for r in main:
    x = X[:, genes.index(r['feature'])]
    y = np.log2(x / x[:4].mean())
    k = (4, 8) if r['contrast'].startswith('DAPT') else (8, 12)
    ci = stats.ttest_ind(y[k[0]:k[1]], y[:4], equal_var=False).confidence_interval(0.95)
    d_welch = max(d_welch, abs(ci.low - float(r['ci95_lo'])), abs(ci.high - float(r['ci95_hi'])), abs(y[k[0]:k[1]].mean() - y[:4].mean() - float(r['log2fc_geo'])))
    fit = sm.OLS(y, design).fit()
    c = fit.conf_int(0.05)[1 if r['contrast'].startswith('DAPT') else 2]
    d_pool = max(d_pool, abs(c[0] - float(r['ci95_lo_pooled3'])), abs(c[1] - float(r['ci95_hi_pooled3'])))


def bh_manual(p):
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    q = np.empty(m)
    prev = 1.0
    for rank in range(m, 0, -1):
        i = o[rank - 1]
        prev = min(prev, p[i] * m / rank)
        q[i] = prev
    return q


d_bh = 0.0
for con in ('DAPT vs CTRL', 'DKK-1 vs CTRL'):
    s = [r for r in main if r['contrast'] == con]
    d_bh = max(d_bh, float(np.max(np.abs(bh_manual([float(r['p_raw_recomputed']) for r in s]) - np.array([float(r['bh_q_17genes']) for r in s])))))
print('max abs difference, Welch interval and effect vs scipy: %.2e' % d_welch)
print('max abs difference, pooled interval vs statsmodels OLS: %.2e' % d_pool)
print('max abs difference, BH q vs hand-written step-up:       %.2e' % d_bh)
