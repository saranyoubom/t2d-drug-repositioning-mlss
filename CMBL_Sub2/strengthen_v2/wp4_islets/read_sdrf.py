"""read_sdrf.py - print the column names and value counts of a MAGE-TAB SDRF file (sample / cell level metadata; no expression values).
Usage: python -I read_sdrf.py <E-MTAB-xxxx.sdrf.txt> [max_unique_to_show]"""
import sys
import pandas as pd
sys.stdout.reconfigure(encoding='utf8')
T = pd.read_csv(sys.argv[1], sep='\t', dtype=str, low_memory=False); nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 20
print('rows', len(T)); print('columns:'); [print('  ', i, c, '| unique', T[c].nunique(dropna=True)) for i, c in enumerate(T.columns)]
for c in T.columns:
    u = T[c].nunique(dropna=True)
    if 1 < u <= nmax: print('\n--', c, '\n', T[c].value_counts(dropna=False).to_string())
