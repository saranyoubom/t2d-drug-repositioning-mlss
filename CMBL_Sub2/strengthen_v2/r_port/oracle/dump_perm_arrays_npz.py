"""Oracle dump: the archived perm_arrays_17.npz of reanalysis_17genes as a CSV (columns in the order of the archive), so that the R port can compare the arrays exactly.
Usage: python -I dump_perm_arrays_npz.py <perm_arrays_17.npz> <out.csv>"""
import sys
import numpy as np

z = np.load(sys.argv[1])
keys = list(z.files)
cols = [z[k] for k in keys]
n = len(cols[0])
with open(sys.argv[2], 'w', encoding='utf8', newline='') as f:
    f.write(','.join(keys) + '\n')
    for i in range(n):
        f.write(','.join(str(int(c[i])) for c in cols) + '\n')
print('wrote', sys.argv[2], keys, n)
