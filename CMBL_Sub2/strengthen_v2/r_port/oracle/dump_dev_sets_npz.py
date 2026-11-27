"""dump_dev_sets_npz.py - oracle for the R port of WP2: the 74 saved 12 x 17 development matrices (wp2_benchmark2/compact/dev_sets.npz, gene sets A and B of the 37 analysed series)
as JSON (floats written with repr, so R reads back the same doubles), next to dev_prepared.json. Run: python -I dump_dev_sets_npz.py"""
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); W = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'wp2_benchmark2', 'compact')
Z = np.load(os.path.join(W, 'dev_sets.npz')); out = {k: Z[k].tolist() for k in Z.files}
json.dump(out, open(os.path.join(HERE, 'dev_sets_archive.json'), 'w'))
print(len(out), 'arrays', Z[Z.files[0]].shape)
