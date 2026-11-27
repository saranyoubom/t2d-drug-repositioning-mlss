"""dump_oracle_runstruct.py - 60 simulated data sets (12 genes, 3 x 4 samples) with the run-structure test of coexp_checks.py (999 permutations, seed = data set index),
for the R port. Writes runstruct_oracle.f64: per data set 144 values of X (samples x genes, row-major) then p_two_sided(r), D(r), null_sd(r), p_two_sided(abs_r)."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', '..'))
import coexp_checks as cc
rs = np.random.default_rng(123); runs = ['A'] * 5 + ['B'] * 4 + ['C'] * 3; out = []
for i in range(60):
    Xs, g = cc.simulate([4, 4, 4], 12, rs); res = cc.run_structure_test(Xs, g, runs, nperm=999, seed=i)
    out.append(np.concatenate([Xs.ravel(), [res['r']['p_two_sided'], res['r']['D'], res['r']['null_sd'], res['abs_r']['p_two_sided']]]))
np.concatenate(out).astype('<f8').tofile(os.path.join(HERE, 'runstruct_oracle.f64'))
print('wrote 60 data sets; mean p(r) %.4f' % np.mean([o[-4] for o in out]))
