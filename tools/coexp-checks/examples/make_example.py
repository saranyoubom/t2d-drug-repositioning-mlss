"""make_example.py - writes a SYNTHETIC example data set (no real data) and runs the checks on it.
3 groups x 4 replicates, 17 genes G01..G17, independent log-normal data with group effects; one latent factor with within-group correlation 0.7 is injected into genes G04 to G09,
so the within-group coupling test should flag it. Seeds are fixed.
Usage (from the folder coexp-checks):
  python -I examples/make_example.py            writes examples/example_data.csv and examples/expected_report.md
  coexp-checks examples/example_data.csv --groups 4,4,4 --nperm 2000 --out result.json --report report.md   (after pip install .)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'src'))
import numpy as np
import coexp_checks as cc

GENES = ['G%02d' % i for i in range(1, 18)]
GROUPS = [0] * 4 + [1] * 4 + [2] * 4


def example():
    rng = np.random.default_rng(2026)
    X, _ = cc.simulate([4, 4, 4], 17, rng, effect_sd=1.0, noise_sd=0.3)
    des = cc.Design(GROUPS)
    return cc.inject_factor(cc.decouple(X, des, rng), des, [3, 4, 5, 6, 7, 8], 0.7, rng)


def write_csv(path):
    X = example(); ids = ['ctrl_%d' % i for i in range(1, 5)] + ['treatA_%d' % i for i in range(1, 5)] + ['treatB_%d' % i for i in range(1, 5)]
    with open(path, 'w', encoding='utf8') as f:
        f.write('sample,' + ','.join(GENES) + '\n')
        for i, row in zip(ids, X): f.write(i + ',' + ','.join('%.6g' % v for v in row) + '\n')


if __name__ == '__main__':
    write_csv(os.path.join(HERE, 'example_data.csv'))
    X = example(); res = cc.run_checks(X, GROUPS, genes=GENES, nperm=2000, run_nperm=2000, which=['observed', 'share', 'coupling', 'consensus_null', 'geometry'])
    open(os.path.join(HERE, 'expected_report.md'), 'w', encoding='utf8').write(cc.to_markdown(res, 'Example (synthetic data)'))
    print(cc.to_markdown(res, 'Example (synthetic data)'))
