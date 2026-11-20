"""Smoke tests of coexp_checks (self-contained: no other file of the study repository is needed).
These are small, fast tests of known answers on simulated data with fixed seeds. They are NOT the validation suite of the study: that suite (56 checks against an independent loop implementation and
against the INS-1 reference numbers) is test_coexp_checks.py in the study repository and needs that repository's other files.
Usage:  python -I tests/test_smoke.py        (about a minute; exit code 1 if a test fails)"""
import os, sys, io, json, tempfile, contextlib, unittest
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'src')); sys.path.insert(0, os.path.join(HERE, '..', 'examples'))     # -I removes the script folder from the path
import numpy as np
import coexp_checks as cc

GROUPS = [0] * 4 + [1] * 4 + [2] * 4


class SmokeTests(unittest.TestCase):
    def test_design(self):
        d = cc.Design(GROUPS)
        self.assertEqual((d.N, d.K, d.df_pooled, d.df_resid), (12, 3, 10, 8))
        self.assertAlmostEqual(d.rc_pooled, 0.576, places=3)

    def test_permutation_keeps_group_means(self):
        rng = np.random.default_rng(1); X, _ = cc.simulate([4, 4, 4], 17, rng); des = cc.Design(GROUPS)
        for Xb in cc.permuted_batches(X, des, 20, seed=3):
            for rows in des.rows:
                self.assertTrue(np.allclose(Xb[:, :, rows].mean(axis=2), X[rows].mean(axis=0)[None, :], rtol=1e-12))

    def test_independent_data_type_I_error_is_not_inflated_much(self):
        des = cc.Design(GROUPS); hits = 0; n = 150
        for k in range(n):
            rng = np.random.default_rng(1000 + k); X, _ = cc.simulate([4, 4, 4], 17, rng, effect_sd=1.0, noise_sd=0.3)
            hits += cc.coupling_test(X, des, nperm=199, seed=k, eig=False)['T_r2']['p_ge_plus1'] < 0.05
        self.assertLessEqual(hits / n, 0.12, 'type I error %.3f' % (hits / n))

    def test_injected_coupling_is_detected(self):
        des = cc.Design(GROUPS); hits = 0; n = 100
        for k in range(n):
            rng = np.random.default_rng(2000 + k); X, _ = cc.simulate([4, 4, 4], 17, rng); X = cc.inject_factor(cc.decouple(X, des, rng), des, [3, 4, 5, 6, 7, 8], 0.7, rng)
            hits += cc.coupling_test(X, des, nperm=199, seed=k, eig=False)['T_r2']['p_ge_plus1'] < 0.05
        self.assertGreaterEqual(hits / n, 0.6, 'detection %.2f' % (hits / n))

    def test_group_mean_share_is_near_one_when_group_effects_dominate(self):
        rng = np.random.default_rng(5); X, _ = cc.simulate([4, 4, 4], 17, rng, effect_sd=2.0, noise_sd=0.05)
        sh = cc.group_mean_share(X, GROUPS)
        self.assertGreater(sh['n_pairs'], 10); self.assertGreater(sh['median'], 0.95)

    def test_same_seed_same_result(self):
        rng = np.random.default_rng(7); X, _ = cc.simulate([4, 4, 4], 17, rng)
        a = cc.coupling_test(X, GROUPS, nperm=99, seed=11, eig=False)['T_r2']['p_ge_plus1']; b = cc.coupling_test(X, GROUPS, nperm=99, seed=11, eig=False)['T_r2']['p_ge_plus1']
        self.assertEqual(a, b)

    def test_constant_gene_is_reported(self):
        # known limitation of version 1.0 and 1.1: a gene without within-group variation gets residual correlations that are numerical noise; run_checks warns (the benchmark's type I error excess came from such genes)
        rng = np.random.default_rng(9); X, _ = cc.simulate([4, 4, 4], 17, rng); X[:, 0] = 5.0
        res = cc.run_checks(X, GROUPS, nperm=99, run_nperm=99, which=['observed', 'coupling'])
        self.assertTrue(any('no within-group variation' in w for w in res['warnings']))

    def test_expected_share_is_returned_and_the_warning_names_it(self):
        rng = np.random.default_rng(5); X, _ = cc.simulate([4, 4, 4], 17, rng, effect_sd=2.0, noise_sd=0.05)
        e = cc.group_mean_share_expected(X, GROUPS, nperm=199, seed=1)
        self.assertGreater(e['observed_median'], 0.95); self.assertGreater(e['null_mean'], 0.95); self.assertFalse(np.isnan(e['excess']))
        res = cc.run_checks(X, GROUPS, nperm=99, run_nperm=99, which=['observed', 'share'])
        self.assertTrue(any('would be expected' in w for w in res['warnings']))

    def test_run_checks_returns_report_parts(self):
        rng = np.random.default_rng(13); X, _ = cc.simulate([4, 4, 4], 17, rng)
        res = cc.run_checks(X, GROUPS, nperm=99, run_nperm=99, which=['observed', 'share', 'coupling'])
        for k in ('meta', 'observed', 'group_mean_share', 'coupling', 'warnings'): self.assertIn(k, res)
        self.assertIn('# Co-expression checks', cc.to_markdown(res))

    def test_command_line_on_the_example(self):
        import make_example
        with tempfile.TemporaryDirectory() as tmp:
            csv = os.path.join(tmp, 'example.csv'); make_example.write_csv(csv)
            out, rep = os.path.join(tmp, 'result.json'), os.path.join(tmp, 'report.md')
            with contextlib.redirect_stdout(io.StringIO()):
                cc.main([csv, '--groups', '4,4,4', '--nperm', '199', '--run-nperm', '199', '--seed', '1', '--out', out, '--report', rep])
            res = json.load(open(out, encoding='utf8'))
            self.assertIn('coupling', res); self.assertTrue(os.path.getsize(rep) > 100)


if __name__ == '__main__':
    unittest.main(verbosity=2)
