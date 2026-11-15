"""run_all_v2.py - re-creates the v2 tables and figures from cached inputs by calling each package's own scripts (no logic is copied here).
Usage: python -I run_all_v2.py [--dry-run] [--only wp1,wp2,wp3,wp4,wp5a] [--skip-heavy] [--list]
Paths are relative to this file; large downloads live under the environment variable V2_DATA_ROOT (default: a scratch folder). One heavy job at a time: packages run in sequence, not in parallel.
--dry-run only checks that every input exists and prints what would run. Every missing input gets its own message. Inputs marked DATA: live under V2_DATA_ROOT (downloads and parsed single-cell data of WP4);
without V2_DATA_ROOT they are not checked and the note says so. WP4 steps: run_datasets (cell-level analysis of the six islet datasets), meta_analysis, gene_graph_null (post hoc add-on), metacell (heavy), perturbation_gse60803,
make_tables, make_fig7_v2; the downloads and parsers that create the DATA: inputs are listed in wp4_islets/README.md section 14."""
import os, sys, subprocess, time
HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE)
ROOT = os.environ.get('V2_DATA_ROOT', '')
PY = [sys.executable, '-I']
# each step: (package, folder, script and arguments, inputs that must exist (relative to V2 unless marked DATA:), what it writes, rough minutes, heavy?)
STEPS = [
    ('wp1', 'wp1_ins1_structure', ['wp1_ins1_analysis.py'], ['wp1_ins1_structure/ins1_sheet_cq.csv', 'wp1_ins1_structure/gene_run_map.csv', 'coexp_checks.py'], 'tbl_2_*.csv, results_ins1.json', 20, True),
    ('wp3', 'wp3_mlss', ['mlss_scenarios.py'], ['wp3_mlss/mlss_port.py', 'wp3_mlss/direction_table.csv', 'wp3_mlss/pic50_sources.csv'], 'mlss_scenarios.csv, mlss_permutation_draws.csv', 5, True),
    ('wp5a', 'wp5_reporting', ['effect_sizes.py'], ['wp5_reporting/effect_sizes.py'], 'effect_sizes.csv/.md', 1, False),
    ('wp2', 'wp2_benchmark2', ['apply_dev.py'], ['wp2_benchmark2/compact/dev_core_results.jsonl', 'wp2_benchmark2/endpoints_lib.py'], 'dev_application_summary.md', 1, False),
    ('wp2', 'wp2_benchmark2', ['summarize_validation.py'], ['wp2_benchmark2/dev_application_per_series.csv', 'wp2_benchmark2/heldout_arrays_results.jsonl', 'wp2_benchmark2/heldout_rnaseq_results.jsonl'], 'validation_summary.md, WP2_summary_table.md', 1, False),
    ('wp2', 'wp2_benchmark2', ['calibrate_t2.py', 'summarize'], ['wp2_benchmark2/compact/calib_series'], 'calibration_tables.csv, calibration_verdict.txt', 1, False),
    ('wp2', 'wp2_benchmark2', ['design_table.py', 'summarize'], ['wp2_benchmark2/compact/design_cells'], 'design_table.md', 1, False),
    ('wp2', 'wp2_benchmark2', ['make_wp2_figure.py'], ['wp2_benchmark2/calibration_tables.csv', 'wp2_benchmark2/design_table.csv'], 'WP2_benchmark2_figure_draft.png/.pdf', 1, False),
    ('wp4', 'wp4_islets', ['run_datasets.py'], ['wp4_islets/wp4_core.py', 'wp4_islets/frozen_beta_cells.csv', 'wp4_islets/ensembl_ids_17genes.csv', 'DATA:v2_data/wp4/parsed/gse81608/beta.npz', 'DATA:v2_data/wp4/parsed/gse86469/beta.npz',
                                               'DATA:v2_data/wp4/parsed/gse81547/beta.npz', 'DATA:v2_data/wp4/parsed/emtab5061/beta.npz', 'DATA:gse84133/parsed/human1_beta_counts.npz'], 'results/ (cell-level analysis, unadjusted and depth-adjusted)', 2, False),
    ('wp4', 'wp4_islets', ['meta_analysis.py'], ['wp4_islets/results/index.csv'], 'meta/ (random-effects pooling, decision rule, pair-level rho)', 1, False),
    ('wp4', 'wp4_islets', ['gene_graph_null.py'], ['wp4_islets/results/index.csv', 'wp4_islets/meta/meta_summary.json'], 'meta/gene_graph_* (post hoc add-on)', 3, False),
    ('wp4', 'wp4_islets', ['metacell.py'], ['wp4_islets/results/index.csv'], 'metacell/ (20 partitions x 3 sizes x 6 datasets)', 25, True),
    ('wp4', 'wp4_islets', ['perturbation_gse60803.py'], ['wp4_islets/perturbation/candidates_decisions.csv', 'wp5_reporting/effect_sizes.csv', 'DATA:v2_data/wp4/GSE60803/GSE60803_series_matrix.txt.gz', 'DATA:v2_data/wp4/GSE60803/GPL11532.annot.gz'], 'perturbation/gse60803_*', 1, False),
    ('wp4', 'wp4_islets', ['positive_control.py'], ['wp4_islets/results/index.csv', 'wp4_islets/rp_ids_human.csv', 'wp4_islets/meta/meta_summary.json', 'DATA:v2_data/wp4/parsed/gse81608/beta.npz', 'DATA:v2_data/wp4/parsed/gse86469/beta.npz', 'DATA:v2_data/wp4/parsed/gse81547/beta.npz',
                                                    'DATA:v2_data/wp4/parsed/emtab5061/beta.npz', 'DATA:gse84133/parsed/human1_beta_counts.npz'], 'positive_control/ (spike-in and ribosomal control, extension A1-A3)', 15, True),
    ('wp4', 'wp4_islets', ['wider_search_analyze.py'], ['wp4_islets/wider_search/frozen_candidates.csv', 'wp4_islets/perturbation/gse60803_fold_changes.csv', 'wp5_reporting/effect_sizes.csv', 'DATA:v2_data/wp4/wider/GSE167073/GSE167073_series_matrix.txt.gz',
                                                       'DATA:v2_data/wp4/wider/GSE167073/GPL17400_data.txt'], 'wider_search/ concordance of the kept perturbation series (extension B4-B5; the live searches wider_search.py, wider_search_meta.py and wider_search_freeze.py are NOT re-run: their results are frozen and fingerprinted)', 2, False),
    ('wp4', 'wp4_islets', ['make_tables.py'], ['wp4_islets/results/index.csv', 'wp4_islets/meta/meta_summary.json', 'wp4_islets/perturbation/gse60803_concordance.json'], 'WP4_summary_tables.md', 1, False),
    ('wp4', 'wp4_islets', ['make_fig7_v2.py'], ['wp4_islets/results/index.csv', 'wp4_islets/meta/meta_summary.json', 'wp4_islets/perturbation/gse60803_fold_changes.csv'], 'figures/Figure_7_v2_islet_crosscheck, Figure_8_draft_islet_meta', 1, False),
]
HEAVY_NOTE = {'wp2 calibration run': 'python -I wp2_benchmark2/calibrate_t2.py run 2 (25 min); design: python -I wp2_benchmark2/design_table.py run 2 300 (10 min); held-out arrays and RNA-seq need the NCBI downloads, see README_v2.md',
              'wp4 downloads and parsing': 'wp4_islets/fetch_expression.py GSE81547 and E-MTAB-5061 (176 MB and 159 MB), fetch_ensembl_ids.py, fetch_gse60803.py, then parse_new_datasets.py and parse_gse81608.py (see wp4_islets/README.md section 14)',
              'wp4 metacell': 'python -I wp4_islets/metacell.py (about 20 to 30 min on this machine; one heavy job at a time)'}


def check(dry):
    missing = []
    for pkg, folder, cmd, inputs, out, mins, heavy in STEPS:
        sp = os.path.join(V2, folder, cmd[0])
        if not os.path.exists(sp): missing.append((pkg, folder + '/' + cmd[0])); print('MISSING SCRIPT [%s]: %s' % (pkg, sp))
        for rel in inputs:
            if rel.startswith('DATA:'):
                if not ROOT: continue                               # not checkable without V2_DATA_ROOT (the note at the end says so)
                p = os.path.join(ROOT, *rel[5:].split('/'))
            else: p = os.path.join(V2, rel)
            if not os.path.exists(p): missing.append((pkg, rel)); print('MISSING INPUT [%s]: %s' % (pkg, p))
        print(('would run' if dry else 'running'), '[%s] %s %s  -> %s (about %d min%s)' % (pkg, folder, ' '.join(cmd), out, mins, ', heavy' if heavy else ''))
    if ROOT and not os.path.isdir(ROOT): print('NOTE: V2_DATA_ROOT does not exist:', ROOT)
    if not ROOT: print('NOTE: V2_DATA_ROOT is not set: only steps that read cached files in this folder are re-run; the download steps need it (README_v2.md)')
    return missing


def main():
    a = sys.argv[1:]; dry = '--dry-run' in a; only = None; skip_heavy = '--skip-heavy' in a
    if '--only' in a: only = set(a[a.index('--only') + 1].split(','))
    if '--list' in a: print(HEAVY_NOTE); return
    missing = check(dry)
    if dry: print('dry run finished;', 'all inputs present' if not missing else '%d inputs missing' % len(missing)); return
    if missing: raise SystemExit('stopping: %d inputs are missing (listed above)' % len(missing))
    for pkg, folder, cmd, inputs, out, mins, heavy in STEPS:
        if only and pkg not in only: continue
        if skip_heavy and heavy: print('[%s] %s skipped (heavy step; run it alone, see --list)' % (pkg, cmd[0])); continue
        t0 = time.time(); r = subprocess.run(PY + [os.path.join(V2, folder, cmd[0])] + cmd[1:], cwd=os.path.join(V2, folder))
        print('[%s] %s exit %d in %.0f s' % (pkg, cmd[0], r.returncode, time.time() - t0))
        if r.returncode: raise SystemExit('step failed: ' + ' '.join(cmd))


if __name__ == '__main__': main()
