"""summarize_validation.py - protocol v2 section 3.6 (WP2 deliverable F): endpoints on the validation sets, pooled and per set, next to the development series.
Reads dev_application_per_series.csv (development, 37 series) and heldout_arrays_results.jsonl (held-out arrays) and, if the RNA-seq arm was run, heldout_rnaseq_results.jsonl.
Writes validation_per_series.csv, validation_endpoints_long.csv, validation_summary.json, validation_summary.md and WP2_summary_table.csv / .md (one table of numbers: arms in columns).
Endpoints (all from endpoints_lib.endpoints): prevalence of within-group coupling (T_r2 P < 0.05, with T_p and T_eig), E1 (consensus count not above its null, plain P as in part 1), E2 on the raw scheme alone
and in the part 1 definitions, E3, the raw pooled-significant share against the geometric expectation 0.609 (2.7), ratio-scheme diagnostics (3.4). Descriptive comparisons: exact binomial P of each T_r2 prevalence
against 0.05, Fisher exact P for held-out versus development. Nothing is re-thresholded. Usage: python -I summarize_validation.py"""
import sys, os, json
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
from scipy import stats
import endpoints_lib as E


def load_jsonl(path, arm_of):
    rows = []
    if not os.path.exists(path): return rows
    for l in open(path, encoding='utf8'):
        r = json.loads(l)
        if r.get('status') != 'analysed': continue
        for tag in ('A', 'B'): rows.append(dict(gse=r['gse'], platform=r.get('platform', r.get('organism')), gene_set=tag, arm=arm_of(r), n_de=r['n_de'], **r[tag]['flat']))
    return rows


dev = pd.read_csv(os.path.join(HERE, 'dev_application_per_series.csv')); dev['arm'] = 'development arrays'
ho = pd.DataFrame(load_jsonl(os.path.join(HERE, 'heldout_arrays_results.jsonl'), lambda r: 'held-out ' + r['platform']))
rs = pd.DataFrame(load_jsonl(os.path.join(HERE, 'heldout_rnaseq_results.jsonl'), lambda r: 'RNA-seq ' + r['platform']))
ALL = pd.concat([dev, ho, rs], ignore_index=True); ALL.to_csv(os.path.join(HERE, 'validation_per_series.csv'), index=False)
ARMS = [('Development arrays (37)', lambda d: d.arm == 'development arrays')]
if len(ho):
    for pl in sorted(ho.platform.unique()): ARMS.append(('Held-out %s' % pl, (lambda p: lambda d: d.arm == 'held-out ' + p)(pl)))
    ARMS.append(('Held-out arrays pooled', lambda d: d.arm.str.startswith('held-out')))
if len(rs):
    for pl in sorted(rs.platform.unique()): ARMS.append(('RNA-seq %s' % pl, (lambda p: lambda d: d.arm == 'RNA-seq ' + p)(pl)))
    ARMS.append(('RNA-seq pooled', lambda d: d.arm.str.startswith('RNA-seq')))
if len(ho) or len(rs): ARMS.append(('All validation sets pooled', lambda d: d.arm.str.startswith(('held-out', 'RNA-seq'))))
# post hoc sensitivity (protocol section 8, WP2 entry): leave out the matrices in which at least one of the 17 genes has no within-group variation
ARMS.append(('Development arrays without degenerate series (post hoc)', lambda d: (d.arm == 'development arrays') & (d.n_zero_var_genes == 0)))
if len(ho) or len(rs): ARMS.append(('All validation sets pooled without degenerate series (post hoc)', lambda d: d.arm.str.startswith(('held-out', 'RNA-seq')) & (d.n_zero_var_genes == 0)))
S = {}; long_rows = []
for gs in ('A', 'B'):
    for label, sel in ARMS:
        d = ALL[(ALL.gene_set == gs) & sel(ALL)]
        if not len(d): continue
        e = E.endpoints(d); S['%s | set %s' % (label, gs)] = e; long_rows += E.long_table(label, gs, e)
pd.DataFrame(long_rows).to_csv(os.path.join(HERE, 'validation_endpoints_long.csv'), index=False); json.dump(S, open(os.path.join(HERE, 'validation_summary.json'), 'w'), indent=1, default=float)


def binom_p(d): return stats.binomtest(d['k'], d['n'], 0.05, alternative='greater').pvalue if d['n'] else float('nan')


def fisher(d1, d2):
    if not d1['n'] or not d2['n']: return float('nan')
    return stats.fisher_exact([[d1['k'], d1['n'] - d1['k']], [d2['k'], d2['n'] - d2['k']]])[1]


md = ['# Validation sets: endpoints of protocol 3.6 (WP2 deliverable F)', '']
for gs in ('A', 'B'):
    md += ['## Gene set %s (%s)' % (gs, '17 random differentially expressed genes' if gs == 'A' else '17 random expressed genes'), '']
    keys = [(l, S.get('%s | set %s' % (l, gs))) for l, _ in ARMS]; keys = [(l, e) for l, e in keys if e]
    md += ['| Endpoint | ' + ' | '.join('%s (n = %d)' % (l, e['n_series']) for l, e in keys) + ' |', '|---|' + '---|' * len(keys)]
    rowdefs = [('Within-group coupling: T_r2 P < 0.05 (primary)', lambda e: E.fmt_prop(e['coupling_T_r2'])), ('T_p P < 0.05', lambda e: E.fmt_prop(e['coupling_T_p'])), ('T_eig P < 0.05', lambda e: E.fmt_prop(e['coupling_T_eig'])),
               ('E1: consensus count not above null (P >= 0.05)', lambda e: E.fmt_prop(e['E1_count_not_above_null'])),
               ('E2 raw scheme alone: consensus pairs with adjusted q < 0.05, median [Q1, Q3]', lambda e: E.fmt_med(e['E2_retained_raw_scheme'], pct=True)),
               ('E2 any scheme (part 1 definition)', lambda e: E.fmt_med(e['E2_retained_any_scheme'], pct=True)), ('E2 all four schemes', lambda e: E.fmt_med(e['E2_retained_all_schemes'], pct=True)),
               ('Pairs significant in all schemes after adjustment', lambda e: E.fmt_med(e['E2_pairs_sig_all_after_adjustment'], nd=1)),
               ('E3 covariance share carried by group means', lambda e: E.fmt_med(e['E3_share_between'], pct=True)),
               ('Raw pooled-significant share among group-effect genes (expectation 0.61)', lambda e: E.fmt_med(e['raw_sig_share_group_effect'])), ('Series with share above 0.61', lambda e: E.fmt_prop(e['raw_sig_share_above_expectation']))]
    for s in E.SCHEMES:
        rowdefs.append(('Adjusted q < 0.05 count above its null (P < 0.05), %s scheme' % E.SCHEME_LABEL[s], (lambda s: lambda e: E.fmt_prop(e['adj_' + s]['share_P_lt_0_05']))(s)))
        rowdefs.append(('Median observed / null mean adjusted count, %s scheme' % E.SCHEME_LABEL[s], (lambda s: lambda e: E.fmt_med(e['adj_' + s]['median_ratio_obs_over_null'], nd=2))(s)))
    for name, f in rowdefs: md.append('| %s | ' % name + ' | '.join(f(e) for _, e in keys) + ' |')
    md += ['', 'T_r2 prevalence against the 5% expected under no coupling (exact binomial, one-sided): ' + '; '.join('%s P = %.2g' % (l, binom_p(e['coupling_T_r2'])) for l, e in keys)]
    d0 = S.get('Development arrays (37) | set %s' % gs); dh = S.get('Held-out arrays pooled | set %s' % gs)
    if d0 and dh: md.append('Held-out arrays versus development (Fisher exact, two-sided): T_r2 prevalence P = %.2g; consensus count not above null P = %.2g.' % (fisher(dh['coupling_T_r2'], d0['coupling_T_r2']), fisher(dh['E1_count_not_above_null'], d0['E1_count_not_above_null'])))
    md.append('')
open(os.path.join(HERE, 'validation_summary.md'), 'w', encoding='utf8').write('\n'.join(md)); print('\n'.join(md))
# ---- one summary table of numbers (arms in columns), gene set A and B stacked
tab = []
for gs in ('A', 'B'):
    for name, f in rowdefs:
        tab.append(dict(gene_set=gs, endpoint=name, **{l: f(S['%s | set %s' % (l, gs)]) for l, _ in ARMS if '%s | set %s' % (l, gs) in S}))
T = pd.DataFrame(tab); T.to_csv(os.path.join(HERE, 'WP2_summary_table.csv'), index=False)
cols = list(T.columns); lines = ['| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)] + ['| ' + ' | '.join(str(v) for v in r) + ' |' for r in T.itertuples(index=False)]
open(os.path.join(HERE, 'WP2_summary_table.md'), 'w', encoding='utf8').write('\n'.join(lines) + '\n')
