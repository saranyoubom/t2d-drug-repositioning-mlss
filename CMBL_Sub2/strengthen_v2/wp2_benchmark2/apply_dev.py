"""apply_dev.py - protocol v2 sections 3.3 and 3.4 on the 37 development series (WP2 deliverable C).
Reads compact/dev_core_results.jsonl (bench2_core.analyse_matrix on the 74 matrices of gate_part1.py: gene sets A and B; coupling tests 2000 permutations, consensus null 2000, calibrated adjusted counts 500).
Writes dev_application_per_series.csv, dev_endpoints_long.csv, dev_application_summary.json and dev_application_summary.md: prevalence of within-group coupling (P < 0.05) with Wilson 95% CIs
for T_r2 (primary), T_p, T_eig; the consensus-count endpoint E1; ratio-scheme diagnostics per scheme (share with P < 0.05, median observed / null); E2 recomputed on the raw scheme alone
(share of consensus pairs with raw-scheme treatment-adjusted q < 0.05), and E3 and the raw pooled-significant share against the geometric expectation 0.609.
Part 1 numbers are shown beside the recomputed ones (per_dataset_results.csv). Usage: python -I apply_dev.py"""
import sys, os, json
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__)); AO = os.path.dirname(os.path.dirname(HERE)); P1 = os.path.join(AO, 'benchmark_route3'); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
import endpoints_lib as E

recs = [json.loads(l) for l in open(os.path.join(HERE, 'compact', 'dev_core_results.jsonl'), encoding='utf8')]
T = pd.DataFrame([dict(gse=r['gse'], platform=r['platform'], gene_set=r['tag'], **r['flat']) for r in recs]); T.to_csv(os.path.join(HERE, 'dev_application_per_series.csv'), index=False)
p1 = pd.read_csv(os.path.join(P1, 'per_dataset_results.csv')); S = {}; long_rows = []; md = ['# Development series: protocol 3.3 and 3.4 (WP2 deliverable C)', '']
arms = [('development, all 37 series', None), ('development, human GPL570', 'GPL570'), ('development, mouse GPL1261', 'GPL1261'),
        ('development, without series with degenerate genes (post hoc sensitivity)', 'ND')]       # ND: no gene with zero within-group variance in the 17-gene set (protocol section 8, WP2 entry)
for gs in ('A', 'B'):
    for label, pl in arms:
        d = T[(T.gene_set == gs) & ((T.n_zero_var_genes == 0) if pl == 'ND' else ((T.platform == pl) if pl else True))]
        if not len(d): continue
        e = E.endpoints(d); S['%s | set %s' % (label, gs)] = e; long_rows += E.long_table(label, gs, e)
    # part 1 values for the same series (loop-based null, 2000 permutations) for comparison
    q = p1[p1.set == gs]; S['part 1 | set %s' % gs] = dict(E1_count_not_above_null=E.prop(q.perm_p >= 0.05), retained_any=E.med(q.retained_any), retained_all4=E.med(q.retained_all4), adj_pairs_all4=E.med(q.adj_pairs_all4), share_between=E.med(q.share_between))
pd.DataFrame(long_rows).to_csv(os.path.join(HERE, 'dev_endpoints_long.csv'), index=False); json.dump(S, open(os.path.join(HERE, 'dev_application_summary.json'), 'w'), indent=1, default=float)
for gs in ('A', 'B'):
    e = S['development, all 37 series | set %s' % gs]; q = S['part 1 | set %s' % gs]
    md += ['## Gene set %s (%s)' % (gs, '17 random differentially expressed genes, FDR < 0.05' if gs == 'A' else '17 random expressed genes'), '',
           '| Endpoint | Result (37 series) |', '|---|---|',
           '| Within-group coupling, T_r2 (primary), P < 0.05 | %s |' % E.fmt_prop(e['coupling_T_r2']), '| T_p, P < 0.05 | %s |' % E.fmt_prop(e['coupling_T_p']), '| T_eig, P < 0.05 | %s |' % E.fmt_prop(e['coupling_T_eig']),
           '| Consensus count not above its permutation null (E1; P >= 0.05), this run | %s |' % E.fmt_prop(e['E1_count_not_above_null']), '| same, part 1 (loop null) | %s |' % E.fmt_prop(q['E1_count_not_above_null']),
           '| E2 raw scheme alone: consensus pairs with raw-scheme adjusted q < 0.05 (series with consensus pairs, n = %d) | median %s |' % (e['E2_n_series_with_consensus_pairs'], E.fmt_med(e['E2_retained_raw_scheme'], pct=True)),
           '| E2 any scheme (part 1 definition), this run / part 1 | %s / %s |' % (E.fmt_med(e['E2_retained_any_scheme'], pct=True), E.fmt_med(q['retained_any'], pct=True)),
           '| E2 all four schemes, this run / part 1 | %s / %s |' % (E.fmt_med(e['E2_retained_all_schemes'], pct=True), E.fmt_med(q['retained_all4'], pct=True)),
           '| Pairs significant in all schemes after adjustment, this run / part 1 | %s / %s |' % (E.fmt_med(e['E2_pairs_sig_all_after_adjustment'], nd=1), E.fmt_med(q['adj_pairs_all4'], nd=1)),
           '| E3 share of pooled covariance carried by group means | %s |' % E.fmt_med(e['E3_share_between'], pct=True),
           '| Raw pooled-significant share among genes with a group effect (expectation 0.61) | %s; above 0.61 in %s |' % (E.fmt_med(e['raw_sig_share_group_effect']), E.fmt_prop(e['raw_sig_share_above_expectation'])), '']
    md += ['Ratio-scheme diagnostics (observed adjusted q < 0.05 count against its null from 500 within-group permutations):', '', '| Scheme | Series with P < 0.05 | Median observed | Median null mean | Median observed / null mean | Median excess |', '|---|---|---|---|---|---|']
    for s in E.SCHEMES:
        a = e['adj_' + s]; md.append('| %s | %s | %s | %s | %s (defined in %d) | %s |' % (E.SCHEME_LABEL[s], E.fmt_prop(a['share_P_lt_0_05']), E.fmt_med(a['median_observed'], nd=1), E.fmt_med(a['median_null_mean'], nd=1), E.fmt_med(a['median_ratio_obs_over_null'], nd=2), a['n_ratio_defined'], E.fmt_med(a['median_excess'], nd=1)))
    md.append('')
    for label, pl in arms[1:3]:
        d = T[(T.gene_set == gs) & (T.platform == pl)]; ee = E.endpoints(d); md.append('By platform, %s: T_r2 P < 0.05 in %s; consensus count not above null in %s.' % (label.split(', ')[1], E.fmt_prop(ee['coupling_T_r2']), E.fmt_prop(ee['E1_count_not_above_null'])))
    ee = S['%s | set %s' % (arms[3][0], gs)]
    md.append('Post hoc sensitivity, %s (n = %d): T_r2 P < 0.05 in %s; T_p in %s; consensus count not above null in %s; E2 raw scheme alone %s; covariance share carried by group means %s.' % (
        arms[3][0].split(', ')[1], ee['n_series'], E.fmt_prop(ee['coupling_T_r2']), E.fmt_prop(ee['coupling_T_p']), E.fmt_prop(ee['E1_count_not_above_null']), E.fmt_med(ee['E2_retained_raw_scheme'], pct=True), E.fmt_med(ee['E3_share_between'], pct=True)))
    md.append('')
open(os.path.join(HERE, 'dev_application_summary.md'), 'w', encoding='utf8').write('\n'.join(md)); print('\n'.join(md))
