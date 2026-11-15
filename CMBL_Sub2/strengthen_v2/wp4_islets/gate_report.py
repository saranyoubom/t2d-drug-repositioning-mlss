"""gate_report.py - WP4 gate (protocol 4.1), stage 1: write gate_output.txt from (a) the file-level comparisons of the re-runs of the existing scripts with the stored results (logs/gate_compare_*.txt, made by
gate_compare.py) and (b) the figures quoted in the READMEs of the two earlier cross-checks and in the WP4 brief, compared with the re-run summaries at the precision at which they were quoted.
Stage 2 (gate_port.py) appends the check of the new module wp4_core.py. Usage: python -I gate_report.py"""
import os, sys, json, time
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__)); G = os.path.join(HERE, 'gate')
s81 = json.load(open(os.path.join(G, 'results_GSE81608', 'summary.json')))
s84 = json.load(open(os.path.join(G, 'results_GSE84133', 'summary_GSE84133.json')))
checks = []


def chk(label, got, rep, nd=None):
    """compare got with the reported value rep; nd = number of decimals at which rep was quoted (None = exact)"""
    ok = (got == rep) if nd is None else (got is not None and abs(round(got, nd) - rep) < 1e-12)
    checks.append((label, got, rep, ok))


a = s81
chk('GSE81608 beta-cells', a['beta_cells'], 503); chk('GSE81608 donors', a['donors'], 18); chk('GSE81608 evaluable genes (of 17)', a['genes_evaluable'], 11)
chk('GSE81608 evaluable consensus pairs (of 21)', a['consensus_pairs_evaluable'], 11); chk('GSE81608 evaluable other pairs (of 115)', a['other_pairs_evaluable'], 44)
chk('GSE81608 median rho, consensus pairs', a['median_rho_consensus'], -0.021, 3); chk('GSE81608 median rho, other pairs', a['median_rho_other'], 0.018, 3)
chk('GSE81608 Mann-Whitney P consensus vs other', a['mannwhitney_p_consensus_vs_other'], 0.16, 2); chk('GSE81608 matched random pairs, mean of the null medians', a['null_median_rho_mean'], 0.020, 3)
chk('GSE81608 95th percentile of the null median', a['null_median_rho_95pct'], 0.049, 3); chk('GSE81608 empirical P, median >= observed', a['empirical_p_median_ge_observed'], 0.985, 3)
chk('GSE81608 consensus pairs with positive rho', a['consensus_positive'], 5); chk('GSE81608 significantly positive (q < 0.05)', a['consensus_sig_positive_q05'], 0); chk('GSE81608 significantly negative (q < 0.05)', a['consensus_sig_negative_q05'], 1)
chk('GSE81608 edges q < 0.05 and abs(rho) >= 0.2', a['human_edges_q05_rho02'], 0); chk('GSE81608 Cacna1d-Glp1r rho', a['Cacna1d_Glp1r']['rho'], 0.002, 3); chk('GSE81608 Cacna1d-Glp1r P', a['Cacna1d_Glp1r']['p'], 0.97, 2)
h = s84['human']
chk('GSE84133 human beta-cells', h['beta_cells'], 2525); chk('GSE84133 human donors', h['donors'], 4); chk('GSE84133 human evaluable genes', h['genes_evaluable'], 6); chk('GSE84133 human evaluable consensus pairs', h['consensus_pairs_evaluable'], 3)
chk('GSE84133 human median rho, consensus', h['median_rho_consensus'], 0.24, 2); chk('GSE84133 human matched random median (mean of null medians)', h['null_median_rho_mean'], 0.05, 2); chk('GSE84133 human 95th percentile of the null median', h['null_median_rho_95pct'], 0.25, 2)
chk('GSE84133 human empirical P (median)', h['empirical_p_median_ge_observed'], 0.071, 3); chk('GSE84133 human empirical P (mean)', h['empirical_p_mean_ge_observed'], 0.014, 3)
chk('GSE84133 human significantly positive / negative', (h['consensus_sig_positive_q05'], h['consensus_sig_negative_q05']), (3, 0)); chk('GSE84133 human edges q < 0.05 and abs(rho) >= 0.2', h['edges_q05_rho02'], 4)
chk('GSE84133 human Cacna1d-Glp1r evaluable', h['Cacna1d_Glp1r']['evaluable'], False)
m = s84['mouse']
chk('GSE84133 mouse beta-cells', m['beta_cells'], 894); chk('GSE84133 mouse evaluable genes', m['genes_evaluable'], 5); chk('GSE84133 mouse evaluable consensus pairs', m['consensus_pairs_evaluable'], 5)
chk('GSE84133 mouse median rho, consensus', m['median_rho_consensus'], -0.04, 2); chk('GSE84133 mouse matched random median', m['null_median_rho_mean'], 0.08, 2); chk('GSE84133 mouse 95th percentile of the null median', m['null_median_rho_95pct'], 0.20, 2)
chk('GSE84133 mouse empirical P (median)', m['empirical_p_median_ge_observed'], 0.805, 3); chk('GSE84133 mouse empirical P (mean)', m['empirical_p_mean_ge_observed'], 0.521, 3)
chk('GSE84133 mouse significantly positive / negative', (m['consensus_sig_positive_q05'], m['consensus_sig_negative_q05']), (2, 2)); chk('GSE84133 mouse edges', m['edges_q05_rho02'], 1)
n_ = s84['human_nondiabetic']; chk('GSE84133 human, non-diabetic donors only, 10%: beta-cells', n_['beta_cells'], 2030); chk('  evaluable genes', n_['genes_evaluable'], 4); chk('  evaluable consensus pairs', n_['consensus_pairs_evaluable'], 0)
d5 = s84['human_nondiabetic_det5']; chk('GSE84133 human, non-diabetic donors only, 5%: evaluable genes', d5['genes_evaluable'], 8); chk('  evaluable consensus pairs', d5['consensus_pairs_evaluable'], 7)
chk('  median rho, consensus', d5['median_rho_consensus'], 0.28, 2); chk('  matched random median', d5['null_median_rho_mean'], 0.08, 2); chk('  95th percentile of the null median', d5['null_median_rho_95pct'], 0.40, 2)
chk('  empirical P (median)', d5['empirical_p_median_ge_observed'], 0.090, 3); chk('  empirical P (mean)', d5['empirical_p_mean_ge_observed'], 0.093, 3)
allok = all(c[3] for c in checks)
lines = ['WP4 GATE (protocol 4.1) - stage 1: the existing scripts re-run on the cached data', 'written %s (+0700) by gate_report.py' % time.strftime('%Y-%m-%d %H:%M:%S'), '',
         'A. Files: every output file of the two earlier cross-checks re-created and compared with the stored result (logs/gate_compare_*.txt)', '']
for f in ('gate_compare_GSE84133.txt', 'gate_compare_GSE81608.txt'): lines += open(os.path.join(HERE, 'logs', f), encoding='utf8').read().splitlines() + ['']
lines += ['  GSE84133: the existing script islet_check_baron.py was run UNCHANGED on the cached sparse beta-cell files (gate/results_GSE84133).',
          '  GSE81608: the existing script human_islet_check.py reads the whole 39851 x 1600 RPKM table as a dense float64 frame (about 510 MB, three to four copies at peak), which the 8 GB machine with 1-2 GB free should not do;',
          '  gate_human_islet_check_lean.py is a COPY in which only the loading block differs (same pandas C parser, chunks of 2500 rows, only the 503 beta-cell columns kept); every statistic line is verbatim.', '',
          'B. Numbers quoted in the READMEs and in the brief, against the re-run summaries (precision as quoted)', '']
for lab, got, rep, ok in checks: lines.append('  %-4s %-72s re-run %-24s reported %s' % ('OK' if ok else 'FAIL', lab, (('%.6g' % got) if isinstance(got, float) else str(got)), rep))
lines += ['', 'STAGE 1 VERDICT: ' + ('PASS (all %d quoted numbers and all files reproduced)' % len(checks) if allok else 'FAIL - see the FAIL lines'), '']
open(os.path.join(HERE, 'gate_output.txt'), 'w', encoding='utf8').write('\n'.join(lines))
print('\n'.join(lines[-(len(checks) + 6):])); sys.exit(0 if allok else 1)
