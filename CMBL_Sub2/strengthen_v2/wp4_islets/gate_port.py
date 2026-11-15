"""gate_port.py - WP4 gate (protocol 4.1), stage 2: the new module wp4_core.py (unadjusted mode) must reproduce the stored results of the two earlier cross-checks before anything new is added to it:
GSE81608, GSE84133 human, mouse, human non-diabetic donors only (10% rule) and (5% rule). Compared: detection of the 17 genes, rho / p / q of all 136 pairs, the pair-level null of every evaluable consensus pair, and
the summary numbers (medians, empirical P of the median and the mean, counts), at tolerance 1e-9 (empirical P exactly). Appends the block to gate_output.txt.
Usage: python -I gate_port.py"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd
import wp4_core as C

HERE = os.path.dirname(os.path.abspath(__file__)); A = os.path.normpath(os.path.join(HERE, '..', '..'))
S1 = os.path.join(A, 'human_islet_GSE81608', 'results'); S2 = os.path.join(A, 'islet_GSE84133', 'results'); BARON = os.path.join(C.DATA_ROOT, 'gse84133', 'parsed')
lines = ['WP4 GATE - stage 2: new module wp4_core.py (unadjusted mode) against the stored results of the earlier scripts', 'written %s (+0700) by gate_port.py' % time.strftime('%Y-%m-%d %H:%M:%S'), '']
fails = []


def close(a, b, tol=1e-9): return np.allclose(np.asarray(a, float), np.asarray(b, float), rtol=tol, atol=tol, equal_nan=True)


def compare(label, res, pairs_csv, rho_col, null_csv, sumd, det_csv=None, det_col=None, edge_key='edges_q05_rho02'):
    P = pd.read_csv(pairs_csv); R = res['pairs']; ok = []
    ok.append(('136 pairs: gene order and evaluable flags', list(P.gene1) == list(R.gene1) and list(P.gene2) == list(R.gene2) and list(P.evaluable) == list(R.evaluable)))
    ok.append(('136 pairs: rho (max abs diff %.2g)' % np.nanmax(np.abs(P[rho_col].to_numpy(float) - R.rho.to_numpy(float))), close(P[rho_col], R.rho)))
    ok.append(('136 pairs: p', close(P.p, R.p))); ok.append(('136 pairs: q', close(P.q, R.q)))
    if det_csv is not None:
        Dd = pd.read_csv(det_csv); ok.append(('detection of the 17 genes (%)', close(Dd[det_col].to_numpy(float), res['detection'].detected_pct.to_numpy(float))))
    if null_csv is not None and os.path.exists(null_csv):
        N = pd.read_csv(null_csv); cn = res['cons']
        ok.append(('pair-level null: pairs in order', list(N.gene1) == list(cn.gene1) and list(N.gene2) == list(cn.gene2)))
        for a_, b_ in (('rho', 'rho'), ('q', 'q'), ('null_mean', 'null_mean'), ('null_95pct', 'null_95pct'), ('empirical_p_ge_observed', 'empirical_p_ge_observed')): ok.append(('pair-level null: ' + a_, close(N[a_], cn[b_])))
    s = res['summary']
    for k in ('beta_cells', 'donors', 'df_for_p', 'genes_evaluable', 'consensus_pairs_evaluable', 'other_pairs_evaluable', 'consensus_positive', 'consensus_sig_positive_q05', 'consensus_sig_negative_q05'):
        if k in sumd: ok.append(('summary ' + k, s[k] == sumd[k]))
    for k in ('median_rho_consensus', 'median_rho_other', 'mannwhitney_p_consensus_vs_other', 'null_median_rho_mean', 'null_median_rho_95pct', 'observed_median_rho', 'empirical_p_median_ge_observed', 'empirical_p_mean_ge_observed'):
        if k in sumd and sumd[k] is not None: ok.append(('summary ' + k + ' (%.6g vs %.6g)' % (s[k], sumd[k]), close(s[k], sumd[k])))
    if edge_key in sumd or 'edges_q05_rho02' in sumd: ok.append(('summary edges q<0.05 & |rho|>=0.2', s['edges_q05_rho02'] == sumd.get(edge_key, sumd.get('edges_q05_rho02'))))
    if 'genes_not_evaluable' in sumd: ok.append(('summary genes_not_evaluable', s['genes_not_evaluable'] == sumd['genes_not_evaluable']))
    lines.append(label)
    for t, f in ok:
        lines.append('   %-4s %s' % ('OK' if f else 'FAIL', t))
        if not f: fails.append(label + ': ' + t)
    lines.append('')


# GSE81608 (RPKM, Entrez ids)
c81 = C.load_parsed('gse81608'); r = C.analyse(c81, 0.10, ('unadjusted',))['unadjusted']
compare('GSE81608 (n = %d beta-cells, %d donors)' % (c81.n, c81.n_donors), r, os.path.join(S1, 'pairs_136_human_beta_coexpression.csv'), 'rho_human_beta_donor_adjusted', os.path.join(S1, 'consensus_pairs_vs_matched_null.csv'),
        json.load(open(os.path.join(S1, 'summary.json'))), os.path.join(S1, 'genes_expression_T2D_vs_ND.csv'), 'detected_pct', 'human_edges_q05_rho02'); del c81
# GSE84133
S = json.load(open(os.path.join(S2, 'summary_GSE84133.json')))
ch = C.load_baron('human', BARON); cm = C.load_baron('mouse', BARON)
for tag, cells, thr, key in (('human', ch, 0.10, 'human'), ('mouse', cm, 0.10, 'mouse'), ('human_nondiabetic', ch.subset(ch.disease == 'non-diabetic'), 0.10, 'human_nondiabetic'), ('human_nondiabetic_det5', ch.subset(ch.disease == 'non-diabetic'), 0.05, 'human_nondiabetic_det5')):
    r = C.analyse(cells, thr, ('unadjusted',))['unadjusted']
    compare('GSE84133 %s (n = %d, %d donors, threshold %.2f)' % (tag, cells.n, cells.n_donors, thr), r, os.path.join(S2, f'{tag}_pairs_136.csv'), 'rho_beta_donor_adjusted', os.path.join(S2, f'{tag}_consensus_pairs_vs_matched_null.csv'), S[key],
            os.path.join(S2, f'{tag}_genes_detection.csv'), 'detected_pct')
lines.append('STAGE 2 VERDICT: ' + ('PASS - wp4_core.py (unadjusted) reproduces every stored number' if not fails else 'FAIL - ' + '; '.join(fails)))
with open(os.path.join(HERE, 'gate_output.txt'), 'a', encoding='utf8') as f: f.write('\n' + '\n'.join(lines) + '\n')
print('\n'.join(lines)); sys.exit(1 if fails else 0)
