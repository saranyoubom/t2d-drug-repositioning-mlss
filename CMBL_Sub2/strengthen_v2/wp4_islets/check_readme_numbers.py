"""check_readme_numbers.py - every headline number of README.md is rebuilt from the result files and searched for, in the format used in the README (so a typing slip is caught). Prints one line per check and
ends with the number of failures (exit code 1 if any). No analysis; it only reads results/, meta/, metacell/, perturbation/ and README.md.   Usage: python -I check_readme_numbers.py"""
import os, sys, json
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); R = open(os.path.join(HERE, 'README.md'), encoding='utf8').read(); fails = []
def has(label, s):
    ok = s in R; print(('ok   ' if ok else 'FAIL ') + label + '   [' + s[:90] + ']')
    if not ok: fails.append(label)
def sg(x, d=3): return '%+.*f' % (d, x)
I = pd.read_csv(os.path.join(HERE, 'results', 'index.csv')); ms = json.load(open(os.path.join(HERE, 'meta', 'meta_summary.json')))
NAME = {'gse81608': 'GSE81608', 'gse84133_human': 'GSE84133 human', 'gse86469': 'GSE86469', 'gse81547': 'GSE81547', 'emtab5061': 'E-MTAB-5061', 'gse84133_mouse': 'GSE84133 mouse'}
for t, n in NAME.items():
    u = I[(I.tag == t) & (I.variant == 'primary') & (I['mode'] == 'unadjusted')].iloc[0]; a = I[(I.tag == t) & (I.variant == 'primary') & (I['mode'] == 'adjusted')].iloc[0]
    has(f'{n}: section 4 row', f'| {n} | {sg(u.median_rho_consensus)}; {sg(u.E)} | {sg(a.median_rho_consensus)}; {sg(a.E)} |')
    if t != 'gse84133_mouse': has(f'{n}: section 5 row', f'| {n} | {int(a.beta_cells)} | {int(a.consensus_pairs_evaluable)} | {sg(a.E)} | {a.SE_E:.3f} |')
    has(f'{n}: section 3/1 evaluable pairs count in Table-1 list', str(int(a.consensus_pairs_evaluable)))
V = {(v['variant'], v['mode']): v for v in ms['variants']}; e = [x for x in V[('primary', 'adjusted')]['estimators'] if x['estimator'].startswith('REML+HK (primary)')][0]
has('pooled primary', '%s [%s, %s]' % (sg(e['mu']), sg(e['lo']), sg(e['hi']))); has('tau2 and I2', 'tau2 %.4f, I2 %.0f%%' % (e['tau2'], 100 * e['I2']))
for (v, m), lab in ((('nondiabetic', 'adjusted'), 'non-diabetic adjusted'), (('primary', 'unadjusted'), 'unadjusted'), (('nondiabetic', 'unadjusted'), 'non-diabetic unadjusted')):
    x = [y for y in V[(v, m)]['estimators'] if y['estimator'].startswith('REML+HK (primary)')][0]; has(lab, '%s [%s, %s]' % (sg(x['mu']), sg(x['lo']), sg(x['hi'])))
for y in V[('primary', 'adjusted')]['estimators']:
    if y['estimator'] in ('DL+HK', 'DL, normal CI', 'fixed effect'): has('sensitivity ' + y['estimator'], '[%s, %s]' % (sg(y['lo']), sg(y['hi'])))
for y in V[('primary', 'adjusted')]['leave_one_out']:
    has('leave out ' + y['left_out'], '%s [%s, %s]' % (sg(y['mu']), sg(y['lo']), sg(y['hi'])))
P = pd.read_csv(os.path.join(HERE, 'meta', 'pair_level_pooled_rho.csv'))
for r in P.itertuples():
    if (r.gene1, r.gene2) in (('Cacna1c', 'Cacna1d'), ('Cacna1d', 'Ptbp1'), ('Cacna1d', 'Glp1r'), ('Tcf7l2', 'Itpr1')): has(f'pair {r.gene1}-{r.gene2}', '%s [%s, %s]' % (sg(r.pooled_rho) if (r.gene1, r.gene2) != ('Tcf7l2', 'Itpr1') else '%.3f' % r.pooled_rho, sg(r.lo), sg(r.hi)) if (r.gene1, r.gene2) != ('Tcf7l2', 'Itpr1') else '%.3f' % r.pooled_rho)
gg = pd.read_csv(os.path.join(HERE, 'meta', 'gene_graph_meta.csv')); g0 = gg[gg.estimator == 'REML+HK'].iloc[0]; has('gene-graph pooled', '%s [%s, %s]' % (sg(g0['mu']), sg(g0['lo']), sg(g0['hi'])))
G = pd.read_csv(os.path.join(HERE, 'meta', 'gene_graph_by_dataset.csv'))
for r in G.itertuples(): has('gene-graph SE ratio ' + r.tag, '%.2f' % r.SE_ratio)
Ms = pd.read_csv(os.path.join(HERE, 'metacell', 'metacell_summary.csv'))
for r in Ms.itertuples():
    if np.isnan(r.consensus_pairs_evaluable_med): continue
    has(f'metacell {r.tag} size {r.size}', ('%d; %d (%d); %s; %s' % (r.metacells_median, r.consensus_pairs_evaluable_med, r.cell_level_pairs, sg(r.median_rho_med), sg(r.E_med))) if r.size == 20 else ('%d; %d; %s; %s' % (r.metacells_median, r.consensus_pairs_evaluable_med, sg(r.median_rho_med), sg(r.E_med))))
pc = json.load(open(os.path.join(HERE, 'perturbation', 'gse60803_concordance.json')))
has('perturbation primary', '| 17 | 7 | %.2f | %s | %.2f |' % (pc['all_donors']['binomial_p'], sg(pc['all_donors']['spearman_rho']), pc['all_donors']['spearman_p'])); has('perturbation S114', '| 17 | 10 | %.2f | %s | %.2f |' % (pc['S114_only']['binomial_p'], sg(pc['S114_only']['spearman_rho']), pc['S114_only']['spearman_p']))
C = pd.read_csv(os.path.join(HERE, 'perturbation', 'candidates_decisions.csv')); has('candidates', '%d unique series' % len(C)); has('dropped', 'Dropped %d' % (C.decision == 'drop').sum())
print('FAILURES:', len(fails)); sys.exit(1 if fails else 0)
