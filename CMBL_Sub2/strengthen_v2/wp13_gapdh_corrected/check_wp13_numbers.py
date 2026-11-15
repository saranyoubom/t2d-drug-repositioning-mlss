"""check_wp13_numbers.py - WP13 (9 Nov 2026): independent check of README.md. Own code: fold changes from the Cq values, categories from the consensus lists and the STRING export, plateau values from the well table.
It reads results/gapdh_versions.json (the stored output of wp13_gapdh_versions.py) and the raw inputs; it does not import the pipeline. Usage: python -I check_wp13_numbers.py (writes check_wp13_last_run.txt)"""
import sys, os, json, csv, re
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.abspath(__file__)); STR = os.path.abspath(os.path.join(HERE, '..'))
RE = os.path.abspath(os.path.join(HERE, '..', '..', 'reanalysis_17genes')); R17 = os.path.abspath(os.path.join(HERE, '..', '..', 'R_17genes'))
TXT = open(os.path.join(HERE, 'README.md'), encoding='utf8').read()
J = json.load(open(os.path.join(HERE, 'results', 'gapdh_versions.json'), encoding='utf8'))
out, fails = [], [0]


def check(name, ok, detail=''):
    out.append('%s  %s %s' % ('PASS' if ok else 'FAIL', name, detail)); fails[0] += 0 if ok else 1


def has(s): return s in TXT


# 1. fold changes from the raw Cq values of the sheet (ins1_sheet_cq.csv): DAPT replicate k uses reference HK_k
cq = list(csv.DictReader(open(os.path.join(STR, 'wp1_ins1_structure', 'ins1_sheet_cq.csv'), encoding='utf8')))
G5 = ['Tcf7l2', 'Kcnj11', 'Cacna1d', 'Cacna1c', 'Ptbp1']
RAW = {'A10': 24.9873996418614, 'B10': 24.1318270358876, 'C10': 22.5685093010354, 'D10': 22.7982726003435}
HK = [RAW['C10'], RAW['D10'], RAW['C10'], RAW['D10']]
refs = {'reported': HK, 'C0': [(RAW['C10'] + RAW['D10']) / 2] * 4, 'C1': [sum(RAW.values()) / 4] * 4}


def fc(g, ref):
    rows = sorted([r for r in cq if r['gene'] == g], key=lambda r: int(r['replicate']))
    ctrl = [2 ** -(float(r['target_cq']) - float(r['reference_cq'])) for r in rows[0:4]]
    dapt = [2 ** -(float(r['target_cq']) - ref[k]) for k, r in enumerate(rows[4:8])]
    return (sum(dapt) / 4) / (sum(ctrl) / 4)


for v in ('reported', 'C0', 'C1'):
    mine = {g: fc(g, refs[v]) for g in G5}
    check('DAPT fold changes %s recomputed from the Cq values' % v, all(abs(mine[g] - J[v]['fc_dapt'][g]) < 1e-9 for g in G5), str({g: round(mine[g], 3) for g in G5}))
check('fold changes quoted: C0 0.78 / 0.40 / 0.76 / 0.72 / 0.62; C1 1.49 / 0.76 / 1.46 / 1.37 / 1.18', has('0.78, 0.40, 0.76, 0.72 and 0.62') and has('1.49, 0.76, 1.46, 1.37 and 1.18'))
# 2. gate and versions from the stored json
r, c0, c1 = J['reported'], J['C0'], J['C1']
check('gate: reported version = 21 consensus pairs, pooled 66/80/84/120, adjusted 0/51/11/15, hubs Glp1r 6 and 5/5/5', r['n_consensus'] == 21 and list(r['sig'].values()) == [66, 80, 84, 120] and list(r['adj_q05'].values()) == [0, 51, 11, 15] and r['degree']['Glp1r'] == 6 and sorted(r['hubs']) == ['Cacna1d', 'Glp1r', 'Ptbp1', 'Tcf7l2'])
check('gate: Fig. 5 with the alias reproduces 5 / 16 / 33 / 66 and 38 in STRING', (r['fig5_alias']['Known Validated'], r['fig5_alias']['Novel Discovery'], r['fig5_alias']['STRING Only'], r['fig5_alias']['Not Validated'], r['in_STRING_alias']) == (5, 16, 33, 66, 38))
check('C0: 24 consensus pairs, 20 of the 21 kept, lost Itpr1-Tcf7l2, 4 new; hubs add Wnt5a (degree 5)', c0['n_consensus'] == 24 and c0['kept'] == 20 and c0['lost'] == ['Itpr1-Tcf7l2'] and len(c0['new']) == 4 and sorted(c0['hubs']) == ['Cacna1d', 'Glp1r', 'Ptbp1', 'Tcf7l2', 'Wnt5a'] and has('24 consensus pairs') and has('Itpr1-Tcf7l2'))
check('C0: pooled 63/80/86/120, adjusted 0/47/36/12, r(Cacna1d, Glp1r) 0.976', list(c0['sig'].values()) == [63, 80, 86, 120] and list(c0['adj_q05'].values()) == [0, 47, 36, 12] and round(c0['r_cacna1d_glp1r'], 3) == 0.976 and has('0/47/36/12'))
check('C1: 21 consensus pairs, only 11 kept (10 lost, 10 new); hubs Glp1r 5 and Cacna1d 5 only; reference = mean of the four wells = 23.6215 (README says so twice, never 23.6225)', c1['n_consensus'] == 21 and c1['kept'] == 11 and len(c1['lost']) == 10 and len(c1['new']) == 10 and sorted(c1['hubs']) == ['Cacna1d', 'Glp1r'] and has('11 of the 21')
      and abs(sum(RAW.values()) / 4 - 23.6215) < 5e-5 and abs(J['ref']['C1'] - 23.6215) < 5e-5 and TXT.count('23.6215') >= 2 and '(mean 23.6225)' not in TXT and '23.6225 for all four' not in TXT)
# 3. independent Fig. 5 categories: consensus list x STRING export (alias) over the 120 pairs of the 16 non-Hes1 genes
edges = {frozenset(('Glut2' if a == 'Slc2a2' else a, 'Glut2' if b == 'Slc2a2' else b)) for a, b in ((x['node1'], x['node2']) for x in csv.DictReader(open(os.path.join(R17, 'string_interactions.csv'), encoding='utf-8-sig', newline='')))}
G17 = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']
allp = {frozenset((a, b)) for i, a in enumerate(G17) for b in G17[i + 1:] if 'Hes1' not in (a, b)}
check('120 pairs among the 16 genes other than Hes1', len(allp) == 120)
for v in ('reported', 'C0', 'C1'):
    cons = {frozenset(p) for p in J[v]['consensus']}
    kv = len(cons & edges & allp); nv = len(cons - edges); so = len((edges & allp) - cons); nt = len(allp - edges - cons)
    f = J[v]['fig5_alias']
    check('Fig. 5 categories %s recomputed: %d / %d / %d / %d' % (v, kv, nv, so, nt), (kv, nv, so, nt) == (f['Known Validated'], f['Novel Discovery'], f['STRING Only'], f['Not Validated']))
check('Fig. 5 with C0: 6 / 18 / 32 / 64 (38 in STRING)', has('6 / 18 / 32 / 64'))
# 4. plateau values
pl = {w['well']: w['plateau_RFU'] for w in J['qc_gapdh_wells']}
others = sorted(v for k, v in pl.items() if k not in ('A10', 'B10'))
check('plateau RFU: A10 513, B10 725; the other ten Gapdh wells 1279 to 1507', pl['A10'] == 513 and pl['B10'] == 725 and others[0] == 1279 and others[-1] == 1507 and len(others) == 10 and has('513 and 725') and has('1279 to 1507'))
med = sorted(others)[4:6]; medv = sum(med) / 2
check('A10 and B10 plateau as a share of the median of the other ten wells: 38% and 54%', round(100 * 513 / medv) == 38 and round(100 * 725 / medv) == 54 and has('38% and 54%'), 'median %.0f' % medv)
sd = lambda x: (sum((a - sum(x) / len(x)) ** 2 for a in x) / (len(x) - 1)) ** 0.5
check('Gapdh replicate SD of Cq in the DAPT group with all four wells 1.14; CTRL 0.10; DKK-1 0.15', abs(sd(list(RAW.values())) - 1.14) < 0.005 and abs(sd([23.626, 23.806, 23.659, 23.814]) - 0.10) < 0.01 and abs(sd([22.519, 22.856, 22.763, 22.819]) - 0.15) < 0.01 and has('1.14'))
out.append('%d checks, %d failures' % (len(out), fails[0]))
open(os.path.join(HERE, 'check_wp13_last_run.txt'), 'w', encoding='utf8').write('\n'.join(out) + '\n')
print('\n'.join(out))
