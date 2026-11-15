"""check_wp14_numbers.py - WP14 (10 Nov 2026): independent check of results/efficiency_sensitivity_analyzed_matrix.csv and of the numbers of README.md.
Own code (numpy arrays, the Cq values read straight from ../wp1_ins1_structure/ins1_sheet_cq.csv, efficiencies typed from Additional file 1 Table S1 through a regular expression of its own).
Usage: python -I check_wp14_numbers.py   (writes check_wp14_last_run.txt)"""
import sys, os, csv, re, zipfile, html
import numpy as np
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.abspath(__file__)); STR = os.path.abspath(os.path.join(HERE, '..')); P = os.path.abspath(os.path.join(STR, '..', '..'))
TXT = open(os.path.join(HERE, 'README.md'), encoding='utf8').read()
out, fails = [], [0]


def check(name, ok, detail=''):
    out.append('%s  %s %s' % ('PASS' if ok else 'FAIL', name, detail)); fails[0] += 0 if ok else 1


x = zipfile.ZipFile(os.path.join(P, 'Supplementary', 'Additional file 1_JTM_working.docx')).read('word/document.xml').decode('utf8')
t0 = re.search(r'<w:tbl>.*?</w:tbl>', x, flags=re.S).group(0)
E = {}
for tr in re.findall(r'<w:tr[ >].*?</w:tr>', t0, flags=re.S)[2:]:
    cells = [''.join(html.unescape(s) for s in re.findall(r'<w:t(?: [^>]*)?>([^<]*)</w:t>', tc)) for tc in re.findall(r'<w:tc>.*?</w:tc>', tr, flags=re.S)]
    if len(cells) >= 7:
        m = re.match(r'[A-Za-z0-9]+', cells[0]); v = re.findall(r'\d+\.\d+|\d+', cells[6])
        if m and v: E[m.group(0)] = float(v[0]) / 100
G17 = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']
check('efficiencies of the 17 genes and Gapdh found; Gapdh 70.1%; range of the 17 genes 53.9 to 74.6%', all(g in E for g in G17 + ['Gapdh']) and abs(E['Gapdh'] - 0.701) < 1e-9 and abs(min(E[g] for g in G17) - 0.539) < 1e-9 and abs(max(E[g] for g in G17) - 0.746) < 1e-9 and 'Gapdh 70.1%' in TXT and '53.9 to 74.6%' in TXT)
cq = list(csv.DictReader(open(os.path.join(STR, 'wp1_ins1_structure', 'ins1_sheet_cq.csv'), encoding='utf8')))
res = {r['gene']: r for r in csv.DictReader(open(os.path.join(HERE, 'results', 'efficiency_sensitivity_analyzed_matrix.csv'), encoding='utf-8-sig'))}
check('17 rows in the result table', len(res) == 17 and set(res) == set(G17))
bad = 0; flips = []; shifts = 0; maxlog = 0; same = 0
for g in G17:
    rows = [r for r in cq if r['gene'] == g]
    T = np.array([float(r['target_cq']) for r in sorted(rows, key=lambda r: int(r['replicate']))]); R = np.array([float(r['reference_cq']) for r in sorted(rows, key=lambda r: int(r['replicate']))])
    for tag, sl in (('DAPT', slice(4, 8)), ('DKK1', slice(8, 12))):
        dT = T[sl].mean() - T[0:4].mean(); dG = R[sl].mean() - R[0:4].mean()
        f100 = 2.0 ** (dG - dT); fE = (1 + E[g]) ** (-dT) / (1 + E['Gapdh']) ** (-dG); rep = np.mean(2.0 ** -(T[sl] - R[sl])) / np.mean(2.0 ** -(T[0:4] - R[0:4]))
        got = (float(res[g][tag + '_FC_reported']), float(res[g][tag + '_FC_E100']), float(res[g][tag + '_FC_Eamplicon']), float(res[g][tag + '_ratio']))
        if not (abs(got[0] - rep) < 6e-4 and abs(got[1] - f100) < 6e-4 and abs(got[2] - fE) < 6e-4 and abs(got[3] - fE / f100) < 6e-4): bad += 1
        if (f100 - 1) * (fE - 1) < 0: flips.append((g, tag, round(f100, 3), round(fE, 3)))
        if abs(np.log2(fE / f100)) > 0.32: shifts += 1
        maxlog = max(maxlog, abs(np.log2(fE / f100))); same += int((rep - 1) * (f100 - 1) > 0)
check('all 17 x 2 fold changes, ratios recomputed independently (tolerance 0.0006)', bad == 0, '%d mismatches' % bad)
check('exactly one direction flip: Tcf7l2 under DKK-1, 1.002 to 0.943', flips == [('Tcf7l2', 'DKK1', 1.002, 0.943)] and 'Tcf7l2 under DKK-1' in TXT and '1.002 without and 0.943' in TXT, str(flips))
check('flagged genes: 5 of 17 (4 shifts above 1.25 and 1 flip)', sum(1 for g in G17 if res[g]['flag']) == 5 and 'Five of 17 genes' in TXT)
check('largest shift 2.0-fold (Hes1 DAPT 0.137 to 0.275)', abs(2 ** maxlog - 2.012) < 0.01 and float(res['Hes1']['DAPT_FC_E100']) == 0.137 and float(res['Hes1']['DAPT_FC_Eamplicon']) == 0.275 and '0.137 to 0.275' in TXT)
check('33 of 34 gene-by-treatment pairs on the same side of 1 as the reported fold change', same == 33 and '33 of the 34' in TXT, str(same))
check('old table quoted for Cacna1d under DAPT: 1.46 against 0.76 reported', abs(float(res['Cacna1d']['DAPT_FC_reported']) - 0.763) < 6e-4 and '1.46 in the old table against 0.76' in TXT)
out.append('%d checks, %d failures' % (len(out), fails[0]))
open(os.path.join(HERE, 'check_wp14_last_run.txt'), 'w', encoding='utf8').write('\n'.join(out) + '\n')
print('\n'.join(out))
