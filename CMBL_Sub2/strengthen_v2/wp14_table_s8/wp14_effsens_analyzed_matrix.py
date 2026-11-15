"""wp14_effsens_analyzed_matrix.py - WP14 (10 Nov 2026), WP7 finding F09: Table S8 (efficiency sensitivity at group level) recomputed from the Cq values of the ANALYZED matrix.
Why: the old Table S8 (reanalysis_17genes/effsens.py) averaged per-plate fold changes over every plate on which a gene occurs in the 2025 re-exports (the 3 June plates included, which supplied no analyzed value; the
26 August DAPT group with all four Gapdh wells). Here the Cq values are the ones the reported analysis used: ../wp1_ins1_structure/ins1_sheet_cq.csv (target_cq and reference_cq of every replicate row of the 17 analyzed genes).
Method (as in effsens.py, flag logic unchanged): per gene and treatment (DAPT, DKK-1) against CTRL: dT = mean target Cq (treated) - mean target Cq (CTRL); dG = the same for the gene's reference (Gapdh) rows;
FC(E = 100%) = 2^-(dT - dG); FC(E amplicon) = (1 + E_gene)^-dT / (1 + E_Gapdh)^-dG with the mean per-amplicon efficiencies of Additional file 1 Table S1; ratio = FC(E amplicon) / FC(E = 100%).
Flag: 'direction flips' if the two fold changes lie on different sides of 1; '>1.5x shift' if |log2 ratio| > 0.585; '>1.25x shift' if > 0.32. The reported fold change (mean of the per-replicate 2^-dCt of the treated group
divided by that of CTRL) is listed for comparison. Read only; usage: python -I wp14_effsens_analyzed_matrix.py"""
import sys, os, re, csv, math, zipfile, html
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.abspath(__file__))
STR = os.path.abspath(os.path.join(HERE, '..'))
P = os.path.abspath(os.path.join(STR, '..', '..'))
AF1 = os.path.join(P, 'Supplementary', 'Additional file 1_JTM_working.docx')
G17 = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']


def table0_rows(path):
    x = zipfile.ZipFile(path).read('word/document.xml').decode('utf8')
    tbl = re.search(r'<w:tbl>.*?</w:tbl>', x, flags=re.S).group(0)
    rows = []
    for tr in re.findall(r'<w:tr[ >].*?</w:tr>', tbl, flags=re.S):
        cells = [''.join(html.unescape(t) for t in re.findall(r'<w:t(?: [^>]*)?>([^<]*)</w:t>', tc)) for tc in re.findall(r'<w:tc>.*?</w:tc>', tr, flags=re.S)]
        rows.append(cells)
    return rows


EFF = {}
rows = table0_rows(AF1)
for c in rows[2:]:
    if len(c) < 7 or c[0] == c[1]: continue
    m = re.findall(r'\d+\.\d+|\d+', c[6])
    if m and re.match(r'[A-Za-z0-9]+', c[0]): EFF[re.match(r'[A-Za-z0-9]+', c[0]).group(0)] = float(m[0]) / 100      # gene names read 'Glp1r(in-house)' once the runs are joined: keep the leading symbol
missing = [g for g in G17 + ['Gapdh'] if g not in EFF]
cq = list(csv.DictReader(open(os.path.join(STR, 'wp1_ins1_structure', 'ins1_sheet_cq.csv'), encoding='utf8')))
by = {}
for r in cq:
    if r['gene'] in G17: by.setdefault(r['gene'], {}).setdefault(r['group'], []).append((float(r['target_cq']), float(r['reference_cq'])))
mean = lambda v: sum(v) / len(v)
out = []
for g in G17:
    rec = {'gene': g, 'E_gene': EFF.get(g), 'E_Gapdh': EFF.get('Gapdh')}
    for t, tag in (('DAPT', 'DAPT'), ('DKK-1', 'DKK1')):
        ct, tt = by[g]['CTRL'], by[g][t]
        dT = mean([a for a, b in tt]) - mean([a for a, b in ct]); dG = mean([b for a, b in tt]) - mean([b for a, b in ct])
        fc100 = 2 ** -(dT - dG); fcE = (1 + EFF[g]) ** (-dT) / (1 + EFF['Gapdh']) ** (-dG)
        rep = mean([2 ** -(a - b) for a, b in tt]) / mean([2 ** -(a - b) for a, b in ct])
        rec[tag + '_FC_reported'] = rep; rec[tag + '_FC_E100'] = fc100; rec[tag + '_FC_Eamplicon'] = fcE; rec[tag + '_ratio'] = fcE / fc100
    flags = []
    for tag in ('DAPT', 'DKK1'):
        A, B = rec[tag + '_FC_E100'], rec[tag + '_FC_Eamplicon']
        if (A - 1) * (B - 1) < 0: flags.append(tag + ': direction flips')
        elif abs(math.log2(B / A)) > 0.585: flags.append(tag + ': >1.5x shift')
        elif abs(math.log2(B / A)) > 0.32: flags.append(tag + ': >1.25x shift')
    rec['flag'] = '; '.join(flags); out.append(rec)
os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
with open(os.path.join(HERE, 'results', 'efficiency_sensitivity_analyzed_matrix.csv'), 'w', encoding='utf-8-sig', newline='') as f:
    w = csv.writer(f); w.writerow(['gene', 'mean_amplicon_efficiency', 'DAPT_FC_reported', 'DAPT_FC_E100', 'DAPT_FC_Eamplicon', 'DAPT_ratio', 'DKK1_FC_reported', 'DKK1_FC_E100', 'DKK1_FC_Eamplicon', 'DKK1_ratio', 'flag'])
    for r in out: w.writerow([r['gene'], round(r['E_gene'], 3), round(r['DAPT_FC_reported'], 3), round(r['DAPT_FC_E100'], 3), round(r['DAPT_FC_Eamplicon'], 3), round(r['DAPT_ratio'], 3), round(r['DKK1_FC_reported'], 3), round(r['DKK1_FC_E100'], 3), round(r['DKK1_FC_Eamplicon'], 3), round(r['DKK1_ratio'], 3), r['flag']])
print('efficiencies missing for:', missing, '| E_Gapdh', EFF.get('Gapdh'))
print('%-8s %5s | DAPT rep / E100 / Eamp / ratio | DKK-1 rep / E100 / Eamp / ratio | flag' % ('gene', 'E'))
for r in out:
    print('%-8s %5.3f | %6.3f %6.3f %6.3f %5.3f | %6.3f %6.3f %6.3f %5.3f | %s' % (r['gene'], r['E_gene'], r['DAPT_FC_reported'], r['DAPT_FC_E100'], r['DAPT_FC_Eamplicon'], r['DAPT_ratio'], r['DKK1_FC_reported'], r['DKK1_FC_E100'], r['DKK1_FC_Eamplicon'], r['DKK1_ratio'], r['flag']))
flips = [r['gene'] for r in out if 'flips' in r['flag']]
print('direction flips:', flips, '| any shift > 1.25x:', sum(1 for r in out if r['flag']), 'of 17 genes | max |log2 ratio|: %.3f' % max(abs(math.log2(r[t + '_ratio'])) for r in out for t in ('DAPT', 'DKK1')))
print('same direction as the reported fold change (E100 group-level vs reported):', sum(1 for r in out for t in ('DAPT', 'DKK1') if (r[t + '_FC_reported'] - 1) * (r[t + '_FC_E100'] - 1) > 0), 'of 34')
