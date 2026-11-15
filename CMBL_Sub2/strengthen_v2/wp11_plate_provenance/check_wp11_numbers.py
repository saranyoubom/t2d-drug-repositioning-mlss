"""check_wp11_numbers.py - WP11 (8 Nov 2026): independent check that the numbers of README.md and WORDING_provenance.md are the numbers of the result files, and that three key facts hold when recomputed from the raw exports with a separate reader.
Part 1 (transcription): every number is looked up in the saved result files (results/*.txt, *.json) and must appear in the README / wording text.
Part 2 (recomputation from raw files, own reader): 26 August DAPT Gapdh wells, well/Cq agreement of the workbook, and the equality of the plate-2 blocks with the sheet.
Read only. Usage: python -I check_wp11_numbers.py   (writes check_wp11_last_run.txt)"""
import sys; sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import os, re, io, csv, json, zipfile
import numpy as np, openpyxl
HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
RES = os.path.join(HERE, 'results')
RT = os.path.abspath(os.path.join(P, '..', '..', 'Undergrad_Glucose uptake signaling', 'RT-PCR'))
W1 = os.path.join(P, 'CMBL_Sub2', 'strengthen_v2', 'wp1_ins1_structure')
README = open(os.path.join(HERE, 'README.md'), encoding='utf8').read()
WORDING = open(os.path.join(HERE, 'WORDING_provenance.md'), encoding='utf8').read() if os.path.exists(os.path.join(HERE, 'WORDING_provenance.md')) else ''
TXT = README + '\n' + WORDING
out = []; fails = 0


def check(name, ok, detail=''):
    global fails
    out.append('%s  %s %s' % ('PASS' if ok else 'FAIL', name, detail)); fails += 0 if ok else 1


def has(s): return s in TXT


def read(f): return open(os.path.join(RES, f), encoding='utf8').read()


# ---------------------------------------------------------------- part 1: transcription
tf = read('threshold_fit.txt'); rv = read('refA_vs_plate1.txt'); ao = read('aug26_nov09_order_check.txt'); aj = json.load(open(os.path.join(RES, 'aug26_nov09_order_check.json'), encoding='utf8')); t150 = read('threshold150_summary.txt')
sj = json.load(open(os.path.join(P, 'CMBL_Sub2', 'strengthen_v2', 'wp10_versions', 'results', 'string_web_export_check.json'), encoding='utf8'))
m = re.search(r'T = ([\d.]+) RFU, mean \|recomputed - sheet\| = ([\d.]+) cycles over (\d+) values', tf)
check('fit: best T 209.91, mean 0.0398 (0.040), 92 values', m and abs(float(m.group(1)) - 209.91) < 0.006 and abs(float(m.group(2)) - 0.0398) < 0.0006 and m.group(3) == '92' and has('209.91') and has('0.040 cycles over 92 values'), m.group(0) if m else '')
check('fit: plate-1 export threshold 45.700', 'threshold of the export = 45.700' in tf and has('45.70'))
for known, fitted in [('150.000', '150.000'), ('138.061', '138.060'), ('36.072', '36.070')]:
    check('fit validation: known %s returns %s' % (known, fitted), re.search(r'known T %s\s+fitted T %s' % (re.escape(known), re.escape(fitted)), tf) is not None)
check('fit: Gapdh-only T 211.59, target-only 209.53, five genes 209.91 (0.0293)', 'Gapdh blocks T = 211.59' in tf and 'target blocks T = 209.53' in tf and 'T = 209.91 (mean 0.0293)' in tf and has('211.6') and has('209.5') and has('0.029'))
check('fit: failed well G04 of Wnt2b DAPT, peak 89.6', 'G04 (wnt2b DAPT, 89.6 RFU)' in tf and has('G04') and has('89.6'))
dm = re.search(r'sheet minus re-export at 45\.70 RFU, 15 gene-group means: mean ([+-][\d.]+), SD ([\d.]+), range ([+-][\d.]+) to ([+-][\d.]+)', rv)
check('refA: dCq re-export diff +0.22 (SD 0.11; +0.06 to +0.51)', dm and dm.groups() == ('+0.22', '0.11', '+0.06', '+0.51') and has('+0.22') and has('+0.06 to +0.51'), str(dm.groups()) if dm else '')
dt = re.search(r'sheet minus curves at 209\.91 RFU, 15 gene-group means: mean ([+-][\d.]+), SD ([\d.]+), range ([+-][\d.]+) to ([+-][\d.]+)', rv)
check('refA: dCq diff at 209.91 -0.03 (SD 0.04; -0.10 to +0.02)', dt and dt.groups() == ('-0.03', '0.04', '-0.10', '+0.02') and has('-0.03') and has('-0.10 to +0.02'), str(dt.groups()) if dt else '')
d1 = re.search(r'sheet minus re-export \(45\.70 RFU\): mean ([+-][\d.]+), SD ([\d.]+), largest absolute ([\d.]+)', rv)
d2 = re.search(r'sheet minus curves at 209\.91 RFU: mean ([+-][\d.]+), SD ([\d.]+), largest absolute ([\d.]+)', rv)
check('refA: ddCq diff re-export -0.14 (SD 0.09; 0.30)', d1 and d1.groups() == ('-0.14', '0.09', '0.30'), str(d1.groups()) if d1 else '')
check('refA: ddCq diff at 209.91 +0.007 (SD 0.025; 0.042)', d2 and d2.groups() == ('+0.007', '0.025', '0.042') and has('+0.007') and has('0.025') and has('0.042'), str(d2.groups()) if d2 else '')
sw = aj['summary_workbook']
check('26 Aug: 32 of 72 rows equal the instrument value of the named well (target 26 of 60, Gapdh 6 of 12)', sw['cq_equal_to_instrument_value_of_named_well'] == 32 and sw['rows'] == 72 and sw['target_rows'] == [26, 60] and sw['gapdh_rows'] == [6, 12] and has('32 of 72') and has('26 of 60') and has('6 of 12'))
check('26 Aug: 40 of 72 rows no longer match', 72 - sw['cq_equal_to_instrument_value_of_named_well'] == 40 and has('40 of 72'))
check('26 Aug: 17 of 18 blocks ascending in the workbook; plate order ascending 1 of 15', '17 of 18' in ao and aj['matrix_aug26']['plate_order_ascending_instrument'] == [1, 15] and has('17 of 18') and has('1 of the 15'))
check('26 Aug: matrix target blocks ascending 15 of 15 and same set 15 of 15', aj['matrix_aug26']['target_blocks_ascending'] == [15, 15] and aj['matrix_aug26']['target_blocks_same_set_as_instrument'] == [15, 15])
check('26 Aug: only the Gapdh DAPT block differs in value set', sw['blocks_different'] == ['Gapdh DAPT'] and sw['blocks_same_value_set_as_instrument'] == 17 and sw['blocks_total'] == 18)
pd = aj['pairing_descriptive']
check('pairing table: as analyzed 0.132 / +0.536 / +0.662', pd['sheet'] == {'within_group_sd_of_dCq': 0.132, 'mean_residual_r_between_genes': 0.536, 'mean_DAPT_minus_CTRL_dCq': 0.662} and has('| 0.132 | +0.536 | +0.662 |'))
check('pairing table: rank instrument 0.397 / +0.975 / -0.276', pd['rank_instrument'] == {'within_group_sd_of_dCq': 0.397, 'mean_residual_r_between_genes': 0.975, 'mean_DAPT_minus_CTRL_dCq': -0.276} and has('| 0.397 | +0.975 | -0.276 |'))
check('pairing table: plate-order instrument 0.481 / +0.938 / -0.276', pd['plate_instrument'] == {'within_group_sd_of_dCq': 0.481, 'mean_residual_r_between_genes': 0.938, 'mean_DAPT_minus_CTRL_dCq': -0.276} and has('| 0.481 | +0.938 | -0.276 |'))
check('9 Nov: 11 blocks, same sequence 11, ascending 0', aj['nov09'] == {'blocks': 11, 'same_sequence_as_export': 11, 'same_set_as_export': 11, 'ascending_in_sheet': 0, 'ascending_in_export': 0} and has('11 gene-by-group blocks'))
check('threshold-150 shifts: Gr1 +2.23, Gr2 +1.97, Round2 +0.16', 'Gr1     n 78  median +2.23' in t150 and 'Gr2     n 78  median +1.97' in t150 and 'Round2  n 72  median +0.16' in t150 and has('+0.16'))
before = t150.split('=== exports with the threshold as found (Before) ===')[1].split('=== shift')[0]
ok5 = all(re.search(r'%s\s+target Cq blocks equal to the export \(CTRL / DAPT / DKK-1\): Round2 / Round2 / Round2' % g, before) for g in ('Hes1', 'Hey1', 'Glut2', 'Ins1', 'Glp1r'))
check('Before: Hes1, Hey1, Glut2, Ins1, Glp1r target and Gapdh blocks = Round2 (plate 2) in all three groups', ok5 and all(re.search(r'%s\s+reference blocks equal to Gapdh of: Round2 / Round2 / Round2' % g, before) for g in ('Hes1', 'Hey1', 'Glut2', 'Ins1', 'Glp1r')))
nomatch = all('- / - / -' in line for line in t150.splitlines() if 'target Cq blocks equal' in line and any(line.strip().startswith(g) for g in ('Wnt2 ', 'Wnt2b', 'Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7')))
check('Wnt-family genes match no block of Gr1, Gr2 or plate 2 at either threshold', nomatch and ('Hes1   target Cq blocks equal to the export (CTRL / DAPT / DKK-1): - / - / -' in t150))
check('STRING: web export = API 12.5 (59 pairs, 0 score difference, same classification)', sj['distinct_pairs']['web_12_5'] == 59 and sj['web_vs_api_12_5']['max_abs_score_difference_in_common_pairs'] == 0.0 and sj['classification_web_equals_api'] and sj['classification_from_web_export_12_5']['alias'] == {'Known Validated': 5, 'Novel Discovery': 16, 'STRING Only': 34, 'Not Validated': 65, 'in_STRING': 39, 'coverage_percent': 32.5} and has('5 / 16 / 34 / 65'))
c = sj['web_12_5_vs_old_12_0']
check('STRING: old 12.0 export, 54 of 55 common pairs changed score, 4 gained, 4 lost', c['pairs_in_both'] == 55 and c['score_change_in_common_pairs']['n_changed'] == 54 and len(c['pairs_gained_in_12_5']) == 4 and len(c['pairs_lost_in_12_5']) == 4 and has('54 of the 55') and has('4 pairs are gained and 4 lost'))
# ---------------------------------------------------------------- part 2: recomputation from the raw files with a separate reader
CAN = {'[content_types].xml': '[Content_Types].xml', 'xl/sharedstrings.xml': 'xl/sharedStrings.xml', 'docprops/core.xml': 'docProps/core.xml', 'docprops/app.xml': 'docProps/app.xml'}


def sheet0(path):
    z = zipfile.ZipFile(path); b = io.BytesIO(); zo = zipfile.ZipFile(b, 'w', zipfile.ZIP_DEFLATED)
    for n in z.namelist():
        l = n.lower()
        nn = CAN.get(l) or ('xl/worksheets/_rels/' + n.split('/')[-1] if l.startswith('xl/worksheets/_rels/') else ('xl/worksheets/' + n.split('/')[-1] if l.startswith('xl/worksheets/') and l.endswith('.xml') else n))
        zo.writestr(nn, z.read(n))
    zo.close(); b.seek(0); wb = openpyxl.load_workbook(b, data_only=True); return wb['0'] if '0' in wb.sheetnames else wb[wb.sheetnames[0]]


ex = sheet0(os.path.join(RT, '2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Cq Results.xlsx'))
rows = list(ex.iter_rows(values_only=True)); h = [str(x) for x in rows[0]]; iw, it, ic = h.index('Well'), h.index('Target'), h.index('Cq')
inst = {}
for r in rows[1:]:
    if r[iw] and isinstance(r[ic], (int, float)):
        t = str(r[it]).split('_')[0].lower(); inst[(str(r[iw]), t)] = float(r[ic])
check('raw: 26 Aug Gapdh DAPT wells A10 24.987, B10 24.132, C10 22.569, D10 22.798', all(abs(inst[(w, 'gapdh')] - v) < 0.0006 for w, v in (('A10', 24.987), ('B10', 24.132), ('C10', 22.569), ('D10', 22.798))), str([round(inst[(w, 'gapdh')], 3) for w in ('A10', 'B10', 'C10', 'D10')]))
sm = sheet0(os.path.join(RT, '2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Summary.xlsx'))
n = ok = 0; okt = nt = 0
for r in sm.iter_rows(min_row=2, values_only=True):
    w, tg, cq = r[1], r[3], r[6]
    if w is None or not isinstance(cq, (int, float)): continue
    n += 1; same = abs(inst[(str(w), str(tg).lower())] - cq) < 0.0051; ok += int(same)
    if str(tg).lower() != 'gapdh': nt += 1; okt += int(same)
check('raw: workbook Cq equals the instrument Cq of the named well in 32 of 72 rows (target 26 of 60)', (n, ok, nt, okt) == (72, 32, 60, 26), str((n, ok, nt, okt)))
dap = [r[6] for r in sm.iter_rows(min_row=2, values_only=True) if r[3] == 'Gapdh' and r[5] == 'DAPT']
check('raw: workbook Gapdh DAPT cells 22.569, 22.798, 22.569, 22.798', [round(x, 3) for x in dap] == [22.569, 22.798, 22.569, 22.798], str([round(x, 3) for x in dap]))
# plate 2 'Before' against the sheet, recomputed here
b = sheet0(os.path.join(P, 'Exported_data', 'Before', '2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 2 -  Quantification Cq Results.xlsx'))
rr = list(b.iter_rows(values_only=True)); hh = [str(x) for x in rr[0]]; iw, it, isamp, ic = hh.index('Well'), hh.index('Target'), hh.index('Sample'), hh.index('Cq')
blk = {}
for r in rr[1:]:
    if r[it] is None or not isinstance(r[ic], (int, float)): continue
    t, s = str(r[it]), r[isamp]
    if '_' in t and not s: t, s = t.rsplit('_', 1)
    g = re.sub(r'[^a-z0-9]', '', re.sub(r'^\d+\.\s*', '', t).lower()); k = 'CTRL' if 'ctrl' in str(s).lower() else 'DAPT' if 'dapt' in str(s).lower() else 'DKK-1'
    blk.setdefault((g, k), []).append(float(r[ic]))
sh = {}
for r in csv.DictReader(open(os.path.join(W1, 'ins1_sheet_cq.csv'), encoding='utf8')):
    sh.setdefault((re.sub(r'[^a-z0-9]', '', r['gene'].lower()), r['group']), []).append((float(r['target_cq']), float(r['reference_cq'])))
eq = 0; eqg = 0
for g in ('hes1', 'hey1', 'glut2', 'ins1', 'glp1r'):
    for k in ('CTRL', 'DAPT', 'DKK-1'):
        eq += int(np.allclose(np.sort([x[0] for x in sh[(g, k)]]), np.sort(blk[(g, k)]), atol=0.0051))
        eqg += int(np.allclose(np.sort([x[1] for x in sh[(g, k)]]), np.sort(blk[('gapdh', k)]), atol=0.0051))
check('raw: plate-2 export (threshold as found) equals the sheet in 15 of 15 target blocks and 15 of 15 Gapdh vectors', (eq, eqg) == (15, 15), str((eq, eqg)))
tabs = list(csv.DictReader(open(os.path.join(RES, 'aug26_workbook_vs_instrument.csv'), encoding='utf8')))
check('export_tables: aug26 csv has 72 rows, 32 with the Cq of the named well, 4 rows outside a matching value set (Gapdh DAPT)', len(tabs) == 72 and sum(int(r['workbook_Cq_equals_instrument_Cq_of_named_well']) for r in tabs) == 32 and sum(1 for r in tabs if r['block_value_set_equals_instrument'] == '0') == 4 and all(r['target'] == 'Gapdh' and r['group_in_workbook'] == 'DAPT' for r in tabs if r['block_value_set_equals_instrument'] == '0'))
vt = list(csv.DictReader(open(os.path.join(RES, 'threshold_validation.csv'), encoding='utf8'))); dq = list(csv.DictReader(open(os.path.join(RES, 'refA_dcq_comparison.csv'), encoding='utf8')))
check('export_tables: 4 validation rows, 15 dCq rows (mean sheet-minus-curves -0.03), 24 residual rows', len(vt) == 4 and len(dq) == 15 and round(sum(float(r['dCq_sheet_minus_curves']) for r in dq) / 15, 2) == -0.03 and len(list(csv.DictReader(open(os.path.join(RES, 'refA_block_residuals.csv'), encoding='utf8')))) == 24)
out.append('%d checks, %d failures' % (len(out), fails))
open(os.path.join(HERE, 'check_wp11_last_run.txt'), 'w', encoding='utf8').write('\n'.join(out) + '\n')
print('\n'.join(out))
