"""export_tables.py - WP11 (9 Nov 2026): structured copies of WP11 results for Additional file 6 (so that no sheet has to be parsed from prose).
Writes, in results/:  aug26_workbook_vs_instrument.csv (72 rows), threshold_validation.csv, refA_block_residuals.csv, refA_dcq_comparison.csv.
Inputs: the 26 August instrument export and the user's 'Quantification Summary' workbook (as in check_aug26_nov09_order.py), results/threshold_fit.txt, results/refA_vs_plate1.txt.  Read only; usage: python -I export_tables.py"""
import sys; sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import os, re, csv
HERE = os.path.dirname(os.path.abspath(__file__))
exec(compile(open(os.path.join(HERE, 'check_threshold150.py'), encoding='utf8').read().split("RUNS = {")[0], 'check_threshold150_header', 'exec'))   # load(), table(), norm(), np, P, RES
RT = os.path.abspath(os.path.join(P, '..', '..', 'Undergrad_Glucose uptake signaling', 'RT-PCR'))
F_EXP = os.path.join(RT, '2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Cq Results.xlsx')
F_SUM = os.path.join(RT, '2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Summary.xlsx')
TOL = 0.0051
ex = table(F_EXP)
inst = {}                                              # (well, gene) -> (group, Cq)
for (g, k), items in ex.items():
    for w, c in items: inst[(w, g)] = (k, c)
ws = load(F_SUM)['0']
rows = []
for r in range(2, ws.max_row + 1):
    w, tg, sm, cq = (ws.cell(row=r, column=c).value for c in (2, 4, 6, 7))
    if w is None or not isinstance(cq, (int, float)): continue
    rows.append((r, str(w), str(tg), grp(sm), float(cq)))
blocks = {}
for (r, w, tg, k, cq) in rows: blocks.setdefault((norm(tg), k), []).append(cq)
iblocks = {key: sorted(c for _, c in items) for key, items in ex.items()}
out = []
for (r, w, tg, k, cq) in rows:
    g = norm(tg); ik, ic = inst[(w, g)]
    bv = sorted(blocks[(g, k)]); same_set = len(bv) == len(iblocks[(g, k)]) and all(abs(a - b) <= TOL for a, b in zip(bv, iblocks[(g, k)]))
    asc = all(blocks[(g, k)][i] <= blocks[(g, k)][i + 1] for i in range(len(blocks[(g, k)]) - 1))
    note = ''
    if g == 'gapdh' and k == 'DAPT':
        note = 'DAPT Gapdh: workbook holds the values of C10 and D10 twice; instrument wells A10 and B10 (24.987 and 24.132) are not in the workbook'
    out.append({'workbook_row': r, 'well': w, 'target': tg, 'group_in_workbook': k, 'workbook_Cq': round(cq, 6), 'instrument_group_of_named_well': ik, 'instrument_Cq_of_named_well': round(ic, 6),
                'abs_difference': round(abs(cq - ic), 6), 'workbook_Cq_equals_instrument_Cq_of_named_well': int(abs(cq - ic) <= TOL), 'block_value_set_equals_instrument': int(same_set),
                'block_ascending_down_the_workbook': int(asc), 'note': note})
with open(os.path.join(RES, 'aug26_workbook_vs_instrument.csv'), 'w', encoding='utf8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
# ---- parse threshold_fit.txt
tf = open(os.path.join(RES, 'threshold_fit.txt'), encoding='utf8').read().splitlines()
val = []
for ln in tf:
    m = re.match(r'\s+(After|Before)\s+(.+?)\s+known T ([\d.]+)\s+fitted T ([\d.]+)\s+mean \|recomputed - exported Cq\| ([\d.]+) \(max ([\d.]+)\)', ln)
    if m: val.append({'export_version': m.group(1), 'plate': m.group(2), 'known_threshold_RFU': m.group(3), 'fitted_threshold_RFU': m.group(4), 'mean_abs_difference_cycles': m.group(5), 'max_abs_difference_cycles': m.group(6)})
with open(os.path.join(RES, 'threshold_validation.csv'), 'w', encoding='utf8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(val[0].keys())); w.writeheader(); w.writerows(val)
res = []
for ln in tf:
    m = re.match(r'\s{2}(Wnt\w+|Lef1|Tcf7|Gapdh\(refA\))\s+(skipped|[\d.]+) / (skipped|[\d.]+) / (skipped|[\d.]+)\s*$', ln)
    if m:
        for k, v in zip(('CTRL', 'DAPT', 'DKK-1'), m.groups()[1:]): res.append({'assay': m.group(1), 'group': k, 'largest_abs_difference_of_4_sorted_values_cycles_at_fitted_threshold': v})
with open(os.path.join(RES, 'refA_block_residuals.csv'), 'w', encoding='utf8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(res[0].keys())); w.writeheader(); w.writerows(res)
rv = open(os.path.join(RES, 'refA_vs_plate1.txt'), encoding='utf8').read().splitlines()
dq = []
for ln in rv:
    m = re.match(r'(\w+)\s+(CTRL|DAPT|DKK-1)\s+([\d.]+) \|\s+([\d.]+)\s+([\d.]+)\s+([+-][\d.]+) \|\s+([\d.]+)\s+([\d.]+)\s+([+-][\d.]+)\s*$', ln)
    if m: dq.append({'gene': m.group(1), 'group': m.group(2), 'mean_Cq_sheet': m.group(3), 'dCq_sheet': m.group(4), 'dCq_reexport_45.70RFU': m.group(5), 'dCq_sheet_minus_reexport': m.group(6), 'mean_Cq_sheet_again': m.group(7), 'dCq_curves_at_209.91RFU': m.group(8), 'dCq_sheet_minus_curves': m.group(9)})
with open(os.path.join(RES, 'refA_dcq_comparison.csv'), 'w', encoding='utf8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(dq[0].keys())); w.writeheader(); w.writerows(dq)
print('rows: aug26 %d, validation %d, residuals %d, dcq %d' % (len(out), len(val), len(res), len(dq)))
print('equal named-well Cq: %d of %d; blocks with the same value set: %d of %d rows' % (sum(r['workbook_Cq_equals_instrument_Cq_of_named_well'] for r in out), len(out), sum(r['block_value_set_equals_instrument'] for r in out), len(out)))
