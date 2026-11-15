"""check_aug26_nov09_order.py - WP11 (8 Nov 2026), author item A13: who put the 26 August rows in ascending Cq order, the instrument export or the workbook?
Inputs (read only): the instrument exports of the .pcrd files made in 2025 ('Quantification Cq Results', wells in plate order), the user's 26 August 'Quantification Summary' workbook (sheet '0': CFX columns B-H plus the user's calculation columns J-O),
the matrix rows of WP1 (../wp1_ins1_structure/ins1_sheet_cq.csv, gene_run_map.csv).
Outputs: results/aug26_nov09_order_check.txt and .json.  Usage: python -I check_aug26_nov09_order.py"""
import sys; sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import os, csv, json
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'check_threshold150.py'), encoding='utf8').read().split("RUNS = {")[0]
exec(compile(src, 'check_threshold150_header', 'exec'))      # load(), table(), norm(), grp(), np, openpyxl, paths (EXD, W1, RES)
RT = os.path.abspath(os.path.join(P, '..', '..', 'Undergrad_Glucose uptake signaling', 'RT-PCR'))
F_EXP = os.path.join(RT, '2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Cq Results.xlsx')
F_SUM = os.path.join(RT, '2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Summary.xlsx')
F_N9 = os.path.join(RT, '2020-11-09_ITPR, RYR1-3 -  Quantification Cq Results.xlsx')
TOL = 0.0051
lines = []; out = {}


def log(s=''):
    lines.append(s); print(s)


def wkey(w): return (int(w[1:]), w[0])                          # plate order: column by column, row A to H


meta = {r['gene']: r['run'] for r in csv.DictReader(open(os.path.join(W1, 'gene_run_map.csv'), encoding='utf8'))}
sheet = {}
for r in csv.DictReader(open(os.path.join(W1, 'ins1_sheet_cq.csv'), encoding='utf8')):
    sheet.setdefault((r['gene'], r['group']), []).append((int(r['replicate_in_group']), float(r['target_cq']), float(r['reference_cq'])))
G = ('CTRL', 'DAPT', 'DKK-1')
genes26 = [g for g in meta if meta[g] == 'Aug26']
genes9 = [g for g in meta if meta[g] == 'Nov09']
log('genes of the 26 August run in the matrix: ' + ', '.join(genes26) + '; of the 9 November run: ' + ', '.join(genes9))
# ------------------------------------------------------------------ A. instrument export, 26 August
ex = table(F_EXP)
blocks = {}
for g in genes26 + ['Gapdh']:
    for k in G:
        items = sorted(ex[(norm(g), k)], key=lambda x: wkey(x[0]))
        blocks[(g, k)] = items
asc = {key: all(items[i][1] <= items[i + 1][1] for i in range(len(items) - 1)) for key, items in blocks.items() if key[0] != 'Gapdh'}
log()
log('A. Instrument export of the 26 August plate (re-export of the .pcrd, wells in plate order)')
log('   target blocks (5 genes x 3 groups) whose four wells are in ascending Cq order in plate order: %d of %d (chance expectation about 1 in 24 per block, i.e. 0.6 of 15)' % (sum(asc.values()), len(asc)))
# ------------------------------------------------------------------ B. the user's Quantification Summary workbook
wsum = load(F_SUM)['0']
rows = []
for r in range(2, wsum.max_row + 1):
    w, tg, sm, cq = (wsum.cell(row=r, column=c).value for c in (2, 4, 6, 7))
    if w is None or not isinstance(cq, (int, float)): continue
    rows.append((r, str(w), norm(tg), grp(sm), float(cq)))
inst = {(w, norm(g)): c for (g, k), items in ex.items() for w, c in items}
tot = len(rows); same = sum(1 for (r, w, g, k, c) in rows if (w, g) in inst and abs(inst[(w, g)] - c) <= TOL)
tgt_rows = [x for x in rows if x[2] != 'gapdh']; ref_rows = [x for x in rows if x[2] == 'gapdh']
same_t = sum(1 for (r, w, g, k, c) in tgt_rows if abs(inst[(w, g)] - c) <= TOL); same_r = sum(1 for (r, w, g, k, c) in ref_rows if (w, g) in inst and abs(inst[(w, g)] - c) <= TOL)
log()
log("B. The user's 26 August 'Quantification Summary' workbook (sheet '0', Well and Cq columns) against the instrument value of the same well")
log('   rows compared: %d (60 target rows, 12 Gapdh rows); Cq equal to the instrument value of the named well: %d of %d (target rows %d of %d, Gapdh rows %d of %d)' % (tot, same, tot, same_t, len(tgt_rows), same_r, len(ref_rows)))
ordered = 0; nblk = 0
for g in genes26 + ['Gapdh']:
    for k in G:
        v = [c for (r, w, gg, kk, c) in rows if gg == norm(g) and kk == k]
        if len(v) == 4: nblk += 1; ordered += int(all(v[i] <= v[i + 1] for i in range(3)))
log('   blocks (gene x group) whose four Cq values run in ascending order down the sheet: %d of %d (Gapdh DAPT is the exception: its four cells hold two values twice)' % (ordered, nblk))
# block-wise comparison with the instrument: same set of values?
log('   block by block, sorted workbook values against sorted instrument values (tolerance 0.005):')
setres = {}
for g in genes26 + ['Gapdh']:
    cells = []
    for k in G:
        v = np.sort([c for (r, w, gg, kk, c) in rows if gg == norm(g) and kk == k]); e = np.sort([c for _, c in blocks[(g, k)]])
        ok = len(v) == len(e) and bool(np.all(np.abs(v - e) <= TOL)); setres[(g, k)] = ok; cells.append('same set' if ok else 'DIFFERENT')
    log('     %-8s %s' % (g, ' / '.join(cells)))
out['summary_workbook'] = {'rows': tot, 'cq_equal_to_instrument_value_of_named_well': same, 'target_rows': [same_t, len(tgt_rows)], 'gapdh_rows': [same_r, len(ref_rows)],
                           'blocks_same_value_set_as_instrument': sum(setres.values()), 'blocks_total': len(setres), 'blocks_different': sorted('%s %s' % key for key, ok in setres.items() if not ok)}
dap = [c for (r, w, gg, kk, c) in rows if gg == 'gapdh' and kk == 'DAPT']
log('   Gapdh DAPT: workbook cells %s; instrument wells %s' % (', '.join('%.3f' % c for c in dap), ', '.join('%s %.3f' % (w, c) for w, c in blocks[('Gapdh', 'DAPT')])))
# ------------------------------------------------------------------ C. the matrix rows of WP1
log()
log('C. The matrix rows (WP1, row order of the sheet) against the instrument export')
asc_m = 0; set_m = 0; nb = 0
for g in genes26:
    for k in G:
        v = [x[1] for x in sorted(sheet[(g, k)])]; e = [c for _, c in blocks[(g, k)]]
        nb += 1; asc_m += int(all(v[i] <= v[i + 1] for i in range(3))); set_m += int(np.all(np.abs(np.sort(v) - np.sort(e)) <= TOL))
log('   26 August target blocks ascending in the sheet row order: %d of %d; same value set as the instrument: %d of %d' % (asc_m, nb, set_m, nb))
refres = []
for k in G:
    v = [x[2] for x in sorted(sheet[(genes26[0], k)])]; e = [c for _, c in blocks[('Gapdh', k)]]
    refres.append((k, v, e))
    log('   Gapdh %-5s matrix reference vector %s | instrument wells %s' % (k, ' '.join('%.3f' % x for x in v), ' '.join('%.3f' % x for x in e)))
out['matrix_aug26'] = {'target_blocks_ascending': [asc_m, nb], 'target_blocks_same_set_as_instrument': [set_m, nb], 'plate_order_ascending_instrument': [sum(asc.values()), len(asc)]}
# ------------------------------------------------------------------ D. pairing of target and Gapdh wells (descriptive, post hoc)
log()
log('D. Pairing of target and Gapdh wells on the 26 August plate (descriptive; the reported analysis keeps the sheet pairing)')


def dcq_matrix(variant):
    M = []
    for g in genes26:
        row = []
        for k in G:
            if variant == 'sheet':                                     # the matrix as analyzed: sheet vectors, k-th smallest with k-th smallest
                v = [x for x in sorted(sheet[(g, k)])]; row += [t - r for (_, t, r) in v]
            else:
                tt = [c for _, c in blocks[(g, k)]]; rr = [c for _, c in blocks[('Gapdh', k)]]
                if variant == 'rank_instrument': tt, rr = sorted(tt), sorted(rr)
                row += list(np.array(tt) - np.array(rr))
        M.append(row)
    return np.array(M)


desc = {}
for name, var in [('as analyzed (sheet vectors, rank pairing, DAPT Gapdh wells A10 and B10 replaced)', 'sheet'), ('instrument values, rank pairing (all four Gapdh wells)', 'rank_instrument'), ('instrument values, plate-order pairing (all four Gapdh wells)', 'plate_instrument')]:
    M = dcq_matrix(var); Gm = M.reshape(len(genes26), 3, 4); res = (Gm - Gm.mean(axis=2, keepdims=True)).reshape(len(genes26), 12)
    C = np.corrcoef(res); iu = np.triu_indices(len(genes26), 1)
    means = Gm.mean(axis=2)
    desc[var] = {'within_group_sd_of_dCq': round(float(Gm.std(axis=2, ddof=1).mean()), 3), 'mean_residual_r_between_genes': round(float(C[iu].mean()), 3), 'mean_DAPT_minus_CTRL_dCq': round(float((means[:, 1] - means[:, 0]).mean()), 3)}
    log('   %-80s within-group SD of dCq %.3f | mean residual r between the 5 genes %+.3f | mean (DAPT - CTRL) dCq %+.3f' % (name, desc[var]['within_group_sd_of_dCq'], desc[var]['mean_residual_r_between_genes'], desc[var]['mean_DAPT_minus_CTRL_dCq']))
out['pairing_descriptive'] = desc
# ------------------------------------------------------------------ E. 9 November
ex9 = table(F_N9)
log()
log('E. 9 November run: matrix rows against the instrument export (exported order)')
same_seq = 0; same_set = 0; nb9 = 0; asc9_sheet = 0; asc9_exp = 0
for g in genes9:
    for k in G:
        v = [x[1] for x in sorted(sheet[(g, k)])]; e = [c for _, c in ex9.get((norm(g), k), [])]
        if len(e) != 4: continue
        nb9 += 1; same_seq += int(np.all(np.abs(np.array(v) - np.array(e)) <= TOL)); same_set += int(np.all(np.abs(np.sort(v) - np.sort(e)) <= TOL))
        asc9_sheet += int(all(v[i] <= v[i + 1] for i in range(3))); asc9_exp += int(all(e[i] <= e[i + 1] for i in range(3)))
log('   blocks compared: %d; same four values in the same row order as the export: %d; same set of values: %d; ascending in the sheet: %d; ascending in the export: %d' % (nb9, same_seq, same_set, asc9_sheet, asc9_exp))
out['nov09'] = {'blocks': nb9, 'same_sequence_as_export': same_seq, 'same_set_as_export': same_set, 'ascending_in_sheet': asc9_sheet, 'ascending_in_export': asc9_exp}
open(os.path.join(RES, 'aug26_nov09_order_check.txt'), 'w', encoding='utf8').write('\n'.join(lines) + '\n')
json.dump(out, open(os.path.join(RES, 'aug26_nov09_order_check.json'), 'w', encoding='utf8'), indent=1)
