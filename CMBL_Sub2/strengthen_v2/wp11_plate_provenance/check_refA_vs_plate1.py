"""check_refA_vs_plate1.py - WP11 (8 Nov 2026), author item A12: the Wnt-family genes of the 24 June sheet (Wnt5a, Wnt5b, Wnt9a, Lef1, Tcf7 and their Gapdh vector 'refA') lie about 2.6 cycles above the 2025 re-export of 24 June plate 1
(file 2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 without ' 2', threshold 45.70 RFU). reconstruct_threshold.py finds that a single threshold of about 210 RFU applied to the same curves reproduces the sheet (mean 0.04 cycles).
This script asks the question that matters for the analysis: how much do the quantities that the analysis uses (dCq = target - Gapdh, and ddCq between groups) differ (a) between the sheet and the 2025 re-export (45.70 RFU) and (b) between the sheet and the curves read at the fitted threshold?
Group means of the four wells are used (the sheet carries no well labels for this run), so the comparison does not depend on how replicates are paired.
Read only. Output: results/refA_vs_plate1.txt.  Usage: python -I check_refA_vs_plate1.py"""
import sys; sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import os, csv
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'reconstruct_threshold.py'), encoding='utf8').read().split("lines = []")[0]
exec(compile(src, 'reconstruct_header', 'exec'))             # load(), norm(), read_plate(), cq_at(), np, SBA, W1, RES
T_FIT = 209.91                                              # fitted by reconstruct_threshold.py (results/threshold_fit.txt)
cyc, curves, tab = read_plate(os.path.join(SBA, '2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1'))
sheet = {}
for r in csv.DictReader(open(os.path.join(W1, 'ins1_sheet_cq.csv'), encoding='utf8')):
    sheet.setdefault((r['gene'], r['group']), []).append((int(r['replicate_in_group']), float(r['target_cq']), float(r['reference_cq'])))
G = ('CTRL', 'DAPT', 'DKK-1'); genesA = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7']
wells = {}
for (w, g, k, q) in tab: wells.setdefault((g, k), []).append((w, q))
lines = []


def log(s=''):
    lines.append(s); print(s)


def mean_export(g, k): return float(np.mean([q for _, q in wells[(g, k)]]))
def mean_T(g, k): return float(np.mean([cq_at(cyc, curves[w], T_FIT) for w, _ in wells[(g, k)]]))


gs = {k: float(np.mean([x[2] for x in sheet[('Wnt5a', k)]])) for k in G}            # refA is the reference vector of Wnt5a (WP1 gene-to-run map)
log('Gapdh group means (4 wells): sheet refA %s | 2025 re-export (45.70 RFU) %s | curves at %.2f RFU %s' % (' '.join('%s %.2f' % (k, gs[k]) for k in G), ' '.join('%s %.2f' % (k, mean_export('gapdh', k)) for k in G), T_FIT, ' '.join('%s %.2f' % (k, mean_T('gapdh', k)) for k in G)))
D = {}
log()
log('%-6s %-6s %9s | %9s %9s %8s | %9s %9s %8s' % ('gene', 'group', 'Cq sheet', 'dCq sheet', 'dCq 45.70', 'diff', 'Cq sheet', 'dCq %.0f' % T_FIT, 'diff'))
for g in genesA:
    for k in G:
        ts = float(np.mean([x[1] for x in sheet[(g, k)]]))
        ds = ts - gs[k]; de = mean_export(norm(g), k) - mean_export('gapdh', k); dt = mean_T(norm(g), k) - mean_T('gapdh', k)
        D[(g, k)] = (ds, de, dt)
        log('%-6s %-6s %9.2f | %9.2f %9.2f %+8.2f | %9.2f %9.2f %+8.2f' % (g, k, ts, ds, de, ds - de, ts, dt, ds - dt))
d1 = np.array([v[0] - v[1] for v in D.values()]); d2 = np.array([v[0] - v[2] for v in D.values()])
log()
log('dCq (target minus Gapdh), sheet minus re-export at 45.70 RFU, 15 gene-group means: mean %+.2f, SD %.2f, range %+.2f to %+.2f' % (d1.mean(), d1.std(ddof=1), d1.min(), d1.max()))
log('dCq (target minus Gapdh), sheet minus curves at %.2f RFU, 15 gene-group means: mean %+.2f, SD %.2f, range %+.2f to %+.2f' % (T_FIT, d2.mean(), d2.std(ddof=1), d2.min(), d2.max()))
log()
log('ddCq (group dCq minus CTRL dCq): sheet | re-export 45.70 | difference || curves at %.0f | difference' % T_FIT)
dd1, dd2 = [], []
for g in genesA:
    cells = []
    for k in ('DAPT', 'DKK-1'):
        s = D[(g, k)][0] - D[(g, 'CTRL')][0]; e = D[(g, k)][1] - D[(g, 'CTRL')][1]; t = D[(g, k)][2] - D[(g, 'CTRL')][2]
        dd1.append(s - e); dd2.append(s - t); cells.append('%s-CTRL %+.2f | %+.2f | %+.2f || %+.2f | %+.2f' % (k, s, e, s - e, t, s - t))
    log('  %-6s %s' % (g, '      '.join(cells)))
dd1, dd2 = np.array(dd1), np.array(dd2)
log('difference of the 10 ddCq values, sheet minus re-export (45.70 RFU): mean %+.2f, SD %.2f, largest absolute %.2f' % (dd1.mean(), dd1.std(ddof=1), np.abs(dd1).max()))
log('difference of the 10 ddCq values, sheet minus curves at %.2f RFU: mean %+.3f, SD %.3f, largest absolute %.3f' % (T_FIT, dd2.mean(), dd2.std(ddof=1), np.abs(dd2).max()))
sd_s = np.mean([np.std([x[1] for x in sheet[(g, k)]], ddof=1) for g in genesA for k in G]); sd_e = np.mean([np.std([q for _, q in wells[(norm(g), k)]], ddof=1) for g in genesA for k in G])
log('mean SD of the four target Cq values per block: sheet %.3f, re-export %.3f' % (sd_s, sd_e))
open(os.path.join(RES, 'refA_vs_plate1.txt'), 'w', encoding='utf8').write('\n'.join(lines) + '\n')
