"""reconstruct_threshold.py - WP11 (8 Nov 2026): which threshold reproduces the values that the calculation sheet holds for the seven Wnt-family genes and their Gapdh (the 'refA' genes)?
Method: Cq is recomputed from the baseline-subtracted amplification curves ('Quantification Amplification Results') of a plate by linear interpolation of the first crossing of a threshold T.
Validation: for plates with a known threshold (the user's exports of 8 Nov: 'Before' as found and 'After' = 150) the fitted T must come back as that threshold, and the recomputed Cq must equal the exported Cq.
Then T is fitted to the sheet values (4 wells per gene and group on the 24 June plate 1, file 2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 without ' 2', 2025 re-export): the T with the smallest mean absolute difference of the sorted blocks.
Blocks in which a well never reaches T are skipped and listed (one well, G04 of the Wnt2b DAPT block, peaks at 89.6 RFU; an earlier version of this script demanded that every well cross T, which wrongly capped T at 89.6 and
made it look as if no single threshold could reproduce the sheet. The sheet values are matched without any well mapping: sorted block against sorted block.)
Read only. Usage: python -I reconstruct_threshold.py   (writes results/threshold_fit.txt)"""
import sys; sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import os, re, io, csv, zipfile
import numpy as np, openpyxl
HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
SBA = os.path.abspath(os.path.join(P, '..', 'SBA', 'Rev2_250626', 'Export_qPCR'))
EXD = os.path.join(P, 'Exported_data')
W1 = os.path.join(P, 'CMBL_Sub2', 'strengthen_v2', 'wp1_ins1_structure')
RES = os.path.join(HERE, 'results')
os.makedirs(RES, exist_ok=True)
CANON = {'[content_types].xml': '[Content_Types].xml', 'xl/sharedstrings.xml': 'xl/sharedStrings.xml', 'docprops/core.xml': 'docProps/core.xml', 'docprops/app.xml': 'docProps/app.xml'}


def fix(n):
    low = n.lower()
    if low in CANON: return CANON[low]
    if low.startswith('xl/worksheets/_rels/') and low.endswith('.xml.rels'): return 'xl/worksheets/_rels/' + n.split('/')[-1]
    if low.startswith('xl/worksheets/') and low.endswith('.xml'): return 'xl/worksheets/' + n.split('/')[-1]
    return n


def load(path):
    z = zipfile.ZipFile(path); buf = io.BytesIO(); zo = zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED)
    for n in z.namelist(): zo.writestr(fix(n), z.read(n))
    zo.close(); buf.seek(0); return openpyxl.load_workbook(buf, data_only=True)


def norm(s): return re.sub(r'[^a-z0-9]', '', str(s).lower())


def grp(s):
    s = str(s).lower(); return 'CTRL' if 'ctrl' in s or 'control' in s else 'DAPT' if 'dapt' in s else 'DKK-1' if 'dkk' in s else None


def read_plate(base):
    """base = path without the ' -  Quantification ...xlsx' suffix -> (cycles, curves {well: array}, cq table [(well, gene, group, cq)])"""
    wsa = load(base + ' -  Quantification Amplification Results.xlsx')
    ws = wsa[wsa.sheetnames[0]]; rows = list(ws.iter_rows(values_only=True)); head = [c for c in rows[0]]
    ci = head.index('Cycle'); cyc = np.array([r[ci] for r in rows[1:] if isinstance(r[ci], (int, float))], float)
    curves = {}
    for j, h in enumerate(head):
        if isinstance(h, str) and re.match(r'^[A-H]\d+$', h):
            w = h[0] + '%02d' % int(h[1:]); curves[w] = np.array([r[j] for r in rows[1:1 + len(cyc)]], float)
    wsc = load(base + ' -  Quantification Cq Results.xlsx'); wsc = wsc[wsc.sheetnames[0]] if '0' not in wsc.sheetnames else wsc['0']
    rr = list(wsc.iter_rows(values_only=True)); h = [str(c) if c is not None else '' for c in rr[0]]
    wi, ti, si, qi = h.index('Well'), h.index('Target'), h.index('Sample'), h.index('Cq'); tab = []
    for r in rr[1:]:
        if r[ti] is None or not isinstance(r[qi], (int, float)): continue
        t, s = str(r[ti]), r[si]
        if '_' in t and (s is None or s == ''): t, s = t.rsplit('_', 1)
        tab.append((str(r[wi]), norm(re.sub(r'^\d+\.\s*', '', t)), grp(s), float(r[qi])))
    return cyc, curves, tab


def cq_at(cyc, y, T, start=8):
    i = np.where((y[:-1] < T) & (y[1:] >= T) & (cyc[1:] > start))[0]
    if len(i) == 0: return np.nan
    i = i[0]
    return cyc[i] + (T - y[i]) / (y[i + 1] - y[i]) * (cyc[i + 1] - cyc[i])


lines = []


def log(s=''):
    lines.append(s); print(s)


def fit_T(cyc, curves, tab, grid):
    best = None
    for T in grid:
        d = [abs(cq_at(cyc, curves[w], T) - q) for (w, g, k, q) in tab if w in curves]
        d = np.array(d)
        if np.isnan(d).any(): continue
        s = d.mean()
        if best is None or s < best[0]: best = (s, T, d.max())
    return best


# 1. validation on the user's exports of 8 Nov (known thresholds)
log('=== validation: the fitted threshold must come back as the known one ===')
for ver, plate, known in [('After', '2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 2', 150.0), ('Before', '2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 2', 138.061), ('Before', '2020-Jun-03_Gr1_13genes-Ctrl-DAPT-DKK1', 36.072), ('After', '2020-Jun-03_Gr1_13genes-Ctrl-DAPT-DKK1', 150.0)]:
    cyc, curves, tab = read_plate(os.path.join(EXD, ver, plate))
    coarse = fit_T(cyc, curves, tab, np.arange(10, 400, 1.0)); fine = fit_T(cyc, curves, tab, np.arange(coarse[1] - 1.5, coarse[1] + 1.5, 0.01))
    log('  %-6s %-46s known T %.3f  fitted T %.3f  mean |recomputed - exported Cq| %.4f (max %.4f)' % (ver, plate, known, fine[1], fine[0], fine[2]))
# 2. the 24 June plate 1 (2025 re-export): threshold of the export, then the threshold of the sheet
base1 = os.path.join(SBA, '2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1')
cyc, curves, tab = read_plate(base1)
coarse = fit_T(cyc, curves, tab, np.arange(10, 400, 1.0)); fine = fit_T(cyc, curves, tab, np.arange(coarse[1] - 1.5, coarse[1] + 1.5, 0.01))
log()
log('=== 24 June plate 1 (file without " 2", 2025 re-export): threshold of the export = %.3f (mean |diff| %.4f, max %.4f) ===' % (fine[1], fine[0], fine[2]))
sheet = {}
for r in csv.DictReader(open(os.path.join(W1, 'ins1_sheet_cq.csv'), encoding='utf8')):
    sheet.setdefault((r['gene'], r['group']), []).append((float(r['target_cq']), float(r['reference_cq'])))
G = ('CTRL', 'DAPT', 'DKK-1')
genesA = ['Wnt2', 'Wnt2b', 'Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7']       # Wnt2 and Wnt2b are not among the 17 analyzed genes; they share the plate and the vector refA (Wnt2b has own DAPT values on the sheet)
wells = {}
for (w, g, k, q) in tab: wells.setdefault((g, k), []).append(w)
peaks = {w: float(np.nanmax(curves[w])) for w in curves}
weak = [(w, g, k) for (w, g, k, q) in tab if peaks[w] < 150]
log('wells per gene and group on this plate: ' + ', '.join('%s:%d' % (g, len(wells[(norm(g), 'CTRL')])) for g in genesA + ['Gapdh'] if (norm(g), 'CTRL') in wells))
log('peak baseline-subtracted RFU over the 96 wells of the seven assays and Gapdh: median %.0f, minimum %.1f; wells with a peak below 150 RFU: %s' % (np.median([peaks[w] for (w, g, k, q) in tab if g in [norm(x) for x in genesA] + ['gapdh']]), min(peaks[w] for (w, g, k, q) in tab if g in [norm(x) for x in genesA] + ['gapdh']), ', '.join('%s (%s %s, %.1f RFU)' % (w, g, k, peaks[w]) for (w, g, k) in weak)))
BL = [((norm(g), k), np.sort([x[0] for x in sheet[(g, k)]]), g) for g in genesA for k in G] + [(('gapdh', k), np.sort([x[1] for x in sheet[('Wnt5a', k)]]), 'Gapdh(refA)') for k in G]


def block_diff(T, sel=lambda name: True):
    d = []; skipped = []
    for (key, s, name) in BL:
        if not sel(name): continue
        e = np.array([cq_at(cyc, curves[w], T) for w in wells[key]])
        if np.isnan(e).any(): skipped.append('%s %s' % (name, key[1])); continue
        d.append(np.abs(np.sort(e) - s))
    return d, skipped


def mean_cost(T, sel=lambda name: True):
    d, sk = block_diff(T, sel)
    return (np.concatenate(d).mean() if d else np.inf), len(sk)


def best_T(sel=lambda name: True, lo=60, hi=420):
    grid = np.arange(lo, hi, 1.0); c = [mean_cost(T, sel) for T in grid]
    ok = [(v, T) for (v, n), T in zip(c, grid) if np.isfinite(v)]
    # a threshold that skips many blocks is not a valid fit: require that at most 2 of the 24 blocks are skipped
    ok = [(v, T) for (v, n), T in zip(c, grid) if np.isfinite(v) and n <= 2]
    v0, T0 = min(ok); fine = np.arange(T0 - 1.5, T0 + 1.5, 0.01); cf = [(mean_cost(T, sel), T) for T in fine]
    (v1, n1), T1 = min([x for x in cf if x[0][1] <= 2], key=lambda x: x[0][0]); return T1, v1, n1


T1, v1, n1 = best_T()
log('threshold that best reproduces the sheet values (seven assays and Gapdh refA, blocks of four wells): T = %.2f RFU, mean |recomputed - sheet| = %.4f cycles over %d values (blocks skipped because a well never reaches T: %d)' % (T1, v1, 4 * (24 - n1), n1))
log('mean |recomputed - sheet| at other thresholds: ' + '; '.join('T %.0f: %.3f' % (T, mean_cost(T)[0]) for T in (45.7, 100, 150, 180, 200, T1, 220, 250)))
log('thresholds fitted separately: Gapdh blocks T = %.2f (mean %.4f); target blocks T = %.2f (mean %.4f); the five analyzed genes only (Wnt5a, Wnt5b, Wnt9a, Lef1, Tcf7) T = %.2f (mean %.4f)' % (*best_T(lambda n: n == 'Gapdh(refA)')[:2], *best_T(lambda n: n != 'Gapdh(refA)')[:2], *best_T(lambda n: n in ('Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7'))[:2]))
log('per block at T = %.2f (largest absolute difference of the four sorted values, in cycles; CTRL / DAPT / DKK-1):' % T1)
for g in genesA + ['Gapdh(refA)']:
    row = []
    for k in G:
        key = [b for b in BL if b[2] == g and b[0][1] == k][0]
        e = np.array([cq_at(cyc, curves[w], T1) for w in wells[key[0]]])
        row.append('skipped' if np.isnan(e).any() else '%.3f' % np.abs(np.sort(e) - key[1]).max())
    log('  %-12s %s' % (g, ' / '.join(row)))
open(os.path.join(RES, 'threshold_fit.txt'), 'w', encoding='utf8').write('\n'.join(lines) + '\n')
