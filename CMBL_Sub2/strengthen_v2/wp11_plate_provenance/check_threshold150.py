"""check_threshold150.py - WP11 (8 Nov 2026): plate provenance of the 24 June / 3 June sheet values with the exports that the user re-made on 8 Nov 2026 from the native .pcrd files
(`Exported_data/Before` = threshold as found: auto 36.072 / auto 48.532 / manual 138.061; `Exported_data/After` = manual threshold 150 RFU for every run).
Question (author items A12): do the target and reference (Gapdh) Cq values that the calculation sheet 'Mix_3+24 June 2020' holds for the nine June genes of the analysis come from the plates
2020-Jun-03 Gr1 and Gr2 (two wells per group each) and 2020-Jun-24 Round2, once the threshold of the export is 150? Read only; writes results/threshold150_matches.csv and results/threshold150_summary.txt.
Sheet values: ../wp1_ins1_structure/ins1_sheet_cq.csv and gene_run_map.csv (the matched rows of the sheets, written by WP1). Usage: python -I check_threshold150.py"""
import os, re, io, sys, csv, zipfile, itertools
import numpy as np, openpyxl
HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
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


def table(path):
    ws = load(path)['0']; rows = list(ws.iter_rows(values_only=True)); h = [str(c) for c in rows[0]]
    wi, ti, si, ci = h.index('Well'), h.index('Target'), h.index('Sample'), h.index('Cq'); out = {}
    for r in rows[1:]:
        if r[ti] is None or not isinstance(r[ci], (int, float)): continue
        t, s = str(r[ti]), r[si]
        if '_' in t and (s is None or s == ''): t, s = t.rsplit('_', 1)
        out.setdefault((norm(re.sub(r'^\d+\.\s*', '', t)), grp(s)), []).append((str(r[wi]), float(r[ci])))
    return out


RUNS = {'Gr1': '2020-Jun-03_Gr1_13genes-Ctrl-DAPT-DKK1 -  Quantification Cq Results.xlsx', 'Gr2': '2020-Jun-03_Gr2_13genes-Ctrl-DAPT-DKK1 -  Quantification Cq Results.xlsx',
        'Round2': '2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1 2 -  Quantification Cq Results.xlsx'}
EX = {v: {r: table(os.path.join(EXD, v, f)) for r, f in RUNS.items()} for v in ('Before', 'After')}
sheet = {}
meta = {}
for r in csv.DictReader(open(os.path.join(W1, 'gene_run_map.csv'), encoding='utf8')):
    meta[r['gene']] = r['run']
for r in csv.DictReader(open(os.path.join(W1, 'ins1_sheet_cq.csv'), encoding='utf8')):
    sheet.setdefault((r['gene'], r['group']), []).append((int(r['replicate_in_group']), float(r['target_cq']), float(r['reference_cq'])))
genes = [g for g in meta if meta[g] == 'Jun24']
out = []
lines = []


def log(s=''):
    lines.append(s); print(s)


def eq(a, b, tol=0.0051):
    a, b = np.sort(np.array(a)), np.sort(np.array(b))
    return len(a) == len(b) and bool(np.all(np.abs(a - b) <= tol))


SOURCES = {'Round2': ['Round2'], 'Gr1+Gr2': ['Gr1', 'Gr2'], 'Gr1': ['Gr1'], 'Gr2': ['Gr2']}
for ver in ('After', 'Before'):
    log('=== exports with the threshold %s ===' % ('150 (After)' if ver == 'After' else 'as found (Before)'))
    for g in genes:
        res = []
        for k in ('CTRL', 'DAPT', 'DKK-1'):
            v = sorted(sheet[(g, k)])
            tv = [x[1] for x in v]
            hit = []
            for name, runs in SOURCES.items():
                ex = [c for rn in runs for (_, c) in EX[ver][rn].get((norm(g), k), [])]
                if len(ex) == 4 and eq(tv, ex): hit.append(name)
            res.append('+'.join(hit) if hit else '-')
            out.append([ver, g, k, ','.join('%.2f' % x for x in tv), ';'.join(hit) if hit else ''])
        log('  %-6s target Cq blocks equal to the export (CTRL / DAPT / DKK-1): %s' % (g, ' / '.join(res)))
    # Gapdh (reference) vectors per gene and group
    log('  reference (Gapdh) vector of each gene against the Gapdh wells of the exports:')
    for g in genes:
        res = []
        for k in ('CTRL', 'DAPT', 'DKK-1'):
            rv = [x[2] for x in sorted(sheet[(g, k)])]
            hit = []
            for name, runs in SOURCES.items():
                ex = [c for rn in runs for (_, c) in EX[ver][rn].get(('gapdh', k), [])]
                if len(ex) == 4 and eq(rv, ex): hit.append(name)
            res.append('+'.join(hit) if hit else '-')
        log('  %-6s reference blocks equal to Gapdh of: %s' % (g, ' / '.join(res)))
    log()
# shift caused by the threshold: After minus Before on the same wells
log('=== shift of the Cq values by the threshold (After minus Before, same wells, all targets) ===')
for rn in RUNS:
    d = []
    for key, wells in EX['After'][rn].items():
        b = dict(EX['Before'][rn].get(key, []))
        for (w, c) in wells:
            if w in b: d.append(c - b[w])
    d = np.array(d); log('  %-7s n %d  median %+.2f  range %+.2f to %+.2f' % (rn, len(d), np.median(d), d.min(), d.max()))
# replicate -> plate and well for the genes with an After match to Gr1+Gr2 or Round2
log('=== which plate and well carries each sheet replicate (After exports; value match within 0.005) ===')
for g in genes:
    parts = []
    for k in ('CTRL', 'DAPT', 'DKK-1'):
        v = sorted(sheet[(g, k)])
        cell = []
        for (rep, tcq, rcq) in v:
            m = [('%s:%s' % (rn, w)) for rn in RUNS for (w, c) in EX['After'][rn].get((norm(g), k), []) if abs(c - tcq) <= 0.0051]
            cell.append('%d=%s' % (rep, m[0] if len(m) == 1 else ('?' if not m else '/'.join(m))))
        parts.append(k + ' ' + ' '.join(cell))
    log('  %-6s %s' % (g, ' | '.join(parts)))
with open(os.path.join(RES, 'threshold150_matches.csv'), 'w', encoding='utf8', newline='') as f:
    w = csv.writer(f); w.writerow(['export_version', 'gene', 'group', 'sheet_target_cq_sorted', 'export_source_with_equal_values']); w.writerows(out)
open(os.path.join(RES, 'threshold150_summary.txt'), 'w', encoding='utf8').write('\n'.join(lines) + '\n')
