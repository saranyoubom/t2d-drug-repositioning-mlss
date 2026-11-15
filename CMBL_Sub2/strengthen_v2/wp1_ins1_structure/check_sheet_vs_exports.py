"""WP1 provenance check (read only): do the target Cq values typed or computed on the authors' calculation sheets appear in the CFX 'Quantification Cq Results' exports that were kept?
Per gene and group the four sheet values are compared (sorted) with the four export values of the same target and sample. Writes sheet_vs_exports.csv.
Sources: 26 Aug and 9 Nov exports in the RT-PCR folder itself; 24 June 'Round2' and 3 June 'Gr1'/'Gr2' exports in SBA/Rev2_250626/Export_qPCR (re-exports made in 2025).
Usage: python -I check_sheet_vs_exports.py     (env RTPCR_DIR, EXPORT_DIR override the paths; about 10 s)"""
import sys, os, re, io, zipfile, glob
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd, openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
RT = os.environ.get('RTPCR_DIR', 'RT-PCR')
EXP = os.environ.get('EXPORT_DIR', 'Export_qPCR')
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


def export_table(path):
    ws = load(path)['0']; rows = list(ws.iter_rows(values_only=True)); h = [str(c) for c in rows[0]]; ti, si, ci = h.index('Target'), h.index('Sample'), h.index('Cq'); out = {}
    for r in rows[1:]:
        if r[ti] is None or not isinstance(r[ci], (int, float)): continue
        t, s = str(r[ti]), r[si]
        if '_' in t and (s is None or s == ''): t, s = t.rsplit('_', 1)                                                       # 24 June export: target is 'Lef1_Ctrl'
        out.setdefault((norm(re.sub(r'^\d+\.\s*', '', t)), grp(s)), []).append(float(r[ci]))
    return out


def sheet_values():
    """(run, gene, group) -> 4 sheet target Cq, from gene_run_map.csv / ins1_sheet_cq.csv (the matched rows of the calculation sheets)"""
    m = pd.read_csv(os.path.join(HERE, 'gene_run_map.csv')).set_index('gene'); cq = pd.read_csv(os.path.join(HERE, 'ins1_sheet_cq.csv')); out = {}
    for (g, k), d in cq.groupby(['gene', 'group']): out[(m.loc[g, 'run'], g, k)] = np.sort(d.target_cq.values)
    return out


EXPORTS = {'Aug26': [os.path.join(RT, '2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Knj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Cq Results.xlsx')], 'Nov09': [os.path.join(RT, '2020-11-09_ITPR, RYR1-3 -  Quantification Cq Results.xlsx')],
           'Jun24': glob.glob(os.path.join(EXP, '2020-Jun-24_Round2*Quantification Cq Results.xlsx')), 'Jun24 vs Gr1+Gr2': glob.glob(os.path.join(EXP, '2020-Jun-03_Gr*Quantification Cq Results.xlsx')),
           'Aug26 2025 re-export': glob.glob(os.path.join(EXP, '2020-08-26*Quantification Cq Results.xlsx')), 'Nov09 2025 re-export': glob.glob(os.path.join(EXP, '2020-11-09*Quantification Cq Results.xlsx'))}


def main():
    SV = sheet_values(); rows = []
    for key, paths in EXPORTS.items():
        run = key.split()[0]
        for p in ([paths] if key == 'Jun24 vs Gr1+Gr2' else [[q] for q in paths]):                                           # the two 3 June plates have 2 wells per group each: compare with their union
            ex = {}
            for q in p:
                for k, v in export_table(q).items(): ex.setdefault(k, []).extend(v)
            tag = 'Jun-03 Gr1+Gr2 (union, 2025 re-export)' if key == 'Jun24 vs Gr1+Gr2' else os.path.basename(p[0]).split(' - ')[0][:34] + (' (2025 re-export)' if ('2025' in key or key == 'Jun24') else ' (RT-PCR folder)')
            for (r, g, k), v in SV.items():
                if r != run: continue
                e = ex.get((norm(g), k))
                if e is None or len(e) < 4: continue
                e = np.sort(np.array(e)); n = min(len(e), len(v)); d = v[:n] - e[:n] if len(e) == 4 else np.nan
                rows.append(dict(run=r, gene=g, group=k, export=tag, n_export_values=len(e), sheet_values=' '.join('%.3f' % x for x in v), export_values=' '.join('%.3f' % x for x in e),
                                 n_values_equal_within_0p005=int(sum(np.isclose(x, e, atol=0.005).any() for x in v)), median_sheet_minus_export=float(np.median(v - e)) if len(e) == 4 else float('nan')))
    t = pd.DataFrame(rows); t.to_csv(os.path.join(HERE, 'sheet_vs_exports.csv'), index=False, float_format='%.4g')
    s = t.groupby(['run', 'export']).agg(blocks=('gene', 'size'), blocks_all4_equal=('n_values_equal_within_0p005', lambda x: int((x == 4).sum())), values_equal=('n_values_equal_within_0p005', 'sum'), median_diff=('median_sheet_minus_export', 'median'),
                                          min_diff=('median_sheet_minus_export', 'min'), max_diff=('median_sheet_minus_export', 'max'))
    pd.set_option('display.width', 220); print(s.round(3).to_string()); return 0


if __name__ == '__main__':
    sys.exit(main())
