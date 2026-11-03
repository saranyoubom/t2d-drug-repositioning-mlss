"""Additional file 2 (per-well efficiencies and melt curves): the data are kept as they are (all assays were measured); the 17-gene analysis only needs the four
excluded assays (Hey1, Ryr1, Wnt2, Wnt2b) to be marked. Writes 'Supplementary/Additional file 2_CMBL_working.xlsx'; the original is not touched."""
import sys
from copy import copy
import openpyxl
sys.stdout.reconfigure(encoding='utf8')
S = r'C:\Users\saran\OneDrive\VSCBIC\C2F-2025\Publication\Submission\CMBL_Sub2_291026'
SRC = S + r'\Supplementary\Additional file 2.xlsx'; OUT = S + r'\Supplementary\Additional file 2_CMBL_working.xlsx'
EXCL = ['Hey1', 'Ryr1', 'Wnt2', 'Wnt2b']
wb = openpyxl.load_workbook(SRC)

# flagged melt wells per gene, from the data
qc = wb['Melt_QC_per_well']; hdr = [c.value for c in qc[1]]; ig, inm, isec = hdr.index('Gene'), hdr.index('No_main_peak'), hdr.index('Secondary_signal')
flag = {}
for r in qc.iter_rows(min_row=2, values_only=True):
    if r[0] is None: continue
    f = flag.setdefault(r[ig], [0, 0]); f[0] += bool(r[inm]); f[1] += bool(r[isec])
flagged = {g: v for g, v in flag.items() if v[0] or v[1]}
print('flagged wells per gene (no main peak, secondary signal):', flagged)
assert set(flagged) <= set(EXCL), flagged

# Per_gene_summary: Analysis_status column
ws = wb['Per_gene_summary']; ncol = ws.max_column; hdr = [c.value for c in ws[1]]; assert 'Analysis_status' not in hdr
c0 = ws.cell(1, ncol + 1, 'Analysis_status'); c0._style = copy(ws.cell(1, ncol)._style)
for r in range(2, ws.max_row + 1):
    g = ws.cell(r, 1).value
    if g is None: continue
    c = ws.cell(r, ncol + 1, 'Excluded from the 17-gene analysis' if g in EXCL else ('Reference gene' if g == 'Gapdh' else 'Analysed')); c._style = copy(ws.cell(r, ncol)._style)
ws.column_dimensions[openpyxl.utils.get_column_letter(ncol + 1)].width = 34

# README: analysis status, and the melt result restated for the analysed assays
rd = wb['README']; rows = [c.value for c in rd['A']]
i = next(k for k, v in enumerate(rows) if isinstance(v, str) and v.startswith('Result: no well of the 20 target assays'))
fl = '; '.join(f'{g}: {v[1]} wells with a secondary signal, {v[0]} without a main peak' for g, v in sorted(flagged.items()))
rd.cell(i + 1, 1).value = ('Result: no well of the 17 analysed assays (and none of Gapdh) is flagged. The flagged wells belong to excluded assays (' + fl +
                        '). Script: R_SuppFig3_melt_curves_npg.R (public deposit, branch Rev3).')
add = ['Analysis status (added November 2026): all assays listed here were measured; Hey1, Ryr1, Wnt2 and Wnt2b are excluded from the analysed 17-gene set (Hey1: primer identity could not be established; Ryr1, Wnt2 and Wnt2b: more than 10% of wells with Cq >= 35 or no melt peak; details in Additional file 1, Excluded assays). Their wells are kept in every sheet so that the data are complete; sheet Per_gene_summary has a column Analysis_status.']
for k, t in enumerate(add):
    c = rd.cell(rd.max_row + 1, 1, t); c._style = copy(rd.cell(i + 1, 1)._style)
wb.save(OUT); print('saved', OUT)
