import sys, numpy as np, pandas as pd
sys.stdout.reconfigure(encoding='utf8')
RE = r'C:\Users\saran\OneDrive\VSCBIC\C2F-2025\Publication\Submission\CMBL_Sub2_291026\Analysis_outputs\reanalysis_17genes'
sys.path.insert(0, RE)
import cfx
d = pd.read_csv('qc26_wells.csv')
d['row'] = d.well.str[0]
t = d[d.gene != 'Gapdh']
print('plateau RFU by plate row (target genes, 26 Aug):'); print(t.groupby('row').plateau.agg(['count', 'median', 'min']).round(0).to_string())
print('\ntarget wells with plateau < 700 RFU:'); print(t[t.plateau < 700][['gene', 'well', 'sample', 'cq', 'plateau']].round(2).to_string(index=False))
# Gapdh Cq SD per sample group per run (all runs)
out = []
for run in cfx.runs():
    try: rows = cfx.sheet_rows(run, 'Quantification Cq Results')['0']
    except Exception as e: continue
    h = rows[0]; iT, iS, iC = h.index('Target'), h.index('Sample'), h.index('Cq')
    df = pd.DataFrame([[r[iT], r[iS], r[iC]] for r in rows[1:] if r[iT]], columns=['t', 's', 'cq']); df = df[df.t.astype(str).str.lower().str.contains('gapdh')]
    df['cq'] = pd.to_numeric(df.cq, errors='coerce')
    for s, g in df.groupby('s'): out.append((run[:28], s, len(g.dropna()), round(g.cq.mean(), 2), round(g.cq.std(ddof=1), 2), round(g.cq.max() - g.cq.min(), 2)))
o = pd.DataFrame(out, columns=['run', 'sample', 'n', 'mean_Cq', 'SD', 'range']); print('\nGapdh replicate wells per group, SD (cycles):'); print(o.to_string(index=False))
print('\nmedian SD of all groups:', o.SD.median(), '| groups with SD > 0.5:', int((o.SD > 0.5).sum()), 'of', len(o))
