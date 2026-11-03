import sys, numpy as np, glob, os, re
from scipy.optimize import curve_fit
sys.stdout.reconfigure(encoding='utf8')
RE = r'C:\Users\saran\OneDrive\VSCBIC\C2F-2025\Publication\Submission\CMBL_Sub2_291026\Analysis_outputs\reanalysis_17genes'
sys.path.insert(0, RE)
import cfx
run = '2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Kcnj11-Cacna1d-Cacna1c-Ptbp1'
amp = cfx.sheet_rows(run, 'Quantification Amplification Results'); cq = cfx.sheet_rows(run, 'Quantification Cq Results')['0']; mp = cfx.sheet_rows(run, 'Melt Curve Peak Results')['0']
cqd = {(r[3], r[1].replace('0', '') if False else r[1]): (r[5], r[7]) for r in cq[1:] if r[1]}
def sig(x, top, mid, k, base): return base + top / (1 + np.exp(-(x - mid) / k))
print('Gapdh wells (26 Aug): sample, Cq, plateau RFU (cycles 36-40 mean), logistic mid / slope-width, fit R2')
rows = []
for gene, data in amp.items():
    if gene == 'Run Information': continue
    hdr = data[0][1:]; cyc = np.array([r[1] for r in data[1:]], float)
    for j, w in enumerate(hdr[1:], start=2):
        y = np.array([r[j] if isinstance(r[j], (int, float)) else np.nan for r in data[1:]], float)
        ok = ~np.isnan(y)
        try:
            p, _ = curve_fit(sig, cyc[ok], y[ok], p0=[y[ok].max(), 25, 1.5, 0], maxfev=20000); pred = sig(cyc[ok], *p); r2 = 1 - ((y[ok] - pred) ** 2).sum() / ((y[ok] - y[ok].mean()) ** 2).sum()
        except Exception: p = [np.nan] * 4; r2 = np.nan
        wl = w if len(w) == 3 else w[0] + '0' + w[1:]
        c = cqd.get((gene, wl), (None, None))
        rows.append(dict(gene=gene, well=w, sample=c[0], cq=c[1], plateau=float(np.nanmean(y[-5:])), top=p[0], mid=p[1], width=p[2], r2=r2))
import pandas as pd
d = pd.DataFrame(rows)
g = d[d.gene == 'Gapdh'].copy(); print(g.round(3).to_string(index=False))
print('\nplateau RFU, all target wells: median %.0f; Gapdh by sample:' % d[d.gene != 'Gapdh'].plateau.median(), g.groupby('sample').plateau.agg(['mean', 'min', 'max']).round(0).to_dict('index'))
# melt
print('\nGapdh melt peaks:'); [print('  ', r[1], r[5], 'Tm', r[6], 'height', r[7]) for r in mp[1:] if r[3] == 'Gapdh']
d.to_csv('qc26_wells.csv', index=False)
