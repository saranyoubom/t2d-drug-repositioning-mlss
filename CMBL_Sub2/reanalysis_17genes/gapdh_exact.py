"""Exact re-computation of the 26 August DAPT values from the author's original Quantification Summary sheet (Undergrad_Glucose uptake signaling/RT-PCR).
Original analysis: DAPT replicate k of every target gene was paired with the Gapdh value HK_k = [22.5685, 22.7983, 22.5685, 22.7983] (raw wells C10 and D10, each used twice;
raw A10 24.9874 and B10 24.1318 were parked in column V and not used). Relative expression = 2^-(Cq_target - Gapdh).
 C0: reference = mean of the two valid wells (22.6834) for all four DAPT replicates (exclusion done without duplicating values)
 C1: reference = mean of all four raw wells (23.6225)   (worst case: failed wells included)
Run from this folder: python gapdh_exact.py"""
import sys, json
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, '.')
import pipeline as P
sys.stdout.reconfigure(encoding='utf8')
G5 = ['Tcf7l2', 'Kcnj11', 'Cacna1d', 'Cacna1c', 'Ptbp1']
HK = np.array([22.5685093010354, 22.7982726003435, 22.5685093010354, 22.7982726003435])      # as used (sheet rows 30-33)
RAW = {'A10': 24.9873996418614, 'B10': 24.1318270358876, 'C10': 22.5685093010354, 'D10': 22.7982726003435}
REF = {'C0': np.mean([RAW['C10'], RAW['D10']]), 'C1': np.mean(list(RAW.values()))}
df0 = P.load(P.GENES17)


def variant(ref):
    d = df0.copy()
    for g in G5: d.loc[4:7, g] = df0.loc[4:7, g].values * 2.0 ** (ref - HK)      # 2^-(Cq - ref) = 2^-(Cq - HK) * 2^(ref - HK)
    return d


def anova_bonf(x):
    g = [x[0:4], x[4:8], x[8:12]]; F_, p = stats.f_oneway(*g); mse = sum(((a - a.mean()) ** 2).sum() for a in g) / 9; out = {}
    for (i, j), nm in zip([(0, 1), (0, 2), (1, 2)], ['CTRL-DAPT', 'CTRL-DKK1', 'DAPT-DKK1']):
        t = (g[i].mean() - g[j].mean()) / np.sqrt(mse * 2 / 4); out[nm] = float(min(1.0, 2 * stats.t.sf(abs(t), 9) * 3))
    return dict(omnibus_p=float(p), fc_dapt=float(g[1].mean() / g[0].mean()), fc_dkk=float(g[2].mean() / g[0].mean()), **out)


def summary(d):
    sch = P.schemes(d); res = P.sig_sets(sch); cons = P.consensus(res, use_sign=False); adj = P.adjusted(sch); cn = P.counts(res); deg = {}
    for a, b in cons: deg[a] = deg.get(a, 0) + 1; deg[b] = deg.get(b, 0) + 1
    pr = ('Cacna1d', 'Glp1r')
    return dict(r_cacna1d_glp1r={k: float(res[k][pr][0]) for k in res}, p_cacna1d_glp1r={k: float(res[k][pr][1]) for k in res}, sig={k: v['sig'] for k, v in cn.items()}, consensus=[list(c) for c in cons], adj_q05={k: int(sum(v[2] < 0.05 for v in dd.values())) for k, dd in adj.items()},
                degree=dict(sorted(deg.items(), key=lambda kv: (-kv[1], kv[0]))))


out = {}
for name, d in [('reported', df0), ('C0', variant(REF['C0'])), ('C1', variant(REF['C1']))]:
    out[name] = dict(per_gene={g: anova_bonf(d[g].values) for g in G5}, **summary(d))
c0 = {tuple(x) for x in out['reported']['consensus']}
print('reference Cq: reported pairs HK', HK.round(3), '| C0 %.4f | C1 %.4f' % (REF['C0'], REF['C1']))
for g in G5:
    print('%-8s DAPT fold change / padj(CTRL-DAPT):  reported %.3f / %.2g | C0 %.3f / %.2g | C1 %.3f / %.2g' % (g, *[v for n in ('reported', 'C0', 'C1') for v in (out[n]['per_gene'][g]['fc_dapt'], out[n]['per_gene'][g]['CTRL-DAPT'])]))
for n in ('reported', 'C0', 'C1'):
    S = out[n]; cs = {tuple(x) for x in S['consensus']}
    print(f"[{n}] sig {S['sig']} | consensus {len(cs)} (kept {len(cs & c0)}, lost {len(c0 - cs)}, new {len(cs - c0)}) | adj q<0.05 {S['adj_q05']} | hubs {[g for g, v in S['degree'].items() if v >= 5]} | Glp1r {S['degree'].get('Glp1r', 0)} Cacna1d {S['degree'].get('Cacna1d', 0)} Ptbp1 {S['degree'].get('Ptbp1', 0)} Tcf7l2 {S['degree'].get('Tcf7l2', 0)}")
pairs = sorted(c0 | {tuple(x) for n in ('C0', 'C1') for x in out[n]['consensus']})
pd.DataFrame([dict(pair='-'.join(pr), involves_G5=bool(set(pr) & set(G5)), reported=pr in c0, C0=pr in {tuple(x) for x in out['C0']['consensus']}, C1=pr in {tuple(x) for x in out['C1']['consensus']}) for pr in pairs]).to_csv('gapdh_exact_pairs.csv', index=False)
json.dump(dict(ref={k: float(v) for k, v in REF.items()}, hk_used=HK.tolist(), raw=RAW, **out), open('gapdh_exact.json', 'w'), indent=1)
