"""wp13_gapdh_versions.py - WP13 (9 Nov 2026), at the user's request ("I do not know why Gapdh wells A10 and B10 were replaced; can you calculate the corrected version"):
the 26 August DAPT reference of the INS-1 matrix in three versions, carried through the pooled analysis and the Fig. 5 classification. No new rule: the versions were defined on 2 Nov 2026 in
Analysis_outputs/reanalysis_17genes/gapdh_exact.py and are used unchanged.
  reported : the DAPT replicates k = 1..4 of the five 26 August genes use the Gapdh values [22.5685, 22.7983, 22.5685, 22.7983] (wells C10 and D10, each entered twice; A10 and B10 not used)
  C0       : A10 and B10 excluded WITHOUT duplicating: the mean of the two valid wells (22.6834) for all four DAPT replicates   ('corrected by exclusion')
  C1       : all four wells (mean 23.6215) for all four DAPT replicates                                                        ('corrected, failed wells included')
Reads only; writes into its own folder (results/). Usage: python -I wp13_gapdh_versions.py"""
import sys, os, json, csv
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.abspath(__file__))
STR = os.path.abspath(os.path.join(HERE, '..'))
RE = os.path.abspath(os.path.join(HERE, '..', '..', 'reanalysis_17genes'))
RES = os.path.join(HERE, 'results'); os.makedirs(RES, exist_ok=True)
os.chdir(RE); sys.path.insert(0, RE); sys.path.insert(0, os.path.join(STR, 'wp10_versions'))      # pipeline.py loads its data with paths relative to its own folder
import numpy as np, pandas as pd
import pipeline as P
import string_sensitivity as ss
G5 = ['Tcf7l2', 'Kcnj11', 'Cacna1d', 'Cacna1c', 'Ptbp1']
HK = np.array([22.5685093010354, 22.7982726003435, 22.5685093010354, 22.7982726003435])
RAW = {'A10': 24.9873996418614, 'B10': 24.1318270358876, 'C10': 22.5685093010354, 'D10': 22.7982726003435}
REF = {'C0': float(np.mean([RAW['C10'], RAW['D10']])), 'C1': float(np.mean(list(RAW.values())))}
df0 = P.load(P.GENES17)


def variant(ref):
    d = df0.copy()
    for g in G5: d.loc[4:7, g] = df0.loc[4:7, g].values * 2.0 ** (ref - HK)          # 2^-(Cq - ref) = 2^-(Cq - HK) * 2^(ref - HK)
    return d


old, cons_rep, pairs_rep = ss.load_inputs()
oedges = [(r['node1'], r['node2']) for r in old]
sset_alias, sset_lit = ss.string_set(oedges, True), ss.string_set(oedges, False)
out = {}
for name, d in [('reported', df0), ('C0', variant(REF['C0'])), ('C1', variant(REF['C1']))]:
    sch = P.schemes(d); res = P.sig_sets(sch); cons = P.consensus(res, use_sign=False); cn = P.counts(res); adj = P.adjusted(sch)
    deg = {}
    for a, b in cons: deg[a] = deg.get(a, 0) + 1; deg[b] = deg.get(b, 0) + 1
    h = res['Hes1-normalized']
    p07 = {'-'.join(pr): r for pr, (r, p) in h.items() if p < 0.05 and abs(r) >= 0.7}
    consset = {'-'.join(sorted(pr)) for pr in cons}
    cl_a = ss.classify(p07, consset, sset_alias); cl_l = ss.classify(p07, consset, sset_lit)
    # fold changes of the DAPT group for the five genes (per-replicate 2^-dCt, relative to the CTRL mean)
    fc = {g: float(d[g].values[4:8].mean() / d[g].values[0:4].mean()) for g in G5}
    out[name] = dict(consensus=[list(c) for c in cons], n_consensus=len(cons), sig={k: v['sig'] for k, v in cn.items()}, adj_q05={k: int(sum(v[2] < 0.05 for v in dd.values())) for k, dd in adj.items()},
                     degree=dict(sorted(deg.items(), key=lambda kv: (-kv[1], kv[0]))), hubs=[g for g, k in deg.items() if k >= 5], n_pairs_abs_r_ge_0_7=len(p07),
                     fig5_alias=ss.counts(cl_a), fig5_literal=ss.counts(cl_l), in_STRING_alias=sum(1 for p in p07 if p in sset_alias), fc_dapt=fc,
                     r_cacna1d_glp1r=float(h[('Cacna1d', 'Glp1r')][0]))
base = {tuple(x) for x in out['reported']['consensus']}
for name in out:
    cs = {tuple(x) for x in out[name]['consensus']}
    out[name]['kept'] = len(cs & base); out[name]['lost'] = sorted('-'.join(p) for p in base - cs); out[name]['new'] = sorted('-'.join(p) for p in cs - base)
# QC evidence from the earlier session's well table (plateau fluorescence of the 12 Gapdh wells)
qc = [r for r in csv.DictReader(open(os.path.join(RE, 'qc26_wells.csv'), encoding='utf8')) if r['gene'] == 'Gapdh']
out['qc_gapdh_wells'] = [dict(well=r['well'], sample=r['sample'], cq=round(float(r['cq']), 3), plateau_RFU=round(float(r['plateau']), 0)) for r in qc]
out['ref'] = REF; out['hk_used'] = HK.tolist(); out['raw'] = RAW
json.dump(out, open(os.path.join(RES, 'gapdh_versions.json'), 'w', encoding='utf8'), indent=1, ensure_ascii=False)
with open(os.path.join(RES, 'gapdh_versions_summary.csv'), 'w', encoding='utf8', newline='') as f:
    w = csv.writer(f); w.writerow(['version', 'DAPT_reference_Gapdh_Cq', 'consensus_pairs', 'kept_of_21', 'lost', 'new', 'hubs_degree_ge_5', 'pooled_sig_raw/Ins1/Kcnj11/Hes1', 'adjusted_q05_raw/Ins1/Kcnj11/Hes1', 'r_Cacna1d_Glp1r_Hes1norm',
                                    'fig5_alias_KnownValidated/Novel/STRINGOnly/NotValidated', 'in_STRING_alias_of_120', 'fc_DAPT_Tcf7l2/Kcnj11/Cacna1d/Cacna1c/Ptbp1'])
    for n, ref in (('reported', 'rows: 22.5685, 22.7983, 22.5685, 22.7983'), ('C0', '%.4f' % REF['C0']), ('C1', '%.4f' % REF['C1'])):
        o = out[n]; c = o['fig5_alias']
        w.writerow([n, ref, o['n_consensus'], o['kept'], '; '.join(o['lost']), '; '.join(o['new']), '; '.join('%s %d' % (g, o['degree'][g]) for g in o['hubs']), '/'.join(str(o['sig'][k]) for k in o['sig']), '/'.join(str(o['adj_q05'][k]) for k in o['adj_q05']),
                    round(o['r_cacna1d_glp1r'], 3), '%d/%d/%d/%d' % (c['Known Validated'], c['Novel Discovery'], c['STRING Only'], c['Not Validated']), o['in_STRING_alias'], '/'.join('%.2f' % o['fc_dapt'][g] for g in G5)])
for n in ('reported', 'C0', 'C1'):
    o = out[n]; c = o['fig5_alias']; l = o['fig5_literal']
    print('[%s] consensus %d (kept %d of 21; lost %s; new %s) | hubs %s | pooled sig %s | adj q<0.05 %s | |r|>=0.7 pairs %d | Fig5 alias %s (in STRING %d) | Fig5 literal %s | r(Cacna1d,Glp1r) %.3f' % (
        n, o['n_consensus'], o['kept'], o['lost'], o['new'], {g: o['degree'][g] for g in o['hubs']}, list(o['sig'].values()), list(o['adj_q05'].values()), o['n_pairs_abs_r_ge_0_7'], c, o['in_STRING_alias'], l, o['r_cacna1d_glp1r']))
    print('     DAPT fold change vs CTRL:', {g: round(v, 3) for g, v in o['fc_dapt'].items()})
print('Gapdh wells (plateau RFU):', [(r['well'], r['sample'], r['cq'], r['plateau_RFU']) for r in out['qc_gapdh_wells']])
