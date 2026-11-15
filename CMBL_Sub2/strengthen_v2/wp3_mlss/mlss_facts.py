"""Facts about the REPORTED MLSS run and about the DrugCentral target mapping that the README quotes.  Nothing is recomputed under an altered input here:
only the reported table (port) and the hub-to-drug mapping are examined.  Usage: python -I mlss_facts.py  (a second; writes mlss_facts.txt next to this file)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np, pandas as pd
from scipy import stats
import mlss_port as M
sys.stdout.reconfigure(encoding='utf8')
LOG = []
def log(s=''):
    LOG.append(str(s)); print(s)

dc = M.read_drugcentral(); pot = M.potency_from(dc)
hubs5, hub_df = M.hubs_from_centrality(min_degree=5); tg = M.targets_from(dc, hubs5)
REP = M.mlss_table(tg, pot, M.corr_from_matrix_csv())
R = pd.read_csv(os.path.join(M.RDIR, 'Fig_6c_MLSS_ALL_Combinations.csv'))
assert np.abs(R['MLSS_v4.0'].values - REP['MLSS_v4.0'].values).max() < 1e-9

def rank_r2(y, cols):
    yr = stats.rankdata(y); A = np.column_stack([np.ones(len(y))] + [np.asarray(c, float) for c in cols])
    beta = np.linalg.lstsq(A, yr, rcond=None)[0]; res = yr - A @ beta
    return float(1 - (res ** 2).sum() / ((yr - yr.mean()) ** 2).sum())

log('A. Reported run (port = R output; see gate_output.txt)')
REP['cls'] = np.where(REP['Complementarity_C'] >= 1 - 1e-12, 'cross-target', REP['Targets_A'] + 'x' + REP['Targets_B'])
for c, g in REP.groupby('cls'):
    log('  class %-18s n %3d  ranks %3d-%3d  C %s  B %s  V %s  MLSS %.4f to %.4f  distinct P %d' % (c, len(g), g['Rank'].min(), g['Rank'].max(), sorted(set(g['Complementarity_C'].round(6))),
        sorted(set(g['Balance_B'].round(6))), sorted(set(g['Network_Coverage_V'].round(6))), g['MLSS_v4.0'].min(), g['MLSS_v4.0'].max(), g['Potency_P'].round(9).nunique()))
s9 = REP['MLSS_v4.0'].round(9); vc = s9.value_counts()
log('  distinct MLSS values %d of 153; combinations sharing a score with another one: %d; largest tie groups: %s' % (s9.nunique(), int((vc[vc > 1]).sum()), ', '.join('%d at %.4f' % (n, v) for v, n in vc.head(4).items())))
drugs = sorted(tg, key=M.r_sort_key); withp = sorted(d for d in drugs if (d, tg[d][0]) in pot); nop = sorted(d for d in drugs if (d, tg[d][0]) not in pot)
log('  drugs %d; with a pIC50 %d; without %d: %s' % (len(drugs), len(withp), len(nop), ', '.join(nop)))
cr = REP[REP['cls'] == 'cross-target']
log('  cross-target combinations with Synergy > 0: %d (range %.4f to %.4f); with zero Synergy: %d' % ((cr['Synergy_Bonus'] > 0).sum(), cr.loc[cr['Synergy_Bonus'] > 0, 'Synergy_Bonus'].min(), cr['Synergy_Bonus'].max(), (cr['Synergy_Bonus'] == 0).sum()))
one = REP[(REP['Mean_pIC50_A'] == 0) ^ (REP['Mean_pIC50_B'] == 0)]
log('  combinations in which exactly one drug has no pIC50: %d; of these Potency_P > 0 in %d (P is computed from half the other drug\'s pIC50, it is not set to zero)' % (len(one), (one['Potency_P'] > 0).sum()))
log('  Antagonism Penalty > 0: %d combinations; negative r in the lookup: %d of %d pairs' % ((REP['Antagonism_Penalty'] > 0).sum(), sum(1 for v in M.corr_from_matrix_csv().values() if v < 0), len(M.corr_from_matrix_csv())))
rho_c = stats.spearmanr(cr['MLSS_v4.0'], cr['Mean_pIC50'])[0]; rho_a = stats.spearmanr(REP['MLSS_v4.0'], REP['Mean_pIC50'])[0]
isx = (REP['cls'] == 'cross-target').astype(float).values
log('  Spearman MLSS vs mean pIC50 of the pair: cross-target %.4f (R^2 %.4f); all 153 %.4f (R^2 %.4f)' % (rho_c, rho_c ** 2, rho_a, rho_a ** 2))
log('  rank R^2 of the MLSS: class (cross vs same target) alone %.4f; potency (rank of pair mean pIC50) alone %.4f; class + potency %.4f; class + potency + interaction %.4f' % (
    rank_r2(REP['MLSS_v4.0'], [isx]), rank_r2(REP['MLSS_v4.0'], [stats.rankdata(REP['Mean_pIC50'])]), rank_r2(REP['MLSS_v4.0'], [isx, stats.rankdata(REP['Mean_pIC50'])]),
    rank_r2(REP['MLSS_v4.0'], [isx, stats.rankdata(REP['Mean_pIC50']), isx * stats.rankdata(REP['Mean_pIC50'])])))
top = REP.head(15); ap = pd.concat([top['Drug_A'], top['Drug_B']]).value_counts()
log('  appearances in the top 15: ' + ', '.join('%s %d' % (k, v) for k, v in ap.items()))

log(''); log('B. DrugCentral rows for the genes of the network (drug.target.interaction.csv, %d rows)' % len(dc))
genes = ['GLP1R', 'CACNA1D', 'CACNA1C', 'TCF7L2', 'PTBP1', 'WNT5A', 'WNT5B', 'WNT9A', 'LEF1', 'TCF7', 'HES1', 'ITPR1', 'RYR2', 'RYR3', 'KCNJ11', 'SLC2A2', 'INS']
tok = dc['G'].fillna('').str.split('|')
for g in genes:
    ex = dc[dc['G'] == g]; inm = dc[tok.apply(lambda t: g in t and len(t) > 1)]
    log('  %-8s exact GENE rows %3d (human %3d); inside a multi-gene GENE field %3d (human %3d)' % (g, len(ex), (ex['ORGANISM'] == 'Homo sapiens').sum(), len(inm), (inm['ORGANISM'] == 'Homo sapiens').sum()))
mg = dc[tok.apply(lambda t: ('CACNA1D' in t or 'CACNA1C' in t) and len(t) > 1) & (dc['ORGANISM'] == 'Homo sapiens')]
log('  human multi-gene rows naming CACNA1D or CACNA1C: %d (GENE strings: %s); distinct drugs %d; among the 18 MLSS drugs %d; these rows are never matched by the exact test toupper(GENE) %%in%% hub genes' % (
    len(mg), ', '.join(sorted(set(mg['GENE']))), mg['DRUG_NAME'].nunique(), mg['DRUG_NAME'].isin(tg.keys()).sum()))
hum = dc[dc['ORGANISM'] == 'Homo sapiens']
only_mg = sorted(set(mg['DRUG_NAME']) - set(hum.loc[hum['G'] == 'CACNA1D', 'DRUG_NAME']), key=M.r_sort_key)
log('  drugs whose only human L-type channel record is the multi-gene field CACNA1D|CACNA1C, hence absent from the ranking (%d): %s' % (len(only_mg), ', '.join(only_mg)))
c1c = set(hum.loc[hum['G'] == 'CACNA1C', 'DRUG_NAME'])
log('  of the 18 MLSS drugs, those that also have an exact human CACNA1C record (CACNA1C has degree 4 and is not a hub): %d (%s)' % (len(c1c & set(tg)), ', '.join(sorted(c1c & set(tg), key=M.r_sort_key))))
log('  action types of the hub-gene human rows of the 18 drugs: ' + str(dc[dc['DRUG_NAME'].isin(tg.keys()) & dc['G'].isin(['CACNA1D', 'GLP1R']) & (dc['ORGANISM'] == 'Homo sapiens')]['ACTION_TYPE'].fillna('(blank)').value_counts().to_dict()))
log(''); log('C. Hub rule on the reported degrees: drugs that enter (hub-to-drug mapping only; the scores under these rules are in mlss_scenarios.csv once that run exists)')
for k in (3, 4, 5, 6):
    h, _ = M.hubs_from_centrality(min_degree=k); t = M.targets_from(dc, h)
    log('  degree >= %d: %d hubs %s; drugs %d; targets %s; combinations %d; drugs with more than one target %d' % (k, len(h), sorted(h), len(t), sorted({g for v in t.values() for g in v}), len(t) * (len(t) - 1) // 2, sum(len(v) > 1 for v in t.values())))
h20 = {'Cacna1d': 6, 'Glp1r': 6, 'Hey1': 6, 'Ptbp1': 6, 'Tcf7l2': 6}
log('  earlier 20-gene network (reanalysis_17genes/pipeline_20genes.json), degree >= 6 hubs: %s -> drug targets %s (same two as now)' % (sorted(h20), sorted({g for v in M.targets_from(dc, [x.upper() for x in h20]).values() for g in v})))
with open(os.path.join(HERE, 'mlss_facts.txt'), 'w', encoding='utf8') as f:
    f.write('\n'.join(LOG) + '\n')
