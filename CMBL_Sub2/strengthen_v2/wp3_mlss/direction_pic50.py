"""Direction-aware table of the MLSS drugs and provenance of their pIC50 values (protocol section 5, last two sentences).
Usage: python -I direction_pic50.py   -> direction_table.csv, pic50_sources.csv (next to this file; inputs are read only).
direction_table.csv: one row per (drug, analysed-gene target in DrugCentral, human). 'used_in_mlss' = the row enters the MLSS target mapping
  (gene is a hub and the GENE field equals it exactly). The mRNA change of the target is the within-experiment fold change (group mean / CTRL mean
  of the relative expression) under DAPT and DKK-1, with the Bonferroni-adjusted p from Additional file 3; the DAPT value is also given
  under the two Gapdh recalculations (C0, C1) for the genes affected. 'relation' compares the direction of the drug action on the target
  (inhibit vs activate, from the DrugCentral ACTION_TYPE) with the direction of the mRNA change when that change has adjusted p < 0.05:
  inhibitor + lower mRNA, or activator + higher mRNA = 'adds to'; the opposite pairs = 'opposes'; otherwise 'no supported mRNA change'.
  The mRNA level is a proxy of target abundance, not of target activity; the table says nothing about the effect of the pair.
pic50_sources.csv: every DrugCentral row of the 18 MLSS drugs that names CACNA1D or GLP1R (any organism), with whether it was used for the pIC50."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np, pandas as pd
import mlss_port as M
sys.stdout.reconfigure(encoding='utf8')

dc = M.read_drugcentral(); pot = M.potency_from(dc)
hubs, _ = M.hubs_from_centrality(min_degree=5)
tg = M.targets_from(dc, hubs)
R = pd.read_csv(os.path.join(M.RDIR, 'Fig_6c_MLSS_ALL_Combinations.csv')).sort_values('Rank')
top15 = R.head(15); appear = pd.concat([top15['Drug_A'], top15['Drug_B']]).value_counts()
HUBDRUGS = set(appear[appear >= 4].index)                                   # same rule as Fig_6e/6f (>= 4 appearances in the top 15)

# ---------------- mRNA fold changes (matrix) and adjusted p (Additional file 3)
X, SRC = M.load_matrix()
def fc(Xv, g, rows):
    return float(Xv.loc[rows, g].mean() / Xv.loc[0:3, g].mean())
FC = {g: dict(dapt=fc(X, g, slice(4, 7)), dkk=fc(X, g, slice(8, 11))) for g in X.columns}
XC0, XC1 = M.gapdh_variant(X, 'C0'), M.gapdh_variant(X, 'C1')
G5 = ['Tcf7l2', 'Kcnj11', 'Cacna1d', 'Cacna1c', 'Ptbp1']
AF3 = os.path.join(os.path.dirname(M.ANALYSIS), 'Submission_package_JTM', 'Additional file 3.xlsx')
padj = {}
try:
    import openpyxl
    ws = openpyxl.load_workbook(AF3, read_only=True, data_only=True)['Analyzed genes (17)']
    rows = list(ws.iter_rows(values_only=True))
    hdr = [i for i, r in enumerate(rows) if r and r[0] == 'Gene'][0]
    for r in rows[hdr + 1:]:
        if r and r[0] in X.columns:
            padj[r[0]] = dict(dapt_mean=float(r[10]), dkk_mean=float(r[16]), p_dapt=float(r[23]), sig_dapt=r[24], p_dkk=float(r[27]), sig_dkk=r[28], test=r[2])
    chk = max(max(abs(padj[g]['dapt_mean'] - FC[g]['dapt']), abs(padj[g]['dkk_mean'] - FC[g]['dkk'])) for g in padj)
    print('Additional file 3 read (%d genes); largest difference between its group means and the matrix fold changes: %.4f' % (len(padj), chk))
except Exception as e:
    print('Additional file 3 not read (%s); adjusted p left empty' % e)

import json
GJ = json.load(open(M.GAPDH_JSON))                                           # per-gene Bonferroni p (CTRL vs DAPT) under the two Gapdh recalculations
def p_gapdh(g2, w):
    return GJ[w]['per_gene'][g2]['CTRL-DAPT'] if g2 in GJ[w]['per_gene'] else padj.get(g2, {}).get('p_dapt', np.nan)

GENE_MAP = {'CACNA1D': 'Cacna1d', 'CACNA1C': 'Cacna1c', 'GLP1R': 'Glp1r', 'KCNJ11': 'Kcnj11', 'RYR2': 'Ryr2', 'RYR3': 'Ryr3', 'SLC2A2': 'Glut2'}
INHIB = {'GATING INHIBITOR', 'BLOCKER', 'INHIBITOR', 'ANTAGONIST', 'NEGATIVE ALLOSTERIC MODULATOR', 'INVERSE AGONIST', 'ANTISENSE INHIBITOR', 'PORE BLOCKER', 'CHANNEL BLOCKER'}
ACTIV = {'AGONIST', 'ACTIVATOR', 'POSITIVE ALLOSTERIC MODULATOR', 'OPENER', 'PARTIAL AGONIST', 'POSITIVE MODULATOR'}
FIG6F_LABEL = {'isradipine': 'Antagonist', '(S)-nitrendipine': 'Antagonist', 'nitrendipine': 'Antagonist', 'albiglutide': 'Agonist', 'liraglutide': 'Agonist'}   # hard-coded in Fig_6f

def direction(actions):
    a = {x.upper() for x in actions}
    if a & INHIB and not a & ACTIV: return 'inhibit'
    if a & ACTIV and not a & INHIB: return 'activate'
    return 'not recorded' if not a else 'mixed'

def relation(drug_dir, p, fcv):
    if p is None or np.isnan(p) or p >= 0.05 or drug_dir not in ('inhibit', 'activate'):
        return 'no supported mRNA change' if (p is None or np.isnan(p) or p >= 0.05) else 'action type not recorded'
    lower = fcv < 1
    return 'adds to (same direction)' if ((drug_dir == 'inhibit') == lower) else 'opposes'

rows = []
cand = dc[(dc['ORGANISM'] == 'Homo sapiens') & dc['DRUG_NAME'].isin(tg.keys())]
for drug in sorted(tg, key=M.r_sort_key):
    sub = cand[cand['DRUG_NAME'] == drug]
    for gene, g2 in GENE_MAP.items():
        s = sub[sub['G'] == gene]
        if len(s) == 0:
            continue
        actions = sorted({a for a in s['ACTION_TYPE'].dropna()})
        used = gene in hubs
        d = direction(actions); f = FC[g2]; pa = padj.get(g2, {})
        rows.append(dict(drug=drug, hub_drug=drug in HUBDRUGS, appearances_in_top15=int(appear.get(drug, 0)), target_gene=gene, rat_gene_measured=g2, used_in_mlss=used,
                         drugcentral_action_type='|'.join(actions) if actions else 'not recorded', fig6f_label=FIG6F_LABEL.get(drug, ''), drug_direction=d,
                         fc_dapt=f['dapt'], p_adj_dapt=pa.get('p_dapt', np.nan), sig_dapt=pa.get('sig_dapt', ''), fc_dkk1=f['dkk'], p_adj_dkk1=pa.get('p_dkk', np.nan), sig_dkk1=pa.get('sig_dkk', ''),
                         fc_dapt_gapdh_C0=fc(XC0, g2, slice(4, 7)), fc_dapt_gapdh_C1=fc(XC1, g2, slice(4, 7)), p_adj_dapt_gapdh_C0=p_gapdh(g2, 'C0'), p_adj_dapt_gapdh_C1=p_gapdh(g2, 'C1'),
                         dapt_value_depends_on_gapdh=g2 in G5,
                         relation_dapt=relation(d, pa.get('p_dapt', np.nan), f['dapt']), relation_dkk1=relation(d, pa.get('p_dkk', np.nan), f['dkk']),
                         relation_dapt_gapdh_C0=relation(d, p_gapdh(g2, 'C0'), fc(XC0, g2, slice(4, 7))), relation_dapt_gapdh_C1=relation(d, p_gapdh(g2, 'C1'), fc(XC1, g2, slice(4, 7))),
                         note=('enters the MLSS' if used else 'not used in the MLSS (gene is not a hub)')))
D = pd.DataFrame(rows)
D = D.sort_values(['used_in_mlss', 'hub_drug', 'appearances_in_top15', 'drug'], ascending=[False, False, False, True]).reset_index(drop=True)
D.to_csv(os.path.join(HERE, 'direction_table.csv'), index=False, float_format='%.4g')
print('direction_table.csv: %d rows (%d used in the MLSS)' % (len(D), int(D['used_in_mlss'].sum())))
print(D[D['hub_drug'] & D['used_in_mlss']][['drug', 'target_gene', 'drugcentral_action_type', 'fc_dapt', 'p_adj_dapt', 'fc_dkk1', 'p_adj_dkk1', 'fc_dapt_gapdh_C1', 'p_adj_dapt_gapdh_C1', 'relation_dapt', 'relation_dkk1', 'relation_dapt_gapdh_C1']].to_string())

# ---------------- pIC50 provenance
HUBTGT = ['CACNA1D', 'GLP1R']
raw = dc[dc['DRUG_NAME'].isin(tg.keys()) & dc['G'].apply(lambda s: isinstance(s, str) and any(t in HUBTGT for t in s.split('|')))].copy()
best = raw[raw['ACT_TYPE'].isin(['IC50', 'Ki', 'EC50']) & (raw['ORGANISM'] == 'Homo sapiens') & (raw['RELATION'] == '=') & raw['ACT_VALUE'].notna() & raw['G'].isin(HUBTGT)]
mx = best.groupby(['DRUG_NAME', 'G'])['ACT_VALUE'].transform('max')
ok_idx = set(best.index[best['ACT_VALUE'] == mx])
def why(i, r):
    if i in ok_idx: return ''
    if r['ORGANISM'] != 'Homo sapiens': return 'target organism is %s, not Homo sapiens' % r['ORGANISM']
    if r['G'] not in HUBTGT: return 'GENE field lists several genes (%s); exact match to a hub gene fails' % r['GENE']
    if pd.isna(r['ACT_VALUE']): return 'no activity value (mechanism-of-action record)'
    if r['ACT_TYPE'] not in ('IC50', 'Ki', 'EC50'): return 'ACT_TYPE %s is not IC50, Ki or EC50' % r['ACT_TYPE']
    if r['RELATION'] != '=': return 'RELATION is not ='
    return 'lower than the maximum for this drug and gene'
def klass(r):
    c = (r['ACT_COMMENT'] if isinstance(r['ACT_COMMENT'], str) else '').lower()
    if pd.isna(r['ACT_VALUE']): return 'no value'
    if r['ACT_TYPE'] == 'EC50': return 'functional potency (EC50)' + ('; cAMP assay' if 'camp' in c else '; assay not described in the record' if not c else '')
    if r['ACT_TYPE'] == 'Ki': return 'binding constant (Ki)' + ('' if c else '; assay not described in the record')
    if 'displacement' in c or 'binding' in c: return 'binding affinity (IC50 of displacement)'
    if 'recombinant' in c or 'expressed' in c or 'flipr' in c: return 'IC50 in a recombinant expression system'
    return 'IC50' + ('' if c else '; assay not described in the record')
P = pd.DataFrame(dict(drug=raw['DRUG_NAME'], in_top15=raw['DRUG_NAME'].isin(set(top15['Drug_A']) | set(top15['Drug_B'])), hub_drug=raw['DRUG_NAME'].isin(HUBDRUGS),
                      gene_field=raw['GENE'], target_organism=raw['ORGANISM'], act_type=raw['ACT_TYPE'].fillna(''), act_value_pActivity=raw['ACT_VALUE'], relation=raw['RELATION'].fillna(''),
                      source_database=raw['ACT_SOURCE'].fillna(''), action_type=raw['ACTION_TYPE'].fillna(''), assay_description_in_record=raw['ACT_COMMENT'].fillna('').str.replace(r'<[^>]+>', '', regex=True),
                      assay_class=[klass(r) for _, r in raw.iterrows()], used_as_pIC50_in_mlss=[i in ok_idx for i in raw.index], reason_not_used=[why(i, r) for i, r in raw.iterrows()],
                      source_url=raw['ACT_SOURCE_URL'].fillna('')))
P['pIC50_used_by_the_R_script'] = P.apply(lambda r: pot.get((r['drug'], r['gene_field'].upper()), np.nan) if r['used_as_pIC50_in_mlss'] else np.nan, axis=1)
P = P.sort_values(['in_top15', 'drug', 'used_as_pIC50_in_mlss', 'target_organism'], ascending=[False, True, False, True]).reset_index(drop=True)
P.to_csv(os.path.join(HERE, 'pic50_sources.csv'), index=False)
U = P[P['used_as_pIC50_in_mlss']]
print('\npic50_sources.csv: %d DrugCentral rows for the 18 drugs; %d rows supply a pIC50 (%d drugs of 18; 7 have none)' % (len(P), len(U), U['drug'].nunique()))
print(U.groupby(['act_type', 'source_database']).size().to_string())
print('used values by drug:'); print(U[['drug', 'in_top15', 'act_type', 'act_value_pActivity', 'source_database', 'assay_class']].to_string(index=False))
nh = P[(~P['used_as_pIC50_in_mlss']) & P['act_value_pActivity'].notna()]
print('\nrows with a value that were NOT used: %d' % len(nh)); print(nh[['drug', 'gene_field', 'target_organism', 'act_type', 'act_value_pActivity', 'reason_not_used']].to_string(index=False))
