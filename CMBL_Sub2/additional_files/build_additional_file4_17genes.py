"""Rebuild 'Supplementary/Additional file 4.xlsx' for the 17-gene analysis (Hey1, Ryr1, Wnt2, Wnt2b excluded).

Input : the 20-gene workbook (kept as 'Additional file 4_before_17gene.xlsx'), the 17-gene CSV outputs in R_17genes/, the 7-gene Fig 1d statistics in R_panels_authors/,
        and the Python pooled-pipeline results in reanalysis_17genes/ (permutation arrays, pipeline.py).
Output: 'Supplementary/Additional file 4_CMBL_working.xlsx' (the 20-gene file is not touched).
Sheets that do not depend on the gene set (Fig 1c, 2b-f, 3a-f, 3g, 6c-6e, Additional file 1 tables, raw data) are copied unchanged; only their wording is adjusted where it named excluded genes or "hub genes".
"""
import sys, os, re, math, shutil, json
from copy import copy
import numpy as np, pandas as pd, openpyxl
from openpyxl.utils import get_column_letter
sys.stdout.reconfigure(encoding='utf8')
S = r'C:\Users\saran\OneDrive\VSCBIC\C2F-2025\Publication\Submission\CMBL_Sub2_291026'
A = S + r'\Analysis_outputs'; R17 = A + r'\R_17genes'; DX = R17 + r'\figures_publication\data_export'; RA = A + r'\R_panels_authors'; RE = A + r'\reanalysis_17genes'
OLD = S + r'\Supplementary\Additional file 4.xlsx'; OUT = S + r'\Supplementary\Additional file 4_CMBL_working.xlsx'
DATE = '2026-11-03'
sys.path.insert(0, RE)
import pipeline as P
from statsmodels.stats.multitest import multipletests

wb = openpyxl.load_workbook(OLD)
tpl = wb['Fig_4b_full_correlation_matrix']
ST = {k: copy(tpl[c]._style) for k, c in dict(title='B1', desc='B2', header='B5', text='B6', num='C6').items()}
EXCL = {'Hey1', 'Wnt2', 'Wnt2b', 'Ryr1'}
REL5 = {'Known Validated': 'Known and analytically consistent', 'Novel Discovery': 'Analytically consistent, database-unregistered pairs',
        'STRING Only': 'STRING Only', 'Not Validated': 'Not analytically consistent'}
log = []


def conv(v):
    if v is None: return None
    if isinstance(v, (int, float, bool)): return v
    s = str(v).replace('\r\n', '\n')
    if s in ('NA', ''): return None
    if s in ('TRUE', 'FALSE'): return s == 'TRUE'
    if re.fullmatch(r'-?\d+', s): return int(s)
    if re.fullmatch(r'-?\d*\.\d+(e-?\d+)?|-?\d+(\.\d*)?e-?\d+', s, re.I): return float(s)
    return s


def csv_rows(path, rename=None):
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    hdr = list(d.columns); hdr[0] = hdr[0] or 'Gene'
    if rename: hdr = [rename.get(h, h) for h in hdr]
    return hdr, [[conv(x) for x in r] for r in d.values.tolist()]


def rebuild(name, desc, source, header, rows, key_width=None, hdr=5):
    idx = wb.sheetnames.index(name); old = wb[name]
    ws = wb.create_sheet(title=f'tmp{idx}', index=idx)
    ncol = len(header); last = get_column_letter(1 + ncol)
    ws['B1'] = name; ws['B1']._style = copy(ST['title'])
    ws['B2'] = desc; ws['B3'] = source
    for r in (2, 3):
        ws.cell(r, 2)._style = copy(ST['desc']); ws.merge_cells(f'B{r}:{last}{r}')
    for j, h in enumerate(header):
        c = ws.cell(hdr, 2 + j, h); c._style = copy(ST['header'])
    for i, row in enumerate(rows):
        for j, v in enumerate(row):
            c = ws.cell(hdr + 1 + i, 2 + j, v)
            if isinstance(v, bool) or v is None or isinstance(v, str): c._style = copy(ST['text'])
            else:
                c._style = copy(ST['num'])
                if isinstance(v, int): c.number_format = 'General'
                elif abs(v) < 1e-4 and v != 0: c.number_format = '0.00E+00'
    widths = []
    for j in range(ncol):
        vals = [str(header[j])] + [str(r[j]) for r in rows[:400] if r[j] is not None]
        w = min(max(max(len(v.split('\n')[0]) for v in vals) + 2, 9), 46)
        widths.append(w)
    if key_width: widths[0] = key_width
    ws.column_dimensions['A'].width = 3.7
    for j, w in enumerate(widths): ws.column_dimensions[get_column_letter(2 + j)].width = w
    tw = sum(widths)
    for r, t in ((2, desc), (3, source)): ws.row_dimensions[r].height = 13.5 * max(1, math.ceil(len(t) / max(tw * 1.1, 40))) + 2
    ws.freeze_panes = f'B{hdr + 1}'
    wb.remove(old); ws.title = name
    log.append((name, len(rows), ncol))
    return ws


def old_rows(name, min_row=5):
    ws = wb[name]; rows = [[c.value for c in r[1:]] for r in ws.iter_rows(min_row=min_row) if any(c.value is not None for c in r[1:])]
    while rows and rows[-1] and all(v is None for v in rows[-1]): rows.pop()
    return rows[0], rows[1:]


def set_desc(name, desc=None, source=None):
    ws = wb[name]
    if desc: ws['B2'] = desc
    if source: ws['B3'] = source


def src(csv): return f'Source: {csv}  |  Generated: {DATE}'


# ------------------------------------------------------------------ Fig 1d (7 genes)
h, rows = csv_rows(RA + r'\Fig_1d_stats.csv')
oh, orow = old_rows('Fig_1d_stats')
def blocks(rs):
    out = {}; cur = None
    for r in rs:
        if r[0] and r[0] != '---' and r[0] != 'Section' and ' - ' in str(r[0]): cur = str(r[0]).split(' - ')[1].split(' ')[0]
        if cur: out.setdefault(cur, []).append([conv(x) for x in r])
    return out
bo, bn = blocks(orow), blocks(rows)
bad = 0
for g in bn:
    a, b = bo[g], bn[g]
    if len(a) != len(b): bad += 1; continue
    for ra, rb in zip(a, b):
        for x, y in zip(ra[1:], rb[1:]):
            if isinstance(x, float) or isinstance(y, float):
                try:
                    if abs(float(x) - float(y)) > 1e-4 * max(1, abs(float(x))): bad += 1
                except Exception: bad += 1
            elif (x or None) != (y or None): bad += 1
print(f'Fig 1d: per-gene blocks of the 7-gene file vs the 10-gene sheet: {len(bn)} genes compared, mismatching cells {bad}')
assert bad == 0
rebuild('Fig_1d_stats', 'Fig 1d: Notch/Wnt gene expression, combined statistical results for the seven analysed genes (Hes1, Wnt5a, Wnt5b, Wnt9a, Lef1, Tcf7, Tcf7l2). Hey1, Wnt2 and Wnt2b were measured but are not analysed (excluded from the 17-gene analysis).',
        src('Fig_1d_stats.csv (R_panels_authors)'), h, rows)
h2, r2 = old_rows('Fig_1d_Statistical_Results'); gi = h2.index('Gene')
r2 = [r for r in r2 if r[gi] not in EXCL]
rebuild('Fig_1d_Statistical_Results', 'Fig 1d: Per-gene statistical results table (seven analysed genes): test statistic, p-value, adjusted p-value, significance label',
        f'Source: Fig_1d_Statistical_Results.csv  |  Generated: 2026-04-23; rows of Hey1, Wnt2 and Wnt2b removed {DATE} (per-gene results are unchanged)', h2, r2)
h3, r3 = old_rows('Fig_1d_Test_Type_Per_Gene'); r3 = [r for r in r3 if r[0] not in EXCL]
rebuild('Fig_1d_Test_Type_Per_Gene', 'Fig 1d: Test type selected per gene (ANOVA or Kruskal-Wallis) based on Shapiro-Wilk normality results (seven analysed genes)',
        f'Source: Fig_1d_Test_Type_Per_Gene.csv  |  Generated: 2026-04-23; rows of Hey1, Wnt2 and Wnt2b removed {DATE}  |  use_param = False: Kruskal-Wallis with Bonferroni-corrected Wilcoxon tests (Wnt5a, Lef1); True: one-way ANOVA with Bonferroni-corrected t-tests. See Additional file 1: Table S4 for both procedures.', h3, r3)

# ------------------------------------------------------------------ Fig 3e description (Ryr1 never in the sheet; Ptbp1 is in Fig 3f)
set_desc('Fig_3e_stats', 'Fig 3e: Calcium channel genes, combined statistics for Ryr2, Ryr3 and Itpr1 (Ptbp1 is in sheet Fig_3f_Ptbp1_stats)')

# ------------------------------------------------------------------ Fig 4b-4e (CSV outputs of the 17-gene R run)
sets = [
 ('Fig_4b_full_correlation_matrix', 'Fig_4b_full_correlation_matrix.csv', 'Fig 4b: Full Pearson correlation matrix (17 × 17 genes), Hes1-normalized expression data, all pairwise r values'),
 ('Fig_4b_pvalue_matrix', 'Fig_4b_pvalue_matrix.csv', 'Fig 4b: Full p-value matrix (17 × 17 genes), corresponding p-values for all Pearson correlations'),
 ('Fig_4b_significant_correlations', 'Fig_4b_significant_correlations.csv', 'Fig 4b: Significant gene pairs only (p ≤ 0.05), r value shown in the matrix (17 × 17 genes)'),
 ('Fig_4c_Raw_expression_data', 'Fig_4c_Raw_expression_data.csv', 'Fig 4c: Raw 2^(-dCt) expression data (17 genes, n = 12) used as input for all four normalization methods'),
 ('Fig_4c_Hes1_normalized_data', 'Fig_4c_Hes1_normalized_data.csv', 'Fig 4c: Expression data normalized to Hes1 (16 genes, n = 12)'),
 ('Fig_4c_Kcnj11_normalized_data', 'Fig_4c_Kcnj11_normalized_data.csv', 'Fig 4c: Expression data normalized to Kcnj11 (16 genes, n = 12)'),
 ('Fig_4c_Ins1_normalized_data', 'Fig_4c_Ins1_normalized_data.csv', 'Fig 4c: Expression data normalized to Ins1 (16 genes, n = 12)'),
 ('Fig_4d_consensus_pairs', 'Fig_4d_consensus_pairs.csv', 'Fig 4d: 21 consensus gene pairs significant (p ≤ 0.05) under all four normalization methods (four-way Venn intersection); none changes sign between methods'),
 ('Fig_4d_Hes1_normalized_pairs', 'Fig_4d_Hes1_normalized_pairs.csv', 'Fig 4d: Significant pairs from Hes1-normalized method (120 of 120 pairs)'),
 ('Fig_4d_Kcnj11_normalized_pairs', 'Fig_4d_Kcnj11_normalized_pairs.csv', 'Fig 4d: Significant pairs from Kcnj11-normalized method (84 of 120 pairs)'),
 ('Fig_4d_Ins1_normalized_pairs', 'Fig_4d_Ins1_normalized_pairs.csv', 'Fig 4d: Significant pairs from Ins1-normalized method (80 of 120 pairs)'),
 ('Fig_4d_Raw_expression_pairs', 'Fig_4d_Raw_expression_pairs.csv', 'Fig 4d: Significant pairs from raw expression method (66 of 136 pairs)'),
 ('Fig_4e_Hes1_significant_correla', 'Fig_4e_Hes1_significant_correlations_only.csv', 'Fig 4e: Pearson correlation matrix (Hes1-normalized), significant pairs only (p ≤ 0.05); 21 consensus pairs annotated'),
 ('Fig_4e_consensus_annotations', 'Fig_4e_consensus_annotations.csv', 'Fig 4e: Consensus annotation key, gene pair labels, r values, consensus status, STRING classification'),
 ('Fig_6b_gene_centrality_summary', 'Fig_6b_Gene_Centrality_Summary.csv', 'Fig 6b: Pooled consensus network of 13 genes (21 edges), node centrality summary (degree, betweenness, closeness, eigenvector) and hub classification (degree ≥ 5, set post hoc after Hey1 was excluded)'),
 ('Fig_6b_edge_correlation_statist', 'Fig_6b_Edge_Correlation_Statistics.csv', 'Fig 6b: Pooled consensus network, edge list with Pearson r, |r|, edge type (positive/negative) and edge width (21 edges)'),
]
for name, csv, desc in sets:
    h, rows = csv_rows(R17 + '\\' + csv)
    if name == 'Fig_4c_normalization_comparison': pass
    rebuild(name, desc, src(csv), h, rows, key_width=30 if name.startswith(('Fig_4b', 'Fig_4e', 'Fig_6b')) else None)
h, rows = csv_rows(R17 + r'\Fig_4c_normalization_comparison_data.csv')
rows = [[(r[0].replace('Normalized', 'normalized') if isinstance(r[0], str) else r[0])] + r[1:] for r in rows]
rebuild('Fig_4c_normalization_comparison', 'Fig 4c: Significant (p ≤ 0.05) / Strong (|r| ≥ 0.7) / Very Strong (|r| ≥ 0.8) gene pair counts across four normalization methods (17 genes: 136 pairs raw, 120 pairs per normalized scheme)', src('Fig_4c_normalization_comparison_data.csv'), h, rows)

# ------------------------------------------------------------------ Fig 5 (category names as in the manuscript)
h, rows = csv_rows(DX + r'\Fig5C_All_Pairs_Classified.csv')
cat = {r[0]: REL5[r[4]] for r in rows}
rows5 = [[r[0], r[1], r[2], r[3], REL5[r[4]] in (REL5['Known Validated'], 'STRING Only'), REL5[r[4]] in (REL5['Known Validated'], REL5['Novel Discovery']), REL5[r[4]]] for r in rows]
rebuild('Fig_5_All_Classified_Pairs', 'Fig 5: Complete classified gene pair list, all 120 pairs (every pair has |r| ≥ 0.7 under Hes1 normalization) with STRING annotation, category (Known and analytically consistent / Analytically consistent, database-unregistered pairs / STRING Only / Not analytically consistent) and r value',
        src('Fig5C_All_Pairs_Classified.csv') + '; category names as used in the manuscript', ['Pair', 'Gene1', 'Gene2', 'r', 'In_STRING', 'Four_Method_Consistent', 'Category'], rows5)
h, rows = csv_rows(DX + r'\Fig5_Summary_Statistics.csv'); rows = [[REL5[r[0]]] + r[1:] for r in rows]
rebuild('Fig_5_Summary_Statistics', 'Fig 5: Summary statistics, pair counts and |r| summary by interaction category', src('Fig5_Summary_Statistics.csv') + '; category names as used in the manuscript', h, rows)
h, rows = csv_rows(DX + r'\Fig5B_Coverage_Summary.csv'); rows = [[REL5[r[0]]] + r[1:] for r in rows]
rebuild('Fig_5b_Coverage_Summary', 'Fig 5b: Network coverage summary, proportion of pairs in each category (120 pairs)', src('Fig5B_Coverage_Summary.csv') + '; category names as used in the manuscript', h, rows)
h, rows = csv_rows(DX + r'\Fig5D_Novel_Discoveries.csv')
rows5d = [[f'{r[0]}-{r[1]}', r[0], r[1], r[2], False, True, REL5['Novel Discovery'], f'{r[0]}-{r[1]}', i + 1] for i, r in enumerate(rows)]
rebuild('Fig_5d_Top_Novel_Discoveries', 'Fig 5d: The 15 strongest of the 16 analytically consistent, database-unregistered pairs (absent from STRING), ranked by |r| (r ≥ 0.93)',
        src('Fig5D_Novel_Discoveries.csv') + '; pair labels, category and rank added', ['Pair', 'Gene1', 'Gene2', 'r', 'In_STRING', 'Four_Method_Consistent', 'Category', 'Pair_Label', 'Rank'], rows5d)

# ------------------------------------------------------------------ Fig 6f (curated pathway-gene map: "mapped genes", not the Fig 6b hubs)
h, rows = csv_rows(R17 + r'\Supplementary_Table_6Layer_Complete_Mapping.csv')
h = ['Approval.Status', 'Hub.Drugs', 'Action.Type', 'Mapped.Genes', 'Pathways', 'Annotated.effects', 'Flow']
rebuild('Fig_6f_six_layer_complete_mappi', 'Fig 6f: Six-layer mapping, approval status → hub drug → action type → mapped gene → pathway → annotated effect (flow = frequency × connection strength × approval confidence). The mapped genes come from a curated pathway-gene table and are not the Fig. 6b hubs; the annotated effects are literature annotations',
        src('Supplementary_Table_6Layer_Complete_Mapping.csv') + '; column "Clinical Outcomes" renamed "Annotated effects", "Hub Genes" renamed "Mapped genes"', h, rows)
wsx = wb['Fig_6f_pathway_gene_connectivit']; wsx['C5'] = 'Mapped_Gene'
set_desc('Fig_6f_pathway_gene_connectivit', 'Fig 6f: Curated pathway-gene table used for the "Mapped genes" layer (not the Fig. 6b hubs): pathway, gene and connection strength',
         f'Source: Fig_6f_pathway_gene_connectivity.csv  |  Generated: 2026-04-23; column "Hub_Gene" renamed "Mapped_Gene" {DATE}')

# ------------------------------------------------------------------ Fig S1, all-pairs statistics, null counts (Python pipeline)
genes = P.GENES17; df = P.load(genes); sch = P.schemes(df); res = P.sig_sets(sch); adj = P.adjusted(sch)
order = ['Raw', 'Hes1-normalized', 'Kcnj11-normalized', 'Ins1-normalized']; lab = dict(zip(order, ['Raw', 'Hes1', 'Kcnj11', 'Ins1']))
oh_, orows = old_rows('R_all_pairs_statistics', 4); old_labels = list(dict.fromkeys(r[0] for r in orows)); print('scheme labels in the old sheet:', old_labels)
lab = {'Raw': old_labels[0]}
for k in order[1:]:
    for o in old_labels:
        if o.lower().startswith(k.split('-')[0].lower()): lab[k] = o
pair_rows = []
for k in order:
    m = sch[k]; g = list(m.columns); r_, p_ = P.corr_p(m.values)
    items = [(g[i], g[j], r_[i, j], p_[i, j]) for i, j in P.pairs(g)]
    q_ = multipletests([it[3] for it in items], method='fdr_bh')[1]
    for it, qq in zip(items, q_):
        a_ = adj[k][tuple(sorted((it[0], it[1])))]
        pair_rows.append([lab[k], it[0], it[1], float(it[2]), float(it[3]), float(qq), float(a_[0]), float(a_[1]), float(a_[2])])
chk = {k: sum(1 for r in pair_rows if r[0] == lab[k] and r[4] < 0.05) for k in order}; print('pooled p<0.05 per scheme:', chk)
assert chk == {'Raw': 66, 'Hes1-normalized': 120, 'Kcnj11-normalized': 84, 'Ins1-normalized': 80}
assert len(pair_rows) == 136 + 3 * 120
rebuild('R_all_pairs_statistics', 'All gene pairs (17 genes) under each normalization scheme: pooled (df 10) and treatment-adjusted (df 8) Pearson r, p and Benjamini-Hochberg q',
        'Source: Analysis_outputs/reanalysis_17genes/pipeline.py (Python re-implementation of the pooled pipeline; reproduces every count of the 20-gene analysis)', oh_, pair_rows, hdr=4)
npz = np.load(RE + r'\perm_arrays_17.npz'); tc = json.load(open(RE + r'\table_counts.json'))['17']
obs = dict(zip(['Raw', 'Hes1-normalized', 'Kcnj11-normalized', 'Ins1-normalized', 'cons', 'cons_sign'], [66, 120, 84, 80, 21, 21]))
names = {'Raw': 'Raw expression', 'Hes1-normalized': 'Hes1-normalized', 'Kcnj11-normalized': 'Kcnj11-normalized', 'Ins1-normalized': 'Ins1-normalized', 'cons': 'Significant under all four schemes', 'cons_sign': 'Same sign under all four schemes'}
s1a = []
for k, o in obs.items():
    a = npz[k].astype(float); s1a.append([names[k], o, round(float(a.mean()), 2), round(float(a.std(ddof=1)), 2), float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5)), round(float((a >= o).mean()), 4)])
rebuild('Fig_S1a_Observed_vs_Null', "Fig S1a: Observed pooled significant-pair counts vs a group-preserving permutation null (17 genes; 2000 permutations; each gene's group means kept, log-scale within-group residuals permuted within group; schemes rebuilt from each permuted matrix)",
        'Source: Analysis_outputs/reanalysis_17genes/perm_arrays.py (seed 2026); per-permutation counts in sheet R_null_counts', ['Scheme', 'Observed', 'Null_mean', 'Null_SD', 'Null_2.5pct', 'Null_97.5pct', 'P_null_ge_observed'], s1a, hdr=4)
s1b = []
for k in order:
    t = tc[k]; s1b.append([names[k], t['pairs'], t['pooled_p'], t['pooled_q'], t['adj_p'], t['adj_q']])
t = tc['all']; s1b.append(['Significant under all four schemes', '-', t['pooled_p'], t['pooled_q'], t['adj_p'], t['adj_q']])
rebuild('Fig_S1b_adjusted_counts', 'Fig S1b: Significant pairs before and after adjustment for treatment (17 genes; partial correlation, df = 8; Benjamini-Hochberg within each scheme)',
        'Source: Analysis_outputs/reanalysis_17genes/table_counts.py; per-pair values in sheet R_all_pairs_statistics', ['Scheme', 'Pairs_tested', 'Pooled_p_lt_0.05', 'Pooled_BH_q_lt_0.05', 'Adjusted_p_lt_0.05', 'Adjusted_BH_q_lt_0.05'], s1b, hdr=4)
oh_n, _ = old_rows('R_null_counts', 4)
nulls = [[i + 1] + [int(npz[k][i]) for k in ['Raw', 'Hes1-normalized', 'Kcnj11-normalized', 'Ins1-normalized', 'cons', 'cons_sign']] for i in range(len(npz['Raw']))]
rebuild('R_null_counts', 'Per-permutation counts of pooled significant pairs (p < 0.05) under the group-preserving null (17 genes; 2000 permutations, seed 2026)',
        'Source: Analysis_outputs/reanalysis_17genes/perm_arrays.py', oh_n, nulls, hdr=4)
del wb['Fig_S1c_Hey1_signflip']                                  # the six sign-flipping Hey1 pairs no longer exist

# ------------------------------------------------------------------ Fig 3g sanity check: the PCA input is the 10 unaffected genes
x = pd.read_excel(R17 + r'\Fig_3f_PCA_Relative_mRNA_expression.xlsx'); X = x.drop(columns=['Group']).values.astype(float); Z = (X - X.mean(0)) / X.std(0, ddof=1)
sd = np.linalg.svd(Z, compute_uv=False) / math.sqrt(len(Z) - 1)
_, ov = old_rows('Fig_3g_PCA_variance_explained'); print('PCA check: PC1 sd now %.4f, sheet %.4f' % (sd[0], ov[0][1]))
assert abs(sd[0] - ov[0][1]) < 1e-3

# ------------------------------------------------------------------ Gapdh-reference sensitivity (26 August run), exact recalculation from the original calculation sheet
gs = json.load(open(RE + r'\gapdh_exact.json')); gp = pd.read_csv(RE + r'\gapdh_exact_pairs.csv')
V = [('reported', 'Reported'), ('C0', 'Mean_of_two_valid_wells'), ('C1', 'Mean_of_all_four_wells')]
c0 = {tuple(x) for x in gs['reported']['consensus']}
sens = []
for g in ['Tcf7l2', 'Kcnj11', 'Cacna1d', 'Cacna1c', 'Ptbp1']:
    for key, lab_ in [('fc_dapt', 'DAPT fold change vs CTRL'), ('CTRL-DAPT', 'Adjusted p, CTRL vs DAPT'), ('DAPT-DKK1', 'Adjusted p, DAPT vs DKK-1')]:
        sens.append([g, lab_] + [round(gs[k]['per_gene'][g][key], 10 if key != 'fc_dapt' else 3) for k, _ in V])
for sch_ in ['Raw', 'Ins1-normalized', 'Kcnj11-normalized', 'Hes1-normalized']:
    sens.append(['Pooled network', f'Significant pairs, {sch_} (p < 0.05)'] + [gs[k]['sig'][sch_] for k, _ in V])
sens.append(['Pooled network', 'Consensus pairs (significant in all four schemes)'] + [len(gs[k]['consensus']) for k, _ in V])
sens.append(['Pooled network', 'Of the 21 reported consensus pairs: still consensus'] + [len(c0 & {tuple(x) for x in gs[k]['consensus']}) for k, _ in V])
sens.append(['Pooled network', 'Of the 21 reported consensus pairs: lost'] + [len(c0 - {tuple(x) for x in gs[k]['consensus']}) for k, _ in V])
sens.append(['Pooled network', 'New consensus pairs'] + [len({tuple(x) for x in gs[k]['consensus']} - c0) for k, _ in V])
for sch_ in ['Raw', 'Ins1-normalized', 'Kcnj11-normalized', 'Hes1-normalized']:
    sens.append(['Treatment-adjusted', f'Pairs with BH q < 0.05, {sch_}'] + [gs[k]['adj_q05'][sch_] for k, _ in V])
for g in ['Glp1r', 'Cacna1d', 'Ptbp1', 'Tcf7l2', 'Wnt5a']:
    sens.append(['Consensus network', f'Degree of {g}'] + [gs[k]['degree'].get(g, 0) for k, _ in V])
sens.append(['Consensus network', 'Hubs (degree >= 5)'] + [', '.join(g for g, v in gs[k]['degree'].items() if v >= 5) for k, _ in V])
sens.append(['Cacna1d-Glp1r', 'Consensus pair (significant in all four schemes)'] + [bool(['Cacna1d', 'Glp1r'] in [sorted(x) for x in gs[k]['consensus']]) for k, _ in V])
for sch_ in ['Hes1-normalized', 'Raw', 'Ins1-normalized', 'Kcnj11-normalized']:
    sens.append(['Cacna1d-Glp1r', f'Pearson r, {sch_}' + (' (input of the MLSS Synergy Bonus)' if sch_ == 'Hes1-normalized' else '')] + [round(gs[k]['r_cacna1d_glp1r'][sch_], 3) for k, _ in V])
desc_g = ('Sensitivity of the 17-gene results to the Gapdh reference of the 26 August 2020 DAPT group. Tcf7l2, Kcnj11, Cacna1d, Cacna1c and Ptbp1 were measured only in that run. Reported = values as analyzed: in the original calculation sheet the Gapdh Cq of wells C10 and D10 '
          '(22.57, 22.80) were entered twice for the four DAPT replicates and those of A10 and B10 (24.99, 24.13) were set aside. Mean_of_two_valid_wells = reference 22.68 for all four replicates. Mean_of_all_four_wells = reference 23.62 (the two wells with low end-point fluorescence included). '
          'Control and DKK-1 values unchanged. Fold change = mean 2^-dCt of the DAPT group / control group; p values: one-way ANOVA with Bonferroni-corrected pooled-variance t tests. Permutation null and MLSS not recomputed.')
for nm in ('Gapdh_sensitivity', 'Gapdh_sensitivity_pairs'): wb.create_sheet(title=nm, index=wb.sheetnames.index('R_all_pairs_statistics'))
rebuild('Gapdh_sensitivity', desc_g, 'Source: Analysis_outputs/reanalysis_17genes/gapdh_exact.py (Python pipeline.py) from Undergrad_Glucose uptake signaling/RT-PCR/2020-08-26_... Quantification Summary.xlsx', ['Section', 'Quantity', 'Reported', 'Mean_of_two_valid_wells', 'Mean_of_all_four_wells'], sens, key_width=20)
wb['Gapdh_sensitivity'].column_dimensions['C'].width = 56
for col in 'DEF': wb['Gapdh_sensitivity'].column_dimensions[col].width = 24
gp_rows = [[r.pair, bool(r.involves_G5), bool(r.reported), bool(r.C0), bool(r.C1)] for r in gp.itertuples()]
rebuild('Gapdh_sensitivity_pairs', 'Consensus pairs (significant in all four normalization schemes) under the reported values and under the two recalculations of the 26 August DAPT reference (see sheet Gapdh_sensitivity); Involves_G5_gene = pair contains Tcf7l2, Kcnj11, Cacna1d, Cacna1c or Ptbp1',
        'Source: Analysis_outputs/reanalysis_17genes/gapdh_exact_pairs.csv', ['Pair', 'Involves_G5_gene', 'Consensus_reported', 'Consensus_two_valid_wells', 'Consensus_all_four_wells'], gp_rows)

# ------------------------------------------------------------------ Index
ws = wb['Index']
ws['B2'] = f'Generated: {DATE}  |  Sheets: {len(wb.sheetnames) - 1}  |  17-gene analysis (Hey1, Ryr1, Wnt2 and Wnt2b excluded)'
rowmap = {}
for r in range(5, ws.max_row + 1):
    n = ws.cell(r, 2).value
    if n: rowmap[n] = r
r_s1c = rowmap.get('Fig_S1c_Hey1_signflip')
if r_s1c: ws.delete_rows(r_s1c)
rowmap = {ws.cell(r, 2).value: r for r in range(5, ws.max_row + 1) if ws.cell(r, 2).value}
changed = {n: (nr, nc) for n, nr, nc in log}
for n, (nr, nc) in changed.items():
    r = rowmap.get(n)
    if r: ws.cell(r, 4).value = nr; ws.cell(r, 5).value = nc; ws.cell(r, 6).value = wb[n]['B2'].value
r = rowmap['Fig_3e_stats']; ws.cell(r, 6).value = wb['Fig_3e_stats']['B2'].value
r = rowmap['Fig_6f_pathway_gene_connectivit']; ws.cell(r, 6).value = wb['Fig_6f_pathway_gene_connectivit']['B2'].value
for n in ('Fig_1d_stats',):
    ws.cell(rowmap[n], 3).value = 'Fig_1d_stats.csv'
r_ins = rowmap['R_all_pairs_statistics']; ws.insert_rows(r_ins, 2)
for k, nm in enumerate(('Gapdh_sensitivity', 'Gapdh_sensitivity_pairs')):
    for c in range(1, 7): ws.cell(r_ins + k, c)._style = copy(ws.cell(r_ins - 1, c)._style)
    ws.cell(r_ins + k, 2).value = nm; ws.cell(r_ins + k, 3).value = '(workbook sheet)'; ws.cell(r_ins + k, 4).value = changed[nm][0]; ws.cell(r_ins + k, 5).value = changed[nm][1]
    ws.cell(r_ins + k, 6).value = wb[nm]['B2'].value
rowmap = {ws.cell(r, 2).value: r for r in range(5, ws.max_row + 1) if ws.cell(r, 2).value}
w5 = wb['Fig_5d_Top_Novel_Discoveries']; w5.title = 'Fig_5d_Top_unregistered_pairs'; w5['B1'] = w5.title      # the manuscript no longer calls these pairs novel
ws.cell(rowmap['Fig_5d_Top_Novel_Discoveries'], 2).value = w5.title; ws.cell(rowmap['Fig_5d_Top_Novel_Discoveries'], 3).value = 'Fig5D_Novel_Discoveries.csv'
wb.save(OUT)
print('saved', OUT, '| sheets', len(wb.sheetnames))
for n, nr, nc in log: print(f'  rebuilt {n:36s} {nr:4d} rows x {nc:2d} cols')
