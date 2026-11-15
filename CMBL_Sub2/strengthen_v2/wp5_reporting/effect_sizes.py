"""WP5a: effect sizes with 95% intervals, across-gene Benjamini-Hochberg q, and the 26 August Gapdh sensitivity for the INS-1 RT-qPCR panel,
plus ratios for DNA content, 2-NBDG and C-peptide per ug DNA where the data support them.

Inputs (read only; project root = three folders above this script, override with the CMBL_ROOT environment variable):
  Submission_package_JTM/Additional file 3.xlsx   reported means, SEMs, replicate fold changes, pairwise p (17 genes)
  Submission_package_JTM/Additional file 4.xlsx   full-precision 2^-dCt matrix, DNA, 2-NBDG and ELISA readings, reported statistics
  Analysis_outputs/reanalysis_17genes/gapdh_exact.json   Gapdh values of the 26 August DAPT group (reported, C0, C1)
  SBA/Figures/Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx   optional cross-check of the matrix
Outputs next to this script: effect_sizes.csv, effect_sizes.md, gate_reproduction.json.
Usage: python -I effect_sizes.py
Gate (written to gate_reproduction.json): the printed means, SEMs, replicate values and pairwise and overall p of Additional file 3, the DNA, 2-NBDG and
C-peptide statistics of Additional file 4 and the Gapdh variants of gapdh_exact.json / Table S1c are recomputed first and compared with what the files print.

Method in one paragraph. Replicate fold change = 2^-dCt of the replicate / arithmetic mean of the four CTRL values of the gene (this is what Additional
file 3 lists as R1-R4). Each replicate fold change is log2-transformed. For treated vs CTRL the effect is the difference of the group means of the log2 values
(= log2 of the ratio of geometric means); the primary interval is the Welch t-interval (n = 4 per group). Sensitivity: pooled-variance interval over all three
groups (df = 9, the model behind the manuscript's pairwise tests). The manuscript's fold change is the ratio of arithmetic means and is listed beside it.
BH q = Benjamini-Hochberg adjustment of the 17 raw pairwise p of one contrast (pooled-variance t-test after ANOVA, Wilcoxon rank-sum for Wnt5a and Lef1),
recomputed at full precision and checked against the printed p. All intervals are within-experiment: one preparation, wells are not independent biological replicates.
Seeds: nothing random is used."""
import os, sys, json, csv, math
from pathlib import Path
import numpy as np
from scipy import stats
import openpyxl
from statsmodels.stats.multitest import multipletests

sys.stdout.reconfigure(encoding='utf8')
HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get('CMBL_ROOT', str(HERE.parents[2])))
PKG = ROOT / 'Submission_package_JTM'
AF3 = PKG / 'Additional file 3.xlsx'
AF4 = PKG / 'Additional file 4.xlsx'
GAPDH_JSON = ROOT / 'CMBL_Sub2' / 'reanalysis_17genes' / 'gapdh_exact.json'
FIG4B = ROOT.parent / 'SBA' / 'Figures' / 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx'
LABEL = 'within-experiment, one preparation, wells are not independent biological replicates'
G5 = ['Tcf7l2', 'Kcnj11', 'Cacna1d', 'Cacna1c', 'Ptbp1']      # DAPT values depend on the 26 August Gapdh wells
CONTRASTS = ['DAPT', 'DKK-1']
GI = {'CTRL': slice(0, 4), 'DAPT': slice(4, 8), 'DKK-1': slice(8, 12)}
GATE = {}


# ---------------------------------------------------------------- loaders
def strip_left(row):
    row = list(row)
    while row and row[0] is None:
        row.pop(0)
    return row


def load_af3():
    wb = openpyxl.load_workbook(AF3, read_only=True, data_only=True)
    rows = list(wb['Analyzed genes (17)'].iter_rows(values_only=True))
    wb.close()
    h = next(i for i, r in enumerate(rows) if r and r[0] == 'Gene')
    out = []
    for r in rows[h + 1:]:
        if r[0] is None or str(r[0]).startswith('Abbreviations'):
            break
        out.append(dict(gene=r[0], category=r[1], test=r[2], overall_p=str(r[3]),
                        mean={'CTRL': r[4], 'DAPT': r[10], 'DKK-1': r[16]}, sem={'CTRL': r[5], 'DAPT': r[11], 'DKK-1': r[17]},
                        rep={'CTRL': list(r[6:10]), 'DAPT': list(r[12:16]), 'DKK-1': list(r[18:22])},
                        p_raw={'DAPT': r[22], 'DKK-1': r[26], 'D-K': r[30]}, p_adj={'DAPT': r[23], 'DKK-1': r[27], 'D-K': r[31]}))
    return out


def af4_rows(wb, name):
    return [strip_left(r) for r in wb[name].iter_rows(values_only=True)]


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def stat_means(rows, start_label=None):
    """descriptive block of an Additional file 4 'stats' sheet (numbers are stored as text there): first row per group -> mean, SEM, n;
    with start_label only rows after the first row starting with it (second block of Fig_3b_stats)"""
    out = {}
    active = start_label is None
    for r in rows:
        if r and start_label and isinstance(r[0], str) and r[0].startswith(start_label):
            active = True
        if active and r and r[0] in ('CTRL', 'DAPT', 'DKK-1') and len(r) >= 4 and num(r[1]) is not None and num(r[2]) is not None and num(r[3]) is not None and r[0] not in out:
            out[r[0]] = dict(mean=float(r[1]), sem=float(r[2]), n=int(float(r[3])))
    return out


def stat_pairwise(rows, start_label=None):
    """CTRL vs DAPT and CTRL vs DKK-1 rows (p_raw, p_adjusted) of a stats sheet; with start_label only rows after the first row starting with it"""
    out = {}
    active = start_label is None
    for r in rows:
        if r and start_label and isinstance(r[0], str) and r[0].startswith(start_label):
            active = True
        if active and r and r[0] == 'CTRL' and len(r) >= 4 and r[1] in ('DAPT', 'DKK-1') and num(r[2]) is not None and num(r[3]) is not None and r[1] not in out:
            out[r[1]] = (float(r[2]), float(r[3]))
    return out


# ---------------------------------------------------------------- statistics
def sem(a):
    a = np.asarray(a, float)
    return float(a.std(ddof=1) / math.sqrt(len(a)))


def welch(a, b, alpha=0.05):
    """treated a minus control b on the log scale: difference, 95% t-interval (Welch-Satterthwaite df), two-sided p"""
    a, b = np.asarray(a, float), np.asarray(b, float)
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    se = math.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    d = float(a.mean() - b.mean())
    t = stats.t.ppf(1 - alpha / 2, df)
    p = float(2 * stats.t.sf(abs(d / se), df))
    return d, d - t * se, d + t * se, float(df), p


def pooled3(groups, i, j, alpha=0.05):
    """difference of means groups[j] - groups[i] with the within-group variance pooled over all groups (df = N - k)"""
    n = [len(g) for g in groups]
    ss = sum(((np.asarray(g) - np.mean(g)) ** 2).sum() for g in groups)
    df = sum(n) - len(groups)
    se = math.sqrt(ss / df * (1 / n[i] + 1 / n[j]))
    d = float(np.mean(groups[j]) - np.mean(groups[i]))
    t = stats.t.ppf(1 - alpha / 2, df)
    return d, d - t * se, d + t * se, df, float(2 * stats.t.sf(abs(d / se), df))


def pairwise_p_anova(groups, i, j):
    """raw p of the manuscript's pairwise test after one-way ANOVA: pooled-variance t-test, MSE over the three groups, df = 9"""
    return pooled3(groups, i, j)[4]


def pairwise_p_wilcox(a, b):
    """Wilcoxon rank-sum, two-sided; exact when there are no ties, normal approximation with continuity correction otherwise (as R wilcox.test)"""
    return float(stats.mannwhitneyu(a, b, alternative='two-sided', method='auto', use_continuity=True).pvalue)


def printed_match(val, s):
    """True if val agrees with the printed string s to the printed precision"""
    s = str(s).strip()
    if s.startswith('<'):
        return val < float(s[1:])
    f = float(s)
    if 'e' in s.lower():
        mant, ex = s.lower().split('e')
        nd = len(mant.split('.')[1]) if '.' in mant else 0
        return abs(val - f) <= 0.5 * 10 ** (int(ex) - nd) + 1e-300
    nd = len(s.split('.')[1]) if '.' in s else 0
    return abs(val - f) <= 0.5 * 10 ** (-nd) + 1e-12


def bh(p):
    return multipletests(list(p), method='fdr_bh')[1]


def fmt(x, d=2):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ''
    return f'{x:.{d}f}'


def fmtp(p):
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return ''
    if p >= 0.001:
        return f'{p:.3f}'
    return f'{p:.2e}'


# ---------------------------------------------------------------- RT-qPCR panel
def rt_panel(af3):
    wb = openpyxl.load_workbook(AF4, read_only=True, data_only=True)
    rows = af4_rows(wb, 'Fig_4c_Raw_expression_data')
    wb.close()
    h = next(i for i, r in enumerate(rows) if r and r[0] == 'Wnt5a')
    genes = [g for g in rows[h]]
    X = np.array([[float(c) for c in r] for r in rows[h + 1:h + 13]])
    assert X.shape == (12, 17), X.shape
    assert [r['gene'] for r in af3] != [] and set(genes) == {r['gene'] for r in af3}
    # optional cross-check against the matrix in SBA/Figures
    if FIG4B.exists():
        import pandas as pd
        m = pd.read_excel(FIG4B)[genes].astype(float).values
        GATE['matrix_AF4_vs_Fig4b_max_abs_rel_diff'] = float(np.max(np.abs(m - X) / np.abs(X)))
    else:
        GATE['matrix_AF4_vs_Fig4b_max_abs_rel_diff'] = None
    res = {}
    over = []
    for g in genes:
        x = X[:, genes.index(g)]
        fc = x / x[GI['CTRL']].mean()
        res[g] = dict(x=x, fc=fc, groups=[fc[GI['CTRL']], fc[GI['DAPT']], fc[GI['DKK-1']]])
        c = x[GI['CTRL']]
        over.append(float(c.mean() / np.exp(np.log(c).mean())))      # CTRL arithmetic mean if the CTRL mean dCt (geometric mean of 2^-dCt) were the reference
    # Methods say "relative to the mean dCt of all CTRL replicates"; the CTRL mean of Additional file 3 is exactly 1, which fits the arithmetic mean of 2^-dCt
    GATE['ctrl_mean_if_normalized_to_mean_dCt_max'] = float(max(over))
    GATE['ctrl_mean_if_normalized_to_mean_dCt_median'] = float(np.median(over))
    return genes, res


def raw_p(groups, i, j, use_kw):
    """raw p of the manuscript's pairwise test between groups i and j"""
    return pairwise_p_wilcox(groups[j], groups[i]) if use_kw else pairwise_p_anova(groups, i, j)


def reproduce_reported(af3, res):
    """gate: Additional file 3 (means, SEMs, replicates, overall and pairwise p) against two recomputations:
    'printed' = from the printed 3-decimal replicate values R1-R4, 'full' = from the 12 full-precision values of Additional file 4"""
    names = ['CTRL', 'DAPT', 'DKK-1']
    G = {}
    for basis in ('printed', 'full'):
        d_mean, d_sem, d_rep, d_ctrl1 = [], [], [], []
        p_ok = p_n = ov_ok = ov_n = 0
        p_bad, ov_bad = [], []
        for r in af3:
            g = r['gene']
            groups = [np.array(r['rep'][k], float) for k in names] if basis == 'printed' else res[g]['groups']
            for k, name in enumerate(names):
                a = groups[k]
                d_mean.append(abs(a.mean() - r['mean'][name]))
                d_sem.append(abs(sem(a) - r['sem'][name]))
                d_rep.extend(np.abs(a - np.array(r['rep'][name], float)))
            d_ctrl1.append(abs(groups[0].mean() - 1))
            use_kw = r['test'].startswith('Kruskal')
            for (i, j), key in [((0, 1), 'DAPT'), ((0, 2), 'DKK-1'), ((1, 2), 'D-K')]:
                praw = raw_p(groups, i, j, use_kw)
                for val, s in ((praw, r['p_raw'][key]), (min(1.0, 3 * praw), r['p_adj'][key])):
                    p_n += 1
                    if printed_match(val, s):
                        p_ok += 1
                    else:
                        p_bad.append((g, key, float(val), str(s)))
            ov = float(stats.kruskal(*groups).pvalue) if use_kw else float(stats.f_oneway(*groups).pvalue)
            ov_n += 1
            if printed_match(ov, r['overall_p']):
                ov_ok += 1
            else:
                ov_bad.append((g, float(ov), str(r['overall_p'])))
        G[basis] = dict(group_means_max_abs_diff=float(max(d_mean)), group_sems_max_abs_diff=float(max(d_sem)), n_group_means=len(d_mean),
                        replicates_max_abs_diff=float(max(d_rep)), n_replicates=len(d_rep), ctrl_mean_is_1_max_abs_dev=float(max(d_ctrl1)),
                        pairwise_p_checked=p_n, pairwise_p_matching_printed=p_ok,
                        pairwise_p_not_matching=p_bad, overall_p_checked=ov_n, overall_p_matching_printed=ov_ok, overall_p_not_matching=ov_bad)
    GATE['n_genes'] = len(af3)
    GATE['from_printed_replicates'] = G['printed']
    GATE['from_full_precision_values'] = G['full']


def gene_rows(af3, genes, res):
    """one record per gene and contrast (reported analysis)"""
    # raw p of the manuscript's pairwise tests, per contrast: from the full-precision values (used) and from the printed 3-decimal replicates (check)
    names = ['CTRL', 'DAPT', 'DKK-1']
    ci_diff = []
    for r in af3:
        R = res[r['gene']]
        use_kw = r['test'].startswith('Kruskal')
        printed = [np.array(r['rep'][k], float) for k in names]
        R['praw'], R['praw_printed'] = {}, {}
        for (i, j), key in [((0, 1), 'DAPT'), ((0, 2), 'DKK-1')]:
            R['praw'][key] = raw_p(R['groups'], i, j, use_kw)
            R['praw_printed'][key] = raw_p(printed, i, j, use_kw)
            a, b = welch(np.log2(R['groups'][j]), np.log2(R['groups'][i])), welch(np.log2(printed[j]), np.log2(printed[i]))
            ci_diff.append(max(abs(a[0] - b[0]), abs(a[1] - b[1]), abs(a[2] - b[2])))
    GATE['ci_limits_full_vs_printed_replicates_max_abs_diff_log2'] = float(max(ci_diff))
    rel = [abs(res[r['gene']]['praw'][k] - res[r['gene']]['praw_printed'][k]) / res[r['gene']]['praw_printed'][k]
           for r in af3 for k in CONTRASTS if res[r['gene']]['praw_printed'][k] > 0]
    GATE['p_raw_full_vs_printed_replicates_max_rel_diff'] = float(max(rel))
    GATE['p_raw_full_vs_printed_replicates_median_rel_diff'] = float(np.median(rel))
    q = {}
    for key in CONTRASTS:
        qq = bh([res[r['gene']]['praw'][key] for r in af3])
        q[key] = {r['gene']: float(v) for r, v in zip(af3, qq)}
        # the same BH on the p of the printed replicates (the manuscript's p) and on the printed p strings
        qp = bh([res[r['gene']]['praw_printed'][key] for r in af3])
        qs = bh([float(r['p_raw'][key]) for r in af3])
        GATE['bh_q_max_abs_diff_full_vs_printed_replicates_' + key] = float(np.max(np.abs(np.array(list(q[key].values())) - qp)))
        GATE['bh_q_max_abs_diff_full_vs_printed_strings_' + key] = float(np.max(np.abs(np.array(list(q[key].values())) - qs)))
        GATE['bh_q05_status_differs_' + key] = [r['gene'] for r, a, b in zip(af3, q[key].values(), qp) if (a < 0.05) != (b < 0.05)]
    recs = []
    for r in af3:
        g = r['gene']
        R = res[g]
        L = [np.log2(a) for a in R['groups']]
        for key, k in (('DAPT', 1), ('DKK-1', 2)):
            d, lo, hi, df, pw = welch(L[k], L[0])
            d3, lo3, hi3, df3, _ = pooled3(L, 0, k)
            recs.append(dict(readout='RT-qPCR', feature=g, category=r['category'], contrast=key + ' vs CTRL', variant='as reported',
                             flag_26Aug_Gapdh=(key == 'DAPT' and g in G5), test_in_manuscript=r['test'], n_ctrl=4, n_trt=4,
                             ctrl_mean_reported=r['mean']['CTRL'], ctrl_sem_reported=r['sem']['CTRL'],
                             trt_mean_reported=r['mean'][key], trt_sem_reported=r['sem'][key],
                             trt_mean_recomputed=float(R['groups'][k].mean()), trt_sem_recomputed=sem(R['groups'][k]),
                             fc_ratio_arith=float(R['groups'][k].mean() / R['groups'][0].mean()),
                             log2fc_geo=d, ci95_lo=lo, ci95_hi=hi, ci_method='Welch t, log2 scale', df=df,
                             fc_geo=2 ** d, fc_geo_lo=2 ** lo, fc_geo_hi=2 ** hi,
                             ci95_lo_pooled3=lo3, ci95_hi_pooled3=hi3, df_pooled3=df3,
                             ci_excludes_0=(lo > 0 or hi < 0), ci_excludes_0_pooled3=(lo3 > 0 or hi3 < 0), p_welch_log2=pw,
                             p_raw_manuscript=str(r['p_raw'][key]), p_adj_bonf_manuscript=str(r['p_adj'][key]),
                             p_raw_recomputed=R['praw'][key], p_adj_bonf_recomputed=min(1.0, 3 * R['praw'][key]), bh_q_17genes=q[key][g], interval_label=LABEL, note=''))
    return recs, q


def variants_g5(af3, res, q_reported):
    """DAPT values of the five 26 August genes under the reported, C0 and C1 Gapdh references"""
    J = json.load(open(GAPDH_JSON, encoding='utf8'))
    HK = np.array(J['hk_used'], float)
    REF = {'reported': None, 'C0': float(J['ref']['C0']), 'C1': float(J['ref']['C1'])}
    out = []
    chk = []
    for var in ('reported', 'C0', 'C1'):
        praw_all = {r['gene']: res[r['gene']]['praw']['DAPT'] for r in af3}     # start from the reported p
        recs = []
        for g in G5:
            R = res[g]
            x = R['x'].copy()
            if REF[var] is not None:
                x[GI['DAPT']] = x[GI['DAPT']] * 2.0 ** (REF[var] - HK)       # 2^-(Cq - ref) = 2^-(Cq - HK) * 2^(ref - HK)
            fc = x / x[GI['CTRL']].mean()
            groups = [fc[GI['CTRL']], fc[GI['DAPT']], fc[GI['DKK-1']]]
            L = [np.log2(a) for a in groups]
            d, lo, hi, df, pw = welch(L[1], L[0])
            d3, lo3, hi3, df3, _ = pooled3(L, 0, 1)
            praw = pairwise_p_anova(groups, 0, 1)
            praw_all[g] = praw
            recs.append(dict(gene=g, variant=var, mean=float(groups[1].mean()), sem=sem(groups[1]), fc_arith=float(groups[1].mean() / groups[0].mean()),
                             log2fc=d, lo=lo, hi=hi, df=df, lo3=lo3, hi3=hi3, p_raw=praw, p_bonf=min(1.0, 3 * praw)))
            # check against gapdh_exact.json (fold change, Bonferroni p CTRL-DAPT)
            ref = J[var]['per_gene'][g]
            chk.append((abs(recs[-1]['fc_arith'] - ref['fc_dapt']), abs(recs[-1]['p_bonf'] - ref['CTRL-DAPT']) / max(ref['CTRL-DAPT'], 1e-300)))
        qv = dict(zip(praw_all.keys(), bh(list(praw_all.values()))))
        for rc in recs:
            rc['bh_q'] = float(qv[rc['gene']])
            rc['sign_vs_reported'] = None
        out.extend(recs)
    # Additional file 1, Table S1c (DAPT fold change as reported / C0 / C1, two decimals as printed there; transcribed)
    s1c = {'Kcnj11': (0.40, 0.40, 0.76), 'Cacna1d': (0.76, 0.76, 1.46), 'Cacna1c': (0.71, 0.72, 1.37), 'Ptbp1': (0.62, 0.62, 1.18), 'Tcf7l2': (0.78, 0.78, 1.49)}
    byv = {(r['gene'], r['variant']): r['fc_arith'] for r in out}
    GATE['gapdh_variants_match_Table_S1c_at_2_decimals'] = all(round(byv[(g, v)], 2) == s1c[g][k] for g in G5 for k, v in enumerate(('reported', 'C0', 'C1')))
    GATE['gapdh_variants_max_abs_diff_fc_vs_gapdh_exact_json'] = float(max(c[0] for c in chk))
    GATE['gapdh_variants_max_rel_diff_bonf_p_vs_gapdh_exact_json'] = float(max(c[1] for c in chk))
    GATE['gapdh_refs'] = dict(C0=REF['C0'], C1=REF['C1'], hk_used=HK.tolist())
    return out


# ---------------------------------------------------------------- other readouts
def ratio_rec(readout, feature, contrast, a_trt, a_ctl, test, p_raw=None, p_adj=None, note='', variant='as reported'):
    a_trt, a_ctl = np.asarray(a_trt, float), np.asarray(a_ctl, float)
    d, lo, hi, df, pw = welch(np.log2(a_trt), np.log2(a_ctl))
    return dict(readout=readout, feature=feature, category='', contrast=contrast, variant=variant, flag_26Aug_Gapdh=False, test_in_manuscript=test,
                n_ctrl=len(a_ctl), n_trt=len(a_trt), ctrl_mean_reported='', ctrl_sem_reported='', trt_mean_reported='', trt_sem_reported='',
                trt_mean_recomputed='', trt_sem_recomputed='', fc_ratio_arith=float(a_trt.mean() / a_ctl.mean()), log2fc_geo=d, ci95_lo=lo, ci95_hi=hi,
                ci_method='Welch t, log2 scale', df=df, fc_geo=2 ** d, fc_geo_lo=2 ** lo, fc_geo_hi=2 ** hi, ci95_lo_pooled3='', ci95_hi_pooled3='', df_pooled3='',
                ci_excludes_0=(lo > 0 or hi < 0), ci_excludes_0_pooled3='', p_welch_log2=pw, p_raw_manuscript='' if p_raw is None else p_raw,
                p_adj_bonf_manuscript='' if p_adj is None else p_adj, p_raw_recomputed='', bh_q_17genes='', interval_label=LABEL, note=note)


def other_readouts():
    wb = openpyxl.load_workbook(AF4, read_only=True, data_only=True)
    recs = []
    extra = {}
    # ---- DNA content: the four values per group listed in Table S7 (wells at 2.8, 5.5, 22 and 44 mM); basal well = sensitivity
    s7 = af4_rows(wb, 'SuppTable7_Cpeptide_all_gluc')
    h = next(i for i, r in enumerate(s7) if r and r[0] == 'group')
    col = {c: k for k, c in enumerate(s7[h])}
    dna = {gname: {} for gname in ('CTRL', 'DAPT', 'DKK-1')}
    for r in s7[h + 1:]:
        if r and r[0] in dna:
            dna[r[0]][float(r[col['condition_mM']])] = float(r[col['well_DNA_ug']])
    four = {g: [v[c] for c in (2.8, 5.5, 22.0, 44.0)] for g, v in dna.items()}
    five = {g: [v[c] for c in (0.0, 2.8, 5.5, 22.0, 44.0)] for g, v in dna.items()}
    extra['dna_four'] = four
    extra['dna_five'] = five
    rows = af4_rows(wb, 'Fig_1c_stats')
    rep = stat_means(rows)
    pw = stat_pairwise(rows)
    GATE['DNA_mean_max_abs_diff_vs_Fig_1c_stats'] = float(max(abs(np.mean(four[g]) - rep[g]['mean']) for g in four))
    GATE['DNA_sem_max_abs_diff_vs_Fig_1c_stats'] = float(max(abs(sem(four[g]) - rep[g]['sem']) for g in four))
    grp = [np.array(four[g]) for g in ('CTRL', 'DAPT', 'DKK-1')]
    GATE['DNA_bonf_p_max_abs_diff_vs_Fig_1c_stats'] = float(max(abs(min(1.0, 3 * pairwise_p_anova(grp, 0, k)) - pw[nm][1]) for k, nm in ((1, 'DAPT'), (2, 'DKK-1'))))
    for nm, k in (('DAPT', 1), ('DKK-1', 2)):
        rc = ratio_rec('DNA content', 'Total DNA (ug), wells at 2.8, 5.5, 22, 44 mM glucose', nm + ' vs CTRL', grp[k], grp[0], 'One-way ANOVA, Bonferroni t',
                       p_raw=pw[nm][0], p_adj=pw[nm][1],
                       note='four wells that differed in their 20-minute glucose exposure; the fifth (basal) well was left out by the authors for an unrecorded reason')
        rc.update(ctrl_mean_reported=rep['CTRL']['mean'], ctrl_sem_reported=rep['CTRL']['sem'], trt_mean_reported=rep[nm]['mean'], trt_sem_reported=rep[nm]['sem'],
                  trt_mean_recomputed=float(grp[k].mean()), trt_sem_recomputed=sem(grp[k]))
        L = [np.log2(a) for a in grp]
        d3, lo3, hi3, df3, p3 = pooled3(L, 0, k)
        rc.update(ci95_lo_pooled3=lo3, ci95_hi_pooled3=hi3, df_pooled3=df3, ci_excludes_0_pooled3=(lo3 > 0 or hi3 < 0))
        recs.append(rc)
        five_g = [np.array(five[g]) for g in ('CTRL', 'DAPT', 'DKK-1')]
        rc5 = ratio_rec('DNA content', 'Total DNA (ug), all five wells incl. basal (sensitivity)', nm + ' vs CTRL', five_g[k], five_g[0], 'not tested by the authors',
                        variant='sensitivity: five wells', note='post hoc sensitivity; the authors analyzed four wells')
        recs.append(rc5)
    # ---- 2-NBDG: 25,000 cells, second reading, as in Fig. 2b (S1), and the 11 other scenarios of Table S6
    nb = af4_rows(wb, 'Fig_2b_2NBDG_raw_blocks')
    h = next(i for i, r in enumerate(nb) if r and r[0] == 'block')
    blocks = {}
    for r in nb[h + 1:]:
        if r and r[0]:
            blocks.setdefault(r[0], []).append(dict(rep=int(r[1]), blank=float(r[2]), CTRL=float(r[3]), DKK1=float(r[4]), DAPT=float(r[5])))
    sc = af4_rows(wb, 'SuppTable6_2NBDG_scenarios')
    h6 = next(i for i, r in enumerate(sc) if r and r[0] == 'scenario')
    sc_rep = {r[0].split(' ')[0]: r for r in sc[h6 + 1:] if r and r[0]}
    sc_col = {c: k for k, c in enumerate(sc[h6])}

    def vals(block, grp_name, ctrl2=None, drop_high_ctrl=False, blank='none'):
        b = blocks[block]
        v = [x[grp_name] for x in b]
        bl = [x['blank'] for x in b]
        if grp_name == 'CTRL' and ctrl2 is not None:
            v[1] = ctrl2
        if blank == 'row':
            v = [a - c for a, c in zip(v, bl)]
        elif blank == 'mean':
            v = [a - float(np.mean(bl)) for a in v]
        if grp_name == 'CTRL' and drop_high_ctrl:
            v = [a for a in v if a != max(v)]
        return np.array(v, float)
    spec = {   # scenario -> (block, CTRL rep2 override, drop highest CTRL, blank correction); from the scenario names of Table S6
        'S1': ('25K_read2', 331.0, True, 'none'), 'S2': ('25K_read2', None, True, 'none'), 'S3': ('25K_read2', 331.0, False, 'none'),
        'S4': ('25K_read2', None, False, 'none'), 'S5': ('25K_read2', None, True, 'row'), 'S6': ('25K_read2', None, False, 'row'),
        'S7': ('25K_read2', None, False, 'mean'), 'S8': ('25K_read1', None, False, 'none'), 'S9': ('25K_read1', None, True, 'none'),
        'S10': ('50K_read1', None, False, 'none'), 'S11': ('50K_read2', None, False, 'none'), 'S12': ('50K_read2', None, False, 'row')}
    nb_rows = []
    gate_nb = []
    for s, (blk, c2, drop, bl) in spec.items():
        c = vals(blk, 'CTRL', c2, drop, bl)
        # drop applies to the highest CTRL value only for the CTRL group, as in the scenario names
        t = {'DAPT': vals(blk, 'DAPT', None, False, bl), 'DKK-1': vals(blk, 'DKK1', None, False, bl)}
        norm = c.mean()
        r6 = sc_rep[s]
        gate_nb.append(max(abs(c.mean() / norm - float(r6[sc_col['CTRL_mean']])), abs(sem(c / norm) - float(r6[sc_col['CTRL_sem']])),
                           abs((t['DAPT'] / norm).mean() - float(r6[sc_col['DAPT_mean']])), abs(sem(t['DAPT'] / norm) - float(r6[sc_col['DAPT_sem']])),
                           abs((t['DKK-1'] / norm).mean() - float(r6[sc_col['DKK1_mean']])), abs(sem(t['DKK-1'] / norm) - float(r6[sc_col['DKK1_sem']]))))
        if min(c.min(), t['DAPT'].min(), t['DKK-1'].min()) <= 0:
            continue
        for nm in ('DAPT', 'DKK-1'):
            rc = ratio_rec('2-NBDG', '2-NBDG fluorescence, scenario ' + s, nm + ' vs CTRL', t[nm], c, 'Kruskal-Wallis, Bonferroni Wilcoxon',
                           variant='published (Fig. 2b)' if s == 'S1' else 'Table S6 scenario ' + s,
                           note='CTRL values used: ' + ', '.join('%g' % v for v in c) + '; published analysis = S1' if s == 'S1' else 'CTRL values used: ' + ', '.join('%g' % v for v in c))
            rc.update(trt_mean_recomputed=float((t[nm] / norm).mean()), trt_sem_recomputed=sem(t[nm] / norm), ctrl_mean_reported=float(r6[sc_col['CTRL_mean']]),
                      ctrl_sem_reported=float(r6[sc_col['CTRL_sem']]),
                      trt_mean_reported=float(r6[sc_col['DAPT_mean' if nm == 'DAPT' else 'DKK1_mean']]), trt_sem_reported=float(r6[sc_col['DAPT_sem' if nm == 'DAPT' else 'DKK1_sem']]))
            nb_rows.append(rc)
    GATE['NBDG_scenarios_max_abs_diff_mean_sem_vs_Table_S6'] = float(max(gate_nb))
    GATE['NBDG_scenarios_checked'] = len(gate_nb)
    # published S1 p values (Fig_2b_stats)
    rows = af4_rows(wb, 'Fig_2b_stats')
    pw2 = stat_pairwise(rows)
    for rc in nb_rows:
        if rc['variant'].startswith('published'):
            nm = rc['contrast'].split(' ')[0]
            rc['p_raw_manuscript'], rc['p_adj_bonf_manuscript'] = pw2[nm]
    recs.extend(nb_rows)
    # ---- C-peptide
    cp = af4_rows(wb, 'Fig_3b_Cpeptide_all_conditions')
    h = next(i for i, r in enumerate(cp) if r and r[0] == 'group')
    cc = {c: k for k, c in enumerate(cp[h])}
    W = {}
    for r in cp[h + 1:]:
        if r and r[0] in ('CTRL', 'DAPT', 'DKK-1'):
            W.setdefault((r[0], float(r[cc['condition_mM']])), []).append(dict(dod=float(r[cc['dOD']]), per=float(r[cc['ng_per_ug_groupmeanDNA']]), od590=float(r[cc['OD590']]),
                                                                              od450=float(r[cc['OD450']]), gdna=float(r[cc['group_mean_DNA_ug']])))
    # spread of the four wells per group at 5.5 mM (all wells above the lowest standard): coefficient of variation of dOD and the OD450 readings
    extra['cpeptide_5.5mM_cv_dod'] = {g: float(np.std([w['dod'] for w in W[(g, 5.5)]], ddof=1) / np.mean([w['dod'] for w in W[(g, 5.5)]])) for g in ('CTRL', 'DAPT', 'DKK-1')}
    extra['cpeptide_5.5mM_od450'] = {g: [w['od450'] for w in W[(g, 5.5)]] for g in ('CTRL', 'DAPT', 'DKK-1')}
    stats_rows = af4_rows(wb, 'Fig_3b_stats')
    pw3 = {'per': stat_pairwise(stats_rows), 'sig': stat_pairwise(stats_rows, 'Signal (dOD)')}
    blocks = {'per': stat_means(stats_rows), 'sig': stat_means(stats_rows, 'Signal (dOD)')}
    supported = []
    below = {}
    for cond in (0.0, 2.8, 5.5, 22.0, 44.0):
        below[cond] = {g: int(sum(w['dod'] < 0.044 - 1e-9 for w in W[(g, cond)])) for g in ('CTRL', 'DAPT', 'DKK-1')}
        if all(v == 0 for v in below[cond].values()):
            supported.append(cond)
    extra['cpeptide_below_floor'] = {str(k): v for k, v in below.items()}
    extra['cpeptide_supported_conditions'] = supported
    # denominators: five DNA values per group (Table S7); group-mean DNA is the workbook's denominator
    gd = {g: np.array(five[g]) for g in five}
    for cond in supported:
        for key, lab in (('per', 'C-peptide per ug group-mean DNA (workbook units), %g mM glucose' % cond), ('dod', 'ELISA signal dOD, %g mM glucose' % cond)):
            ctl = np.array([w[key] for w in W[('CTRL', cond)]])
            for nm in ('DAPT', 'DKK-1'):
                trt = np.array([w[key] for w in W[(nm, cond)]])
                pp = pw3['per' if key == 'per' else 'sig'].get(nm)
                rc = ratio_rec('C-peptide', lab, nm + ' vs CTRL', trt, ctl, 'Kruskal-Wallis, Bonferroni Wilcoxon', p_raw=None if pp is None else pp[0], p_adj=None if pp is None else pp[1],
                               note='wells only; the DNA denominator is a constant' if key == 'per' else '')
                rc['trt_mean_recomputed'] = float(trt.mean())
                rc['trt_sem_recomputed'] = sem(trt)
                rc['ctrl_mean_reported'] = float(ctl.mean())
                rc['ctrl_sem_reported'] = sem(ctl)
                recs.append(rc)
                if key == 'per':
                    # interval that also carries the uncertainty of the DNA denominator (mean of five DNA values per group), delta method on the log scale
                    lt, lc = np.log(trt), np.log(ctl)
                    comp = [lt.var(ddof=1) / 4, lc.var(ddof=1) / 4, np.log(gd[nm]).var(ddof=1) / 5, np.log(gd['CTRL']).var(ddof=1) / 5]
                    dfc = [3, 3, 4, 4]
                    se = math.sqrt(sum(comp))
                    dfw = se ** 4 / sum(c ** 2 / d for c, d in zip(comp, dfc))
                    dd = (lt.mean() - lc.mean()) / math.log(2)
                    t = stats.t.ppf(0.975, dfw)
                    recs.append(dict(rc, variant='with DNA-denominator uncertainty', log2fc_geo=dd, ci95_lo=dd - t * se / math.log(2), ci95_hi=dd + t * se / math.log(2), df=dfw,
                                     fc_geo=2 ** dd, fc_geo_lo=2 ** (dd - t * se / math.log(2)), fc_geo_hi=2 ** (dd + t * se / math.log(2)),
                                     ci_excludes_0=(dd - t * se / math.log(2) > 0 or dd + t * se / math.log(2) < 0), p_welch_log2='', p_raw_manuscript='', p_adj_bonf_manuscript='',
                                     note='post hoc: Welch-Satterthwaite combination of the log variances of the four wells per group and of the five DNA values per group (the group means used as denominators); '
                                          'the DNA values come from 5 different wells, one per glucose condition'))
                    # sensitivity: CTRL well with the optical artefact left out (OD590 0.105 against 0.036-0.040), post hoc
                    if cond == 5.5:
                        keep = np.array([w['od590'] < 0.08 for w in W[('CTRL', cond)]])
                        if (~keep).sum() == 1:
                            recs.append(dict(ratio_rec('C-peptide', lab, nm + ' vs CTRL', trt, ctl[keep], 'not tested by the authors', variant='sensitivity: CTRL artefact well left out',
                                                       note='post hoc; the authors kept the well (Table S7 note); CTRL n = 3'), ctrl_mean_reported='', ctrl_sem_reported=''))
    # gate: 5.5 mM means and SEMs of the four wells against Fig_3b_stats (per ug group-mean DNA and dOD)
    dm = []
    for sect, key in (('per', 'per'), ('sig', 'dod')):
        for g in ('CTRL', 'DAPT', 'DKK-1'):
            v = np.array([w[key] for w in W[(g, 5.5)]])
            dm.append(abs(v.mean() - blocks[sect][g]['mean']))
            dm.append(abs(sem(v) - blocks[sect][g]['sem']))
    GATE['Cpeptide_5.5mM_mean_sem_max_abs_diff_vs_Fig_3b_stats'] = float(max(dm))
    wb.close()
    return recs, extra, W, below, pw3


# ---------------------------------------------------------------- main
def main():
    af3 = load_af3()
    genes, res = rt_panel(af3)
    reproduce_reported(af3, res)
    recs, q = gene_rows(af3, genes, res)
    var = variants_g5(af3, res, q)
    oth, extra, W, below, pw3 = other_readouts()
    # CSV
    allrecs = recs + oth
    cols = ['readout', 'feature', 'category', 'contrast', 'variant', 'flag_26Aug_Gapdh', 'test_in_manuscript', 'n_ctrl', 'n_trt',
            'ctrl_mean_reported', 'ctrl_sem_reported', 'trt_mean_reported', 'trt_sem_reported', 'trt_mean_recomputed', 'trt_sem_recomputed',
            'fc_ratio_arith', 'log2fc_geo', 'ci95_lo', 'ci95_hi', 'ci_method', 'df', 'fc_geo', 'fc_geo_lo', 'fc_geo_hi',
            'ci95_lo_pooled3', 'ci95_hi_pooled3', 'df_pooled3', 'ci_excludes_0', 'ci_excludes_0_pooled3', 'p_welch_log2',
            'p_raw_manuscript', 'p_adj_bonf_manuscript', 'p_raw_recomputed', 'p_adj_bonf_recomputed', 'bh_q_17genes', 'interval_label', 'note']
    # 26 August variants as additional rows
    vrecs = []
    for v in var:
        r0 = next(r for r in af3 if r['gene'] == v['gene'])
        vrecs.append(dict(readout='RT-qPCR', feature=v['gene'], category=r0['category'], contrast='DAPT vs CTRL', variant='Gapdh ' + v['variant'] if v['variant'] != 'reported' else 'as reported (26 Aug Gapdh)',
                          flag_26Aug_Gapdh=True, test_in_manuscript='One-way ANOVA', n_ctrl=4, n_trt=4, trt_mean_recomputed=v['mean'], trt_sem_recomputed=v['sem'],
                          fc_ratio_arith=v['fc_arith'], log2fc_geo=v['log2fc'], ci95_lo=v['lo'], ci95_hi=v['hi'], ci_method='Welch t, log2 scale', df=v['df'],
                          fc_geo=2 ** v['log2fc'], fc_geo_lo=2 ** v['lo'], fc_geo_hi=2 ** v['hi'], ci95_lo_pooled3=v['lo3'], ci95_hi_pooled3=v['hi3'], df_pooled3=9,
                          ci_excludes_0=(v['lo'] > 0 or v['hi'] < 0), ci_excludes_0_pooled3=(v['lo3'] > 0 or v['hi3'] < 0),
                          p_raw_recomputed=v['p_raw'], p_adj_bonf_recomputed=v['p_bonf'], bh_q_17genes=v['bh_q'], interval_label=LABEL,
                          note={'reported': 'Gapdh reference of the 26 Aug DAPT group as analyzed (C10 and D10 entered twice)',
                                'C0': 'Gapdh = mean of the two valid wells C10 and D10 (22.68) for all four DAPT replicates',
                                'C1': 'Gapdh = mean of all four wells A10-D10 (23.62), failed wells included (worst case)'}[v['variant']]))
    with open(HERE / 'effect_sizes.csv', 'w', newline='', encoding='utf8') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in allrecs + [x for x in vrecs if not x['variant'].startswith('as reported')]:
            w.writerow({k: (('%.6g' % v) if isinstance(v, float) else v) for k, v in r.items()})
    json.dump(dict(gate=GATE, extra=extra, below_floor=below), open(HERE / 'gate_reproduction.json', 'w', encoding='utf8'), indent=1, default=str)
    write_md(af3, genes, res, recs, var, oth, extra, below)
    print(json.dumps(GATE, indent=1, default=str))


def write_md(af3, genes, res, recs, var, oth, extra, below):
    """effect_sizes.md"""
    L = []
    A = L.append

    def gi(g):
        return '*%s*' % g
    A('# Effect sizes with 95% intervals, across-gene BH q and Gapdh sensitivity (WP5a)')
    A('')
    A('Produced by `effect_sizes.py` from Additional files 3 and 4 and `gapdh_exact.json`. Every interval below is **%s**.' % LABEL)
    A('')
    A('## How to read the tables')
    A('')
    A('* **Values used.** For each gene, the 12 relative expression values (2^-dCt, one per well-level replicate; Additional file 4, sheet `Fig_4c_Raw_expression_data`, full precision) divided by the arithmetic mean of the four CTRL values of that gene. These are the replicate fold changes R1-R4 of Additional file 3 (reproduced to rounding, see README gate). Rows 1-4 are CTRL, 5-8 DAPT, 9-12 DKK-1. DAPT values of the five 26 August genes are as reported unless a variant is named.')
    A('* **log2 FC [95% CI].** Each replicate fold change is log2-transformed. Effect = mean log2 (treated) minus mean log2 (CTRL), that is log2 of the ratio of geometric means. Interval = Welch t-interval (n = 4 vs 4; Welch-Satterthwaite df, range in the summary below), so no equal-variance assumption is made. Not adjusted for the three pairwise contrasts per gene or for the 17 genes. '
      'The same values rounded to the 3 decimals printed in Additional file 3 move the interval limits by at most ' + '%.4f' % GATE['ci_limits_full_vs_printed_replicates_max_abs_diff_log2'] + ' log2 units.')
    A('* **Mean +/- SEM.** As printed in Additional file 3 (fold change of the arithmetic mean relative to the CTRL arithmetic mean; SEM of the four replicate fold changes). The manuscript fold change is the ratio of arithmetic means, so it differs slightly from 2^(log2 FC).')
    A('* **p (Bonf.).** The manuscript\'s Bonferroni-corrected pairwise p (Additional file 3): pooled-variance t-tests after one-way ANOVA, or Wilcoxon rank-sum tests after Kruskal-Wallis for Wnt5a and Lef1 (marked KW; the smallest attainable adjusted p is 0.086).')
    G1 = GATE['from_printed_replicates']
    A('* **BH q.** Benjamini-Hochberg adjustment over the 17 raw pairwise p of the same contrast (CTRL vs DAPT, or CTRL vs DKK-1). The raw p are those of the manuscript\'s tests (pooled-variance t-test, or Wilcoxon rank-sum for the two KW genes) recomputed from the 12 full-precision values. The p printed in Additional file 3 were computed from its 3-decimal replicate values (all %d of %d printed pairwise p are reproduced from those values); the full-precision p differ from them by a median of %.1f%% and at most %.0f%% (relative), and q changes by at most %.4f with no change in whether q < 0.05. q is reported beside the manuscript\'s tests and does not replace them; it answers a different question (17 genes) from the Bonferroni p (3 contrasts per gene). q < 0.05 is marked in bold.' % (
        G1['pairwise_p_matching_printed'], G1['pairwise_p_checked'], 100 * GATE['p_raw_full_vs_printed_replicates_median_rel_diff'], 100 * GATE['p_raw_full_vs_printed_replicates_max_rel_diff'],
        max(GATE['bh_q_max_abs_diff_full_vs_printed_replicates_DAPT'], GATE['bh_q_max_abs_diff_full_vs_printed_replicates_DKK-1'])))
    A('* **Excludes 0 (pooled).** Sensitivity: whether the interval excludes 0 when the within-group variance is pooled over the three groups (df = 9), the model behind the manuscript\'s pairwise tests.')
    A('* **†** The DAPT value depends on the Gapdh wells of the 26 August run (Tcf7l2, Kcnj11, Cacna1d, Cacna1c, Ptbp1); see the section on Gapdh variants.')
    A('')
    # summary counts
    A('## Summary of the 17-gene panel')
    A('')
    allp = [x for x in recs]
    dfs = [x['df'] for x in allp]
    widths = [x['ci95_hi'] - x['ci95_lo'] for x in allp]
    A('Welch-Satterthwaite df range %.1f to %.1f; median width of the 95%% interval %.2f log2 units (a %.1f-fold range). Counts below are out of 17 genes per contrast.' % (
        min(dfs), max(dfs), float(np.median(widths)), 2 ** float(np.median(widths))))
    A('')
    A('| Contrast | Manuscript Bonferroni p < 0.05 | BH q < 0.05 | Welch interval excludes 0 | Pooled interval excludes 0 | Genes where Bonferroni and BH disagree | Genes where Welch interval and Bonferroni disagree |')
    A('|---|---|---|---|---|---|---|')
    for key in CONTRASTS:
        rr = [x for x in recs if x['contrast'] == key + ' vs CTRL']
        bon = {x['feature']: float(x['p_adj_bonf_manuscript']) < 0.05 for x in rr}
        qq = {x['feature']: x['bh_q_17genes'] < 0.05 for x in rr}
        cw = {x['feature']: bool(x['ci_excludes_0']) for x in rr}
        cp = {x['feature']: bool(x['ci_excludes_0_pooled3']) for x in rr}
        d1 = ['%s (%s)' % (g, 'q only' if qq[g] else 'Bonferroni only') for g in bon if bon[g] != qq[g]]
        d2 = ['%s (%s)' % (g, 'interval only' if cw[g] else 'Bonferroni only') for g in bon if bon[g] != cw[g]]
        A('| %s vs CTRL | %d | %d | %d | %d | %s | %s |' % (key, sum(bon.values()), sum(qq.values()), sum(cw.values()), sum(cp.values()),
                                                         ', '.join('*%s' % d.replace(' (', '* (', 1) for d in d1) or 'none', ', '.join('*%s' % d.replace(' (', '* (', 1) for d in d2) or 'none'))
    A('')
    A('Bonferroni p < 0.05 is the manuscript\'s own criterion. The intervals are unadjusted, so a few genes whose Bonferroni p is just above 0.05 have intervals that exclude 0.')
    A('')
    for key in CONTRASTS:
        A('## RT-qPCR: %s vs CTRL (17 genes)' % key)
        A('')
        A('Interval label: %s.' % LABEL)
        A('')
        A('| Gene | Test | CTRL mean +/- SEM | ' + key + ' mean +/- SEM | log2 FC [95% CI] | 2^log2 FC [95% CI] | p (Bonf.), manuscript | BH q | Excludes 0 (Welch / pooled) |')
        A('|---|---|---|---|---|---|---|---|---|')
        for r in [x for x in recs if x['contrast'] == key + ' vs CTRL']:
            g = r['feature']
            dag = '†' if r['flag_26Aug_Gapdh'] else ''
            qv = r['bh_q_17genes']
            qs = ('**%s**' % fmtp(qv)) if qv < 0.05 else fmtp(qv)
            A('| %s%s | %s | %s +/- %s | %s +/- %s | %s [%s, %s] | %s [%s, %s] | %s | %s | %s / %s |' % (
                gi(g), dag, 'KW' if r['test_in_manuscript'].startswith('Kruskal') else 'ANOVA', '%g' % r['ctrl_mean_reported'],
                '%g' % r['ctrl_sem_reported'], '%g' % r['trt_mean_reported'], '%g' % r['trt_sem_reported'],
                fmt(r['log2fc_geo'], 2), fmt(r['ci95_lo'], 2), fmt(r['ci95_hi'], 2), fmt(r['fc_geo'], 2), fmt(r['fc_geo_lo'], 2), fmt(r['fc_geo_hi'], 2),
                r['p_adj_bonf_manuscript'], qs, 'yes' if r['ci_excludes_0'] else 'no', 'yes' if r['ci_excludes_0_pooled3'] else 'no'))
        A('')
    # Gapdh variants
    A('## The five genes whose DAPT values depend on the 26 August Gapdh wells')
    A('')
    A('DAPT vs CTRL; CTRL and DKK-1 values are identical in all three variants. Reported = Gapdh values as analyzed in the authors\' calculation sheet (C10 and D10 entered twice, A10 and B10 set aside). C0 = reference is the mean of the two valid wells (Cq 22.68) for all four DAPT replicates. C1 = reference is the mean of all four wells (Cq 23.62), failed wells included (worst case). Same interval method as above, label: %s. BH q is recomputed over the 17 genes of the DAPT contrast with the five p replaced by the variant values.' % LABEL)
    A('')
    A('| Gene | Variant | Mean +/- SEM | log2 FC [95% CI] | 2^log2 FC [95% CI] | p (Bonf.) | BH q | Sign of log2 FC |')
    A('|---|---|---|---|---|---|---|---|')
    for v in var:
        qs = ('**%s**' % fmtp(v['bh_q'])) if v['bh_q'] < 0.05 else fmtp(v['bh_q'])
        sgn = 'down' if v['log2fc'] < 0 else 'up'
        A('| %s | %s | %.3f +/- %.3f | %s [%s, %s] | %s [%s, %s] | %s | %s | %s |' % (
            gi(v['gene']), {'reported': 'reported', 'C0': 'C0 (two valid wells)', 'C1': 'C1 (all four wells)'}[v['variant']], v['mean'], v['sem'],
            fmt(v['log2fc'], 2), fmt(v['lo'], 2), fmt(v['hi'], 2), fmt(2 ** v['log2fc'], 2), fmt(2 ** v['lo'], 2), fmt(2 ** v['hi'], 2), fmtp(v['p_bonf']), qs, sgn))
    A('')
    vr = {v['gene']: v for v in var if v['variant'] == 'reported'}
    v0 = {v['gene']: v for v in var if v['variant'] == 'C0'}
    v1 = {v['gene']: v for v in var if v['variant'] == 'C1'}
    flip = [g for g in G5 if (v1[g]['log2fc'] > 0) != (vr[g]['log2fc'] > 0)]
    inc0_c1 = [g for g in G5 if v1[g]['lo'] <= 0 <= v1[g]['hi']]
    inc0_c0 = [g for g in G5 if v0[g]['lo'] <= 0 <= v0[g]['hi']]
    inc0_rep = [g for g in G5 if vr[g]['lo'] <= 0 <= vr[g]['hi']]
    A('The "reported" rows are the main-table values; all p and q in this table are recomputed from the full-precision values (pooled-variance t-tests, df = 9), so their third digit can differ from the printed p of Additional file 3. Fold changes and Bonferroni p of all three variants match `gapdh_exact.json` (max absolute difference of fold changes %.1e, max relative difference of p %.1e). '
      'With C1 the sign of the log2 FC reverses for %s. The 95%% interval includes 0 under C1 for %s, under C0 for %s and under the reported values for %s.' % (
        GATE['gapdh_variants_max_abs_diff_fc_vs_gapdh_exact_json'], GATE['gapdh_variants_max_rel_diff_bonf_p_vs_gapdh_exact_json'],
        ', '.join('*%s*' % g for g in flip) or 'none', ', '.join('*%s*' % g for g in inc0_c1) or 'none', ', '.join('*%s*' % g for g in inc0_c0) or 'none', ', '.join('*%s*' % g for g in inc0_rep) or 'none'))
    A('')
    # other readouts
    A('## DNA content, 2-NBDG and C-peptide per ug DNA')
    A('')
    A('Interval label: %s. Intervals are Welch t-intervals on the log2 scale of the replicate values, as above.' % LABEL)
    A('')
    A('| Readout | Contrast | Variant | n (CTRL / treated) | Mean ratio (arith.) | log2 ratio [95% CI] | 2^log2 ratio [95% CI] | p (Bonf.), manuscript |')
    A('|---|---|---|---|---|---|---|---|')
    for r in oth:
        if r['readout'] == '2-NBDG' and not r['variant'].startswith('published'):
            continue
        pj = r['p_adj_bonf_manuscript']
        pj = '' if pj == '' else (fmtp(float(pj)) if isinstance(pj, (int, float)) else str(pj))
        A('| %s | %s | %s | %d / %d | %s | %s [%s, %s] | %s [%s, %s] | %s |' % (
            r['feature'], r['contrast'], r['variant'], r['n_ctrl'], r['n_trt'], fmt(r['fc_ratio_arith'], 3), fmt(r['log2fc_geo'], 2), fmt(r['ci95_lo'], 2), fmt(r['ci95_hi'], 2),
            fmt(r['fc_geo'], 2), fmt(r['fc_geo_lo'], 2), fmt(r['fc_geo_hi'], 2), pj))
    A('')

    def pick(readout, feat, contrast, variant):
        return next(x for x in oth if x['readout'] == readout and feat in x['feature'] and x['contrast'] == contrast and x['variant'] == variant)
    dna_d, dna_k = pick('DNA content', 'wells at 2.8', 'DAPT vs CTRL', 'as reported'), pick('DNA content', 'wells at 2.8', 'DKK-1 vs CTRL', 'as reported')
    dna5_d, dna5_k = pick('DNA content', 'five wells', 'DAPT vs CTRL', 'sensitivity: five wells'), pick('DNA content', 'five wells', 'DKK-1 vs CTRL', 'sensitivity: five wells')
    nbs = {}
    for x in oth:
        if x['readout'] == '2-NBDG':
            nbs.setdefault(x['contrast'].split(' ')[0], []).append(x)
    n_ex = sum(1 for k in nbs for x in nbs[k] if x['ci_excludes_0'])
    n_all = sum(len(v) for v in nbs.values())
    cp_d, cp_k = pick('C-peptide', 'per ug', 'DAPT vs CTRL', 'as reported'), pick('C-peptide', 'per ug', 'DKK-1 vs CTRL', 'as reported')
    cs_d, cs_k = pick('C-peptide', 'dOD', 'DAPT vs CTRL', 'as reported'), pick('C-peptide', 'dOD', 'DKK-1 vs CTRL', 'as reported')
    A('Notes on this table.')
    A('')
    A('* **DNA content.** The four values per group are wells that differed in their 20-minute glucose exposure (Methods), not four replicates of one condition, and the CTRL wells range from %s to %s ug. The lower interval limit is %s (DAPT) and %s (DKK-1) log2 units above 0, that is, barely above no change. Adding the fifth (basal) well, which the authors left out for an unrecorded reason, gives %s and %s log2 units with lower limits %s and %s (post hoc sensitivity, not a replacement).' % (
        '%.3f' % min(extra['dna_four']['CTRL']), '%.3f' % max(extra['dna_four']['CTRL']),
        fmt(dna_d['ci95_lo'], 2), fmt(dna_k['ci95_lo'], 2), fmt(dna5_d['log2fc_geo'], 2), fmt(dna5_k['log2fc_geo'], 2), fmt(dna5_d['ci95_lo'], 2), fmt(dna5_k['ci95_lo'], 2)))
    A('* **2-NBDG: not supported.** Across the 12 scenarios of Table S6, %d of %d intervals exclude 0, the log2 ratio of DAPT runs from %s to %s and that of DKK-1 from %s to %s, so even the direction depends on the analysis choice. The published analysis (S1) uses the CTRL values 195, 331 and 331 (the value 331 appears twice because the analysis block of the workbook holds 331 for CTRL replicate 2, where the reading block holds 184; the reason is not recorded). Only S1 is shown in the table above; all twelve are in the scenario table below and in the CSV.' % (
        n_ex, n_all, fmt(min(x['log2fc_geo'] for x in nbs['DAPT']), 2), fmt(max(x['log2fc_geo'] for x in nbs['DAPT']), 2), fmt(min(x['log2fc_geo'] for x in nbs['DKK-1']), 2), fmt(max(x['log2fc_geo'] for x in nbs['DKK-1']), 2)))
    A('* **C-peptide at 5.5 mM glucose.** The raw signal (dOD) is higher than CTRL in both treated groups (mean ratio %s and %s) whereas the value per ug of group-mean DNA is lower (%s and %s): the sign of the difference is set by the DNA denominator, which is 2.3 and 2.6 times larger in the treated groups (Table S7). Three of the four CTRL wells read 0.047 to 0.048, only 0.003 to 0.004 above the lowest standard (0.044), and the fourth (0.083) carries the optical artefact noted in Table S7. The "wells only" interval treats the DNA denominator as an exact constant; the row "with DNA-denominator uncertainty" adds the variance of the five DNA values per group (one well per glucose condition) and is much wider. The narrow intervals of the row "CTRL artefact well left out" come from near-identical readings (the readings are 0.001 apart in dOD) and describe the spread of the readings, not the uncertainty of the biology.' % (
        fmt(cs_d['fc_ratio_arith'], 2), fmt(cs_k['fc_ratio_arith'], 2), fmt(cp_d['fc_ratio_arith'], 2), fmt(cp_k['fc_ratio_arith'], 2)))
    A('')
    A('### 2-NBDG: range over the twelve analysis scenarios of Table S6')
    A('')
    A('| Scenario | CTRL values used | DAPT mean ratio | DAPT log2 ratio [95% CI] | DKK-1 mean ratio | DKK-1 log2 ratio [95% CI] |')
    A('|---|---|---|---|---|---|')
    sc = {}
    for r in oth:
        if r['readout'] == '2-NBDG':
            s = r['feature'].split('scenario ')[1]
            sc.setdefault(s, {})[r['contrast'].split(' ')[0]] = r
    for s, d in sc.items():
        A('| %s | %s | %s | %s [%s, %s] | %s | %s [%s, %s] |' % (s, d['DAPT']['note'].split('CTRL values used: ')[1].split(';')[0], fmt(d['DAPT']['fc_ratio_arith'], 2),
                                                               fmt(d['DAPT']['log2fc_geo'], 2), fmt(d['DAPT']['ci95_lo'], 2), fmt(d['DAPT']['ci95_hi'], 2),
                                                               fmt(d['DKK-1']['fc_ratio_arith'], 2), fmt(d['DKK-1']['log2fc_geo'], 2), fmt(d['DKK-1']['ci95_lo'], 2), fmt(d['DKK-1']['ci95_hi'], 2)))
    A('')
    A('### C-peptide: which glucose conditions the data support')
    A('')
    A('A comparison is estimated only where no well of any of the three groups lies below the lowest standard (dOD < 0.044). Wells below the lowest standard per group (CTRL / DAPT / DKK-1):')
    A('')
    A('| Glucose (mM) | Below lowest standard | Estimated here |')
    A('|---|---|---|')
    for cond, d in below.items():
        A('| %s | %d / %d / %d | %s |' % (cond, d['CTRL'], d['DAPT'], d['DKK-1'], 'yes' if all(v == 0 for v in d.values()) else 'no: at least one group has wells below the lowest standard'))
    A('')
    open(HERE / 'effect_sizes.md', 'w', encoding='utf8').write('\n'.join(L))


if __name__ == '__main__':
    main()
