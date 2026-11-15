"""WP1 / protocol v2 section 2 (items 2.2-2.8) on the INS-1 17-gene matrix. POST HOC and exploratory: one cell preparation, n = 12 (3 groups x 4 wells); nothing here is confirmatory.
Needs gene_run_map.csv and ins1_sheet_cq.csv from build_gene_run_map.py (the verified map replaces the provisional map of ../recon_run_structure.py) and ../coexp_checks.py.
Writes results_ins1.json and tbl_*.csv next to this script. Seed 2026 throughout; 20000 label permutations (2.2, 2.3), 2000 adjusted-count permutations (2.5), 10000 for the other nulls (2.6, 2.8).
Items marked ADD-ON are not in the protocol (a diagnostic added after the verified map showed the sheet layout); they are reported as such.
Usage: python -I wp1_ins1_analysis.py        (about 3-4 minutes; env INS1_MATRIX and BENCH_DIR override the matrix and benchmark_route3 paths)"""
import sys, os, json, time
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))          # -I removes the script folder from sys.path
import numpy as np, pandas as pd
from scipy import stats
import coexp_checks as cc

MATRIX = os.environ.get('INS1_MATRIX', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'R_17genes', 'Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx'))
BENCH = os.environ.get('BENCH_DIR', os.path.join(HERE, '..', '..', 'benchmark_route3'))
G17 = ['Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']; DEN = ['Ins1', 'Kcnj11', 'Hes1']; SEED = 2026
PROVISIONAL = {'Glp1r': 'Jun03', 'Glut2': 'Jun03', 'Hes1': 'Jun03', 'Ins1': 'Jun03', 'Lef1': 'Jun24', 'Tcf7': 'Jun24', 'Wnt5a': 'Jun24', 'Wnt5b': 'Jun24', 'Wnt9a': 'Jun24', 'Cacna1c': 'Aug26', 'Cacna1d': 'Aug26',
               'Kcnj11': 'Aug26', 'Ptbp1': 'Aug26', 'Tcf7l2': 'Aug26', 'Itpr1': 'Nov09', 'Ryr2': 'Nov09', 'Ryr3': 'Nov09'}        # recon_run_structure.py (planning guess from Table S1c), SUPERSEDED by the verified map
OUT = {}


def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.floating, float)): return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.bool_,)): return bool(o)
    if isinstance(o, np.ndarray): return clean(o.tolist())
    return o


def save(df, name, nd=4):
    df.to_csv(os.path.join(HERE, name), index=False, float_format='%.' + str(nd + 2) + 'g'); return df


def corr_stack(D):
    """Pearson correlation matrices of the columns of D (B, N, n) over the sample axis"""
    Z = D - D.mean(1, keepdims=True); Z = Z / np.sqrt((Z ** 2).sum(1, keepdims=True)); return np.swapaxes(Z, 1, 2) @ Z


def upper_mean(R): iu = np.triu_indices(R.shape[-1], 1); return R[..., iu[0], iu[1]].mean(-1)


# ------------------------------------------------------------------------------------------------------------------ load
X = pd.read_excel(MATRIX)[G17].astype(float).values; grp = np.repeat(['CTRL', 'DAPT', 'DKK-1'], 4); des = cc.Design(grp); t_start = time.time()
mp = pd.read_csv(os.path.join(HERE, 'gene_run_map.csv')).set_index('gene'); cq = pd.read_csv(os.path.join(HERE, 'ins1_sheet_cq.csv'))
assert all(mp.loc[g, 'status'] == 'resolved' for g in G17), 'a gene is unresolved: exclude it from the run tests before running this script'
runs = [mp.loc[g, 'run'] for g in G17]; refs = [mp.loc[g, 'reference_set'] for g in G17]
T = np.column_stack([cq[cq.gene == g].sort_values('replicate').target_cq.values for g in G17]); RF = np.column_stack([cq[cq.gene == g].sort_values('replicate').reference_cq.values for g in G17])
assert np.allclose(np.exp2(-(T - RF)), X, rtol=1e-12, atol=0), 'sheet Cq values do not reproduce the matrix'
base = cc.observed_summary(X, grp, DEN, G17); iu17 = np.triu_indices(17, 1); NAMES = base['scheme_names']
OUT['status'] = 'POST HOC exploratory analysis of one cell preparation (3 groups x 4 wells); protocol v2 section 2; not confirmatory'
OUT['map'] = {g: dict(run=mp.loc[g, 'run'], reference_set=mp.loc[g, 'reference_set'], sheet=mp.loc[g, 'sheet'], rows=mp.loc[g, 'sheet_rows']) for g in G17}
OUT['run_sizes_17genes'] = pd.Series(runs).value_counts().to_dict(); OUT['reference_set_sizes_17genes'] = pd.Series(refs).value_counts().to_dict()
OUT['baseline'] = dict(consensus_any=base['consensus_any'], consensus_same_sign=base['consensus_same_sign'], pooled_sig_per_scheme=base['pooled_sig_per_scheme'], adj_q05_per_scheme=base['adj_q_lt_alpha_per_scheme'],
                       cons_adj_any=base['cons_adj_any'], degree=base['degree'], consensus_pairs=base['consensus_pairs'], scheme_names=NAMES)
print('loaded; baseline consensus %d, pooled %s, adjusted %s' % (base['consensus_any'], base['pooled_sig_per_scheme'], base['adj_q_lt_alpha_per_scheme']))

# ------------------------------------------------------------------------------------------------------------------ 2.2 run-structure test
print('2.2 run structure'); rows = []


def flatten(res, label):
    for tag in ('r', 'abs_r'):
        t = res[tag]; rows.append(dict(labelling=label, statistic=tag, scope='all runs', n_genes=res['n_genes'], n_inside_pairs=res['n_same_pairs'], n_cross_pairs=res['n_between_pairs'], mean_inside=t['mean_same'], mean_cross=t['mean_between'], D=t['D'], null_sd=t['null_sd'], p_two_sided=t['p_two_sided']))
        for run, v in t['per_run'].items(): rows.append(dict(labelling=label, statistic=tag, scope='run ' + run, n_genes=v['n_genes'], n_inside_pairs=v['n_pairs_inside'], mean_inside=v['mean_inside'], mean_cross=v['mean_cross'], D=v['D'], p_two_sided=v['p_two_sided']))


rs_primary = cc.run_structure_test(X, grp, runs, 20000, SEED, G17); flatten(rs_primary, 'verified map (run = sheet): primary')
rs_refset = cc.run_structure_test(X, grp, refs, 20000, SEED, G17); flatten(rs_refset, 'verified map, run = reference set (ADD-ON secondary: Jun24 splits in two Gapdh vectors)')
rs_prov = cc.run_structure_test(X, grp, [PROVISIONAL[g] for g in G17], 20000, SEED, G17); flatten(rs_prov, 'provisional map of the planning reconnaissance (superseded)')
t22 = save(pd.DataFrame(rows), 'tbl_2_2_run_structure.csv'); OUT['2_2_run_structure'] = dict(primary=rs_primary, reference_set=rs_refset, provisional_map=rs_prov)
r_log = cc.residual_corr(X, grp)                                                                                           # within-group r on the log scale (the statistic of 2.2)
r_lin = cc.residual_corr(X, grp, log=False)                                                                                # treatment-adjusted partial r on the linear scale (the manuscript's adjusted analysis)
print(t22[t22.statistic == 'r'][['labelling', 'scope', 'n_inside_pairs', 'mean_inside', 'mean_cross', 'D', 'p_two_sided']].round(3).to_string())

# ------------------------------------------------------------------------------------------------------------------ 2.3 consensus pairs and hubs by run
print('2.3 consensus by run')
elig = [i for i, g in enumerate(G17) if g not in DEN]; ge = [G17[i] for i in elig]; ne = len(elig); iue = np.triu_indices(ne, 1)
lab_e = [runs[i] for i in elig]; uniq_e = list(dict.fromkeys(lab_e)); code_e = np.array([uniq_e.index(x) for x in lab_e])
consv = base['_cons_any'][np.ix_(elig, elig)][iue].astype(float); same_e = code_e[iue[0]] == code_e[iue[1]]
a, b, c_, d = int((same_e & (consv == 1)).sum()), int((same_e & (consv == 0)).sum()), int((~same_e & (consv == 1)).sum()), int((~same_e & (consv == 0)).sum())
fo, fp = stats.fisher_exact([[a, b], [c_, d]]); fpg = stats.fisher_exact([[a, b], [c_, d]], alternative='greater')[1]
rng = np.random.default_rng(SEED); cp = code_e[np.argsort(rng.random((20000, ne)), axis=1)]; null_same = (cp[:, iue[0]] == cp[:, iue[1]]).astype(float) @ consv; obs_same = float(same_e.astype(float) @ consv)
p_ge = float((1 + (null_same >= obs_same).sum()) / 20001); p_le = float((1 + (null_same <= obs_same).sum()) / 20001); p_two = float((1 + (np.abs(null_same - null_same.mean()) >= abs(obs_same - null_same.mean()) - 1e-12).sum()) / 20001)
deg = np.array([base['degree'].get(g, 0) for g in ge]); nrun = pd.Series(lab_e).value_counts()
hub_rows = []
for i, g in enumerate(ge):
    k = lab_e[i]; s_possible = int(nrun[k] - 1); partners = [j for j in range(ne) if j != i and (base['_cons_any'][elig[min(i, j)], elig[max(i, j)]])]
    same_edges = sum(lab_e[j] == k for j in partners); hub_rows.append(dict(gene=g, run=k, reference_set=refs[elig[i]], genes_in_run_17=int(sum(r == k for r in runs)), eligible_genes_in_run=int(nrun[k]), degree=int(deg[i]), is_hub_deg_ge_5=bool(deg[i] >= 5),
                                                                     same_run_edges=int(same_edges), between_run_edges=int(deg[i] - same_edges), possible_same_run_partners=s_possible, expected_same_run_edges=float(deg[i] * s_possible / (ne - 1))))
hubs = save(pd.DataFrame(hub_rows).sort_values(['degree', 'gene'], ascending=[False, True]), 'tbl_2_3_hubs_vs_run.csv')
runrows = []
for k in uniq_e:
    m = np.array([l == k for l in lab_e]); n_k = int(m.sum()); n_hub = int((deg[m] >= 5).sum()); n_all_hub = int((deg >= 5).sum())
    mean_deg = float(deg[m].mean()); null_mean_deg = np.array([deg[rng.permutation(ne)[:n_k]].mean() for _ in range(5000)])
    runrows.append(dict(run=k, eligible_genes=n_k, mean_degree=mean_deg, hubs_deg_ge_5=n_hub, hubs_expected=n_all_hub * n_k / ne, hypergeometric_p_ge=float(stats.hypergeom.sf(n_hub - 1, ne, n_all_hub, n_k)), p_mean_degree_ge=float((1 + (null_mean_deg >= mean_deg - 1e-12).sum()) / 5001)))
runtab = save(pd.DataFrame(runrows), 'tbl_2_3_degree_by_run.csv'); H = stats.kruskal(*[deg[code_e == k] for k in range(len(uniq_e))])
nullH = np.array([stats.kruskal(*[deg[code_e[rng.permutation(ne)] == k] for k in range(len(uniq_e))]).statistic for _ in range(3000)]); pH = float((1 + (nullH >= H.statistic - 1e-12).sum()) / 3001)
inv = []                                                                                                                    # consensus pairs that involve at least one gene of a run (descriptive, post hoc)
for k in uniq_e:
    m_inv = np.array([lab_e[i] == k or lab_e[j] == k for i, j in zip(*iue)]); n_inv = int(m_inv.sum()); n_c = int((m_inv & (consv == 1)).sum())
    inv.append(dict(run=k, eligible_pairs_involving_run=n_inv, consensus_pairs_involving_run=n_c, consensus_expected=float(consv.sum() * n_inv / len(consv)), hypergeometric_p_ge=float(stats.hypergeom.sf(n_c - 1, len(consv), n_inv, int(consv.sum())))))
save(pd.DataFrame(inv), 'tbl_2_3_pairs_involving_run.csv')
pairs_e = pd.DataFrame(dict(gene_a=[ge[i] for i in iue[0]], gene_b=[ge[j] for j in iue[1]], run_a=[lab_e[i] for i in iue[0]], run_b=[lab_e[j] for j in iue[1]], same_run=same_e, consensus=consv.astype(bool)))
save(pairs_e, 'tbl_2_3_eligible_pairs.csv')
OUT['2_3_consensus_by_run'] = dict(eligible_pairs=int(len(consv)), same_run_pairs=int(same_e.sum()), between_run_pairs=int((~same_e).sum()), consensus_same_run=a, same_run_not_consensus=b, consensus_between_run=c_, between_run_not_consensus=d,
                                   share_consensus_same_run=a / max(a + b, 1), share_consensus_between_run=c_ / max(c_ + d, 1), fisher_odds_ratio=float(fo), fisher_p_two_sided=float(fp), fisher_p_greater=float(fpg),
                                   label_permutation_null_mean=float(null_same.mean()), label_permutation_null_sd=float(null_same.std()), label_permutation_p_ge=p_ge, label_permutation_p_le=p_le, label_permutation_p_two_sided=p_two,
                                   pairs_involving_run=inv, hubs=[g for g, k in zip(ge, deg) if k >= 5], hub_table=hub_rows, degree_by_run=runrows, kruskal_H_degree_by_run=float(H.statistic), kruskal_p_asymptotic=float(H.pvalue), kruskal_p_label_permutation=pH)
print('  consensus same-run %d/%d, between-run %d/%d; Fisher OR %.2f P %.3f; label permutation P(>=) %.4f (null mean %.1f); degree by run: %s' % (a, a + b, c_, c_ + d, fo, fp, p_ge, null_same.mean(), {r['run']: round(r['mean_degree'], 2) for r in runrows}))

# ------------------------------------------------------------------------------------------------------------------ 2.4 shared-reference contribution
print('2.4 shared reference')
r_cq = cc.residual_corr(T, grp, log=False); r_dcq = cc.residual_corr(X, grp); r_ref = cc.residual_corr(RF, grp, log=False)
assert np.allclose(r_dcq, r_log)
# analytic expectation of the within-group r of two dCq series if the TARGET Cq values were independent and only the (observed) reference residuals were shared
Qm = des.Q; E = Qm @ RF; Tr = Qm @ T; dfq = des.N - des.K; svar = (Tr ** 2).sum(0)                                      # residual sums of squares of the targets
num = E.T @ E; den_ = np.sqrt(np.outer(svar + (E ** 2).sum(0), svar + (E ** 2).sum(0))); r_exp = num / den_
P = pd.DataFrame(dict(gene_a=[G17[i] for i in iu17[0]], gene_b=[G17[j] for j in iu17[1]], run_a=[runs[i] for i in iu17[0]], run_b=[runs[j] for j in iu17[1]], ref_a=[refs[i] for i in iu17[0]], ref_b=[refs[j] for j in iu17[1]]))
P['same_run'] = P.run_a == P.run_b; P['same_reference_set'] = P.ref_a == P.ref_b; P['r_target_cq'] = r_cq[iu17]; P['r_dcq'] = r_dcq[iu17]; P['r_reference_vectors'] = r_ref[iu17]; P['r_expected_shared_reference_only'] = r_exp[iu17]
save(P, 'tbl_2_4_pairs.csv'); rc8 = des.rc_resid


def summ(mask, label):
    s = P[mask]
    if len(s) == 0: return dict(scope=label, n_pairs=0)
    return dict(scope=label, n_pairs=int(len(s)), mean_r_target_cq=float(s.r_target_cq.mean()), mean_r_dcq=float(s.r_dcq.mean()), mean_r_expected_shared_reference_only=float(s.r_expected_shared_reference_only.mean()),
                mean_abs_r_target_cq=float(s.r_target_cq.abs().mean()), mean_abs_r_dcq=float(s.r_dcq.abs().mean()), n_abs_r_ge_0632_target_cq=int((s.r_target_cq.abs() >= rc8).sum()), n_abs_r_ge_0632_dcq=int((s.r_dcq.abs() >= rc8).sum()),
                wilcoxon_p_dcq_vs_target_cq=float(stats.wilcoxon(s.r_dcq, s.r_target_cq).pvalue) if len(s) >= 6 else float('nan'))


srows = [summ(P.same_run, 'all same-run pairs'), summ(~P.same_run, 'all between-run pairs')] + [summ((P.run_a == k) & (P.run_b == k), 'same run: ' + k) for k in sorted(set(runs))] + \
        [summ((P.ref_a == k) & (P.ref_b == k), 'same reference set: ' + k) for k in sorted(set(refs))] + [summ(P.same_run & ~P.same_reference_set, 'same run, different reference sets (Jun24 refA x refC)')]
t24 = save(pd.DataFrame(srows), 'tbl_2_4_shared_reference_summary.csv'); print(t24.round(3).to_string())
# ADD-ON: layout control for the 26 Aug run, whose rows are in ascending Cq order within group (rank pairing), versus a random pairing; simulated from the observed group means and SDs
aug = [G17.index(g) for g in G17 if runs[G17.index(g)] == 'Aug26']; B = 5000; rngs = np.random.default_rng(SEED); sims = {}
for layout in ('ascending rows (as in the sheet)', 'random row order'):
    T_ = np.empty((B, 12, len(aug))); R_ = np.empty((B, 12))
    for k, rows_ in enumerate(des.rows):
        m = T[np.ix_(rows_, aug)].mean(0); s = T[np.ix_(rows_, aug)].std(0, ddof=1); draw = rngs.normal(m, s, (B, 4, len(aug)))
        if layout.startswith('ascending'): draw = np.sort(draw, axis=1)
        T_[:, rows_, :] = draw; ref_obs = RF[rows_, aug[0]]
        if k == 1: R_[:, rows_] = ref_obs                                                                                    # DAPT reference as entered in the sheet (C10, D10 twice), same in every simulation
        else:
            rd = rngs.normal(ref_obs.mean(), ref_obs.std(ddof=1), (B, 4)); R_[:, rows_] = np.sort(rd, axis=1) if layout.startswith('ascending') else rd
    sims[layout] = dict(mean_r_target_cq=upper_mean(corr_stack(des.Q @ T_)), mean_r_dcq=upper_mean(corr_stack(des.Q @ (T_ - R_[:, :, None]))))
aug_obs = P[(P.run_a == 'Aug26') & (P.run_b == 'Aug26')]
OUT['2_4_addon_layout_control_Aug26'] = dict(note='ADD-ON (not in the protocol): simulated Aug26 data with independent genes, observed group means/SDs; compare with the observed mean within-group r of the 10 Aug26 pairs',
                                             observed=dict(mean_r_target_cq=float(aug_obs.r_target_cq.mean()), mean_r_dcq=float(aug_obs.r_dcq.mean())),
                                             simulated={k: dict(mean_r_target_cq=dict(mean=float(v['mean_r_target_cq'].mean()), q025=float(np.percentile(v['mean_r_target_cq'], 2.5)), q975=float(np.percentile(v['mean_r_target_cq'], 97.5))),
                                                                 mean_r_dcq=dict(mean=float(v['mean_r_dcq'].mean()), q025=float(np.percentile(v['mean_r_dcq'], 2.5)), q975=float(np.percentile(v['mean_r_dcq'], 97.5)),
                                                                                 p_ge_observed=float((v['mean_r_dcq'] >= aug_obs.r_dcq.mean()).mean()))) for k, v in sims.items()})
print('  Aug26 layout control:', json.dumps(clean(OUT['2_4_addon_layout_control_Aug26']['simulated']), indent=None)[:600], '| observed', OUT['2_4_addon_layout_control_Aug26']['observed'])
# sensitivity versions of the pooled analysis
var_rows = []; base_set = {frozenset(p) for p in base['consensus_pairs']}


def variant(kind):
    R2 = RF.copy()
    if kind == 'a':                                                                                                         # (a) per-replicate reference -> mean of that gene's reference within group
        for rows_ in des.rows: R2[rows_] = RF[rows_].mean(0)
    elif kind == 'b': R2[:] = RF.mean(0)                                                                                  # (b) -> run-wide mean of the reference (= no normalization; a per-gene constant does not change r)
    return np.exp2(-(T - R2))


for name, Xv in (('reported (primary)', X), ('(a) reference = run-and-group mean', variant('a')), ('(b) reference = run-wide mean (no normalization)', variant('b'))):
    ov = cc.observed_summary(Xv, grp, DEN, G17); cs = {frozenset(p) for p in ov['consensus_pairs']}; cnv = cc.consensus_null(Xv, grp, DEN, 10000, SEED, genes=G17); ctv = cc.coupling_test(Xv, grp, 10000, SEED)
    rsv = cc.run_structure_test(Xv, grp, runs, 20000, SEED, G17); hubs_v = [g for g, k in ov['degree'].items() if k >= 5]; sh = cc.group_mean_share(Xv, grp, np.triu(ov['_S'][0], 1))
    rn = lambda p: '-'.join(sorted(p)) + ' (' + '/'.join(sorted({runs[G17.index(g)] for g in p})) + ')'; Q0 = ov['_Q'][0]
    extra = {'%s_%s' % (s, k): v[s2] for k, v in rsv['r']['per_run'].items() for s, s2 in (('mean_inside', 'mean_inside'), ('D', 'D'), ('p', 'p_two_sided'))}
    var_rows.append(dict(version=name, consensus_pairs=ov['consensus_any'], kept=len(cs & base_set), lost=len(base_set - cs), new=len(cs - base_set), null_mean=cnv['null_mean'], p_null_ge=cnv['p_ge'], pooled_raw=ov['pooled_sig_per_scheme'][0], pooled_Ins1=ov['pooled_sig_per_scheme'][1],
                         pooled_Kcnj11=ov['pooled_sig_per_scheme'][2], pooled_Hes1=ov['pooled_sig_per_scheme'][3], adj_raw=ov['adj_q_lt_alpha_per_scheme'][0], adj_Ins1=ov['adj_q_lt_alpha_per_scheme'][1], adj_Kcnj11=ov['adj_q_lt_alpha_per_scheme'][2],
                         adj_Hes1=ov['adj_q_lt_alpha_per_scheme'][3], hubs_deg_ge_5='; '.join('%s %d' % (g, ov['degree'][g]) for g in hubs_v), max_degree=max(ov['degree'].values()) if ov['degree'] else 0, T_r2=ctv['T_r2']['observed'], T_r2_null_mean=ctv['T_r2']['null_mean'],
                         T_r2_p=ctv['T_r2']['p_ge_plus1'], T_p=ctv['T_p']['observed'], T_p_null_mean=ctv['T_p']['null_mean'], T_p_p=ctv['T_p']['p_ge_plus1'], share_median=sh['median'], D_same_minus_between=rsv['r']['D'], D_p=rsv['r']['p_two_sided'],
                         mean_r_same_run=rsv['r']['mean_same'], mean_r_between_run=rsv['r']['mean_between'], **extra,
                         raw_scheme_pairs_q_lt_005='; '.join('%s-%s (%s) q=%.3f' % (G17[i], G17[j], '/'.join(sorted({runs[i], runs[j]})), Q0[i, j]) for i, j in zip(*iu17) if Q0[i, j] < 0.05),
                         consensus_pairs_lost='; '.join(sorted(rn(p) for p in base_set - cs)), consensus_pairs_new='; '.join(sorted(rn(p) for p in cs - base_set))))
t24v = save(pd.DataFrame(var_rows), 'tbl_2_4_sensitivity_versions.csv'); OUT['2_4_shared_reference'] = dict(summary=srows, versions=var_rows)
print(t24v[['version', 'consensus_pairs', 'kept', 'lost', 'new', 'null_mean', 'p_null_ge', 'adj_Ins1', 'adj_Kcnj11', 'adj_Hes1', 'hubs_deg_ge_5', 'T_r2', 'T_r2_p', 'D_same_minus_between', 'D_p']].round(3).to_string())

# ------------------------------------------------------------------------------------------------------------------ 2.5 calibrated adjusted counts
print('2.5 calibrated adjusted counts')
an = cc.adjusted_null(X, grp, DEN, 2000, SEED, genes=G17); ac = []
for n in an['scheme_names']:
    v = an['schemes'][n]; ac.append(dict(scheme=n, observed_q_lt_005=v['observed'], null_mean=v['null_mean'], null_sd=v['null_sd'], null_q95=v['null_q95'], excess=v['excess'], p_null_ge_obs=v['p_ge'], p_plus1=v['p_ge_plus1'], calibrated=(n == 'raw'),
                     reading='calibrated: no pair has q < 0.05 (null mean ~0)' if n == 'raw' else 'artefact-prone: shared denominator, observed count not above its null' if v['p_ge'] >= 0.05 else 'artefact-prone: shared denominator; observed above null at P < 0.05 (still not calibrated)'))
t25 = save(pd.DataFrame(ac), 'tbl_2_5_adjusted_calibrated.csv'); Q = base['_Q']
rq = P[['gene_a', 'gene_b', 'run_a', 'run_b', 'same_run']].copy(); rq['r_pooled_raw_linear'] = np.corrcoef(X.T)[iu17]; rq['r_adjusted_partial_linear'] = r_lin[iu17]
rq['p_adjusted'] = 2 * stats.t.sf(np.abs(rq.r_adjusted_partial_linear) * np.sqrt(des.df_resid / (1 - rq.r_adjusted_partial_linear ** 2)), des.df_resid)
for s, n in enumerate(NAMES): rq['q_' + n] = Q[s][iu17]
rq['consensus_pair'] = [bool(base['_cons_any'][i, j]) for i, j in zip(*iu17)]; save(rq.sort_values('q_raw'), 'tbl_2_5_q_per_pair.csv')
OUT['2_5_adjusted_calibrated'] = dict(schemes=ac, raw_scheme_min_q=float(np.nanmin(Q[0])), raw_scheme_pairs_q_lt_005=int((Q[0][iu17] < 0.05).sum()), smallest_raw_q_pairs=rq.sort_values('q_raw').head(5)[['gene_a', 'gene_b', 'r_adjusted_partial_linear', 'p_adjusted', 'q_raw']].to_dict('records'))
print(t25.round(3).to_string())

# ------------------------------------------------------------------------------------------------------------------ 2.6 within-group coupling and run-block null
print('2.6 coupling and run-block null')
sm, bm = cc.same_run_mask(runs); masks = dict(same_run=sm, between_run=bm)
ct_i = cc.coupling_test(X, grp, 10000, SEED, pair_masks=masks); ct_b = cc.coupling_test(X, grp, 10000, SEED, blocks=runs, pair_masks=masks)
cn_i = cc.consensus_null(X, grp, DEN, 10000, SEED, genes=G17); cn_b = cc.consensus_null(X, grp, DEN, 10000, SEED, blocks=runs, genes=G17)
ctrows = []
for nm, ct, nul in (('independent genes (protocol null)', ct_i, 'within-group permutation per gene'), ('run-block (replicate labels permuted jointly within each run)', ct_b, 'within-group permutation per run')):
    for k, v in ct.items(): ctrows.append(dict(null=nm, statistic=k, observed=v['observed'], null_mean=v['null_mean'], null_sd=v['null_sd'], null_q95=v['null_q95'], p_ge_plus1=v['p_ge_plus1']))
for nm, cn in (('independent genes (protocol null)', cn_i), ('run-block (replicate labels permuted jointly within each run)', cn_b)):
    ctrows.append(dict(null=nm, statistic='consensus_count', observed=cn['observed'], null_mean=cn['null_mean'], null_sd=cn['null_sd'], null_q95=cn['null_q95'], p_ge_plus1=cn['p_ge_plus1'], p_ge=cn['p_ge']))
    if nm.startswith('run-block'):                                                                                          # same-run pairs keep their correlation under this null, so only the between-run part is random
        ctrows.append(dict(null=nm, statistic='consensus_count_between_run_only', observed=cn['observed'] - a, null_mean=cn['null_mean'] - a, null_sd=cn['null_sd'], null_q95=cn['null_q95'] - a, p_ge_plus1=cn['p_ge_plus1'], p_ge=cn['p_ge']))
t26 = save(pd.DataFrame(ctrows), 'tbl_2_6_coupling.csv'); OUT['2_6_coupling'] = dict(independent=ct_i, run_block=ct_b, consensus_independent={k: v for k, v in cn_i.items() if not k.startswith('_')}, consensus_run_block={k: v for k, v in cn_b.items() if not k.startswith('_')},
                                                                                 consensus_same_run_pairs_observed=a, consensus_between_run_pairs_observed=c_)
print(t26.round(4).to_string())

# ------------------------------------------------------------------------------------------------------------------ 2.7 geometric expectation
print('2.7 geometry')
gz = {}
for nm, lg in (('ANOVA on log scale', True), ('ANOVA on linear scale', False)): gz[nm] = cc.geometry_check(X, grp, log=lg)
bm_ = pd.read_csv(os.path.join(BENCH, 'per_dataset_results.csv')); A = bm_[bm_['set'] == 'A'].copy(); A['share_raw'] = A.pooled_raw / 136.0
rcB = cc.rcrit(10); expB = cc.geometric_expectation(rcB, 3)
bsum = dict(n_series=int(len(A)), median_share=float(A.share_raw.median()), q25=float(A.share_raw.quantile(.25)), q75=float(A.share_raw.quantile(.75)), mean_share=float(A.share_raw.mean()), expected=expB, n_above_expected=int((A.share_raw > expB).sum()),
            by_platform={p: dict(n=int(len(g)), median=float(g.share_raw.median())) for p, g in A.groupby('platform')})
geo = pd.DataFrame([dict(data='INS-1, genes with ANOVA group effect (log scale, BH q<0.05): %d genes' % gz['ANOVA on log scale']['n_group_effect_genes'], n_pairs=gz['ANOVA on log scale']['n_pairs_group_effect'], n_sig=gz['ANOVA on log scale']['n_sig_group_effect'], share=gz['ANOVA on log scale']['share_sig_group_effect'], expected=gz['ANOVA on log scale']['expected']),
                    dict(data='INS-1, all 17 genes', n_pairs=136, n_sig=gz['ANOVA on log scale']['n_sig_all'], share=gz['ANOVA on log scale']['share_sig_all'], expected=gz['ANOVA on log scale']['expected']),
                    dict(data='INS-1, genes with ANOVA group effect (linear scale): %d genes' % gz['ANOVA on linear scale']['n_group_effect_genes'], n_pairs=gz['ANOVA on linear scale']['n_pairs_group_effect'], n_sig=gz['ANOVA on linear scale']['n_sig_group_effect'], share=gz['ANOVA on linear scale']['share_sig_group_effect'], expected=gz['ANOVA on linear scale']['expected']),
                    dict(data='benchmark part 1, gene set A (17 DE genes), median of %d series' % len(A), n_pairs=136, n_sig=float(A.pooled_raw.median()), share=bsum['median_share'], expected=expB)])
save(geo, 'tbl_2_7_geometry.csv'); OUT['2_7_geometry'] = dict(ins1=gz, benchmark_set_A=bsum, note='anchor, not a test; assumes the pooled correlation is carried by the group means with uniform orientation of the centred group-mean vectors')
print(geo.round(3).to_string())

# ------------------------------------------------------------------------------------------------------------------ 2.8 scale sensitivity
print('2.8 scale sensitivity'); srows8 = []
for kind, lab in (('pearson', 'linear-scale Pearson (the manuscript)'), ('log', 'log-scale Pearson'), ('spearman', 'Spearman')):
    ov = cc.observed_summary(X, grp, DEN, G17, kind=kind); cnk = cc.consensus_null(X, grp, DEN, 10000, SEED, kind=kind, genes=G17); sh = cc.group_mean_share(X, grp, np.triu(ov['_S'][0], 1), kind=kind)
    srows8.append(dict(scale=lab, pooled_raw=ov['pooled_sig_per_scheme'][0], pooled_Ins1=ov['pooled_sig_per_scheme'][1], pooled_Kcnj11=ov['pooled_sig_per_scheme'][2], pooled_Hes1=ov['pooled_sig_per_scheme'][3], consensus_pairs=ov['consensus_any'], same_sign=ov['consensus_same_sign'],
                       null_mean=cnk['null_mean'], null_sd=cnk['null_sd'], p_null_ge=cnk['p_ge'], covariance_share_median=sh['median'], covariance_share_n_pairs=sh['n_pairs'], adj_q_raw=ov['adj_q_lt_alpha_per_scheme'][0], adj_q_Ins1=ov['adj_q_lt_alpha_per_scheme'][1],
                       adj_q_Kcnj11=ov['adj_q_lt_alpha_per_scheme'][2], adj_q_Hes1=ov['adj_q_lt_alpha_per_scheme'][3], consensus_with_adj_q_any=ov['cons_adj_any'], hubs_deg_ge_5='; '.join('%s %d' % (g, k) for g, k in ov['degree'].items() if k >= 5)))
t28 = save(pd.DataFrame(srows8), 'tbl_2_8_scale_sensitivity.csv'); OUT['2_8_scale_sensitivity'] = srows8; print(t28.round(3).to_string())

OUT['runtime_seconds'] = round(time.time() - t_start, 1); json.dump(clean(OUT), open(os.path.join(HERE, 'results_ins1.json'), 'w'), indent=1); print('done in %.0fs; wrote results_ins1.json and tbl_*.csv' % (time.time() - t_start))
