"""WP3 / protocol section 5: recompute the MLSS with the independent port (mlss_port.py) under the pre-specified scenarios and report what
drives the ranking.  Usage: python -I mlss_scenarios.py   (CPU-bound, a few minutes; memory use is small, about 200-300 MB; writes mlss_scenarios.csv, mlss_scenario_scores.csv,
mlss_permutation_draws.csv and mlss_scenarios_log.txt next to this file after every family of scenarios, so a killed run keeps what it finished;
nothing outside this folder is written).  The first run on 4 Nov 2026 was reaped by the harness for system memory pressure before it wrote anything.
Scenarios of protocol section 5:
  gapdh_C0 / gapdh_C1  the 26 August DAPT reference recalculated (reanalysis_17genes/gapdh_exact.py) and carried through the four
                       schemes, consensus pairs, degree, hubs (degree >= 5), drug targets, Hes1-normalized r lookup and scores;
  r_*                  r of the Synergy Bonus/Antagonism Penalty from Kcnj11-, Ins1-normalized or raw pooled r, from the raw
                       treatment-adjusted r (linear and log scale, with and without the 'p < 0.05 only' filter of the R lookup), or none;
                       '_all' variants also feed V from the replaced lookup (the R script uses one lookup for V, Synergy, Antagonism);
  hub_ge4/5/6          hub rule on the reported degrees (Fig_6b_Gene_Centrality_Summary.csv);
  perm_potency         potency values permuted across the 18 drugs (missing values travel with the permutation), 1,000 draws, seed 2026;
  perm_target          drug-to-target mapping permuted across drugs (each drug keeps its potency), 1,000 draws, seed 2026.
Extras, labelled as such (not in the protocol): r_sweep (r of Cacna1d-Glp1r set to fixed values), r_any_pair (that r taken from each of the
120 Hes1-normalized pairs), weights_all153 / ablation_all153 (the Fig_6d_S2 weight scenarios and leave-one-out ablation on all 153
combinations; the R script uses the top 15 only) and a reproduction of the Fig_6d numbers on the top 15 for the port check.
Endpoints are computed against the reported scores in Fig_6c_MLSS_ALL_Combinations.csv (see column notes in README.md)."""
import os, sys, io
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np, pandas as pd
from scipy import stats
import mlss_port as M
sys.stdout.reconfigure(encoding='utf8')

SEED, NPERM = 2026, 1000
LOG = []
def log(s=''):
    LOG.append(str(s)); print(s)

# ------------------------------------------------------------------ inputs and the reported reference
dc = M.read_drugcentral(); POT = M.potency_from(dc)
HUBS5, HUB_DF = M.hubs_from_centrality(min_degree=5)
CORR_REP = M.corr_from_matrix_csv()
TG_REP = M.targets_from(dc, HUBS5)
X, SRC = M.load_matrix(); NET_REP = M.network_from_matrix(X)
R_CSV = pd.read_csv(os.path.join(M.RDIR, 'Fig_6c_MLSS_ALL_Combinations.csv'))
fs = lambda a, b: frozenset((a, b))
REP_SCORE = {fs(a, b): s for a, b, s in zip(R_CSV['Drug_A'], R_CSV['Drug_B'], R_CSV['MLSS_v4.0'])}
_R15 = R_CSV.sort_values('Rank').iloc[:15]
REP_TOP15 = {fs(a, b) for a, b in zip(_R15['Drug_A'], _R15['Drug_B'])}
REP = M.mlss_table(TG_REP, POT, CORR_REP)
assert np.abs(REP['MLSS_v4.0'].values - R_CSV.sort_values('Rank')['MLSS_v4.0'].values).max() < 1e-9      # the gate, repeated


def keys_of(T):
    return [fs(a, b) for a, b in zip(T['Drug_A'], T['Drug_B'])]


def rank_r2(y, cols):
    """R^2 of the average ranks of y regressed on the given columns (plus intercept)."""
    yr = stats.rankdata(y); A = np.column_stack([np.ones(len(y))] + [np.asarray(c, float) for c in cols])
    beta = np.linalg.lstsq(A, yr, rcond=None)[0]; res = yr - A @ beta
    return float(1 - (res ** 2).sum() / ((yr - yr.mean()) ** 2).sum()) if ((yr - yr.mean()) ** 2).sum() > 0 else np.nan


def safe_rho(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(stats.spearmanr(a, b)[0]) if len(a) > 2 and np.ptp(a) > 0 and np.ptp(b) > 0 else np.nan


def endpoints(T):
    """Endpoints of a scenario table T against the reported scores (see README for definitions)."""
    if len(T) == 0:
        return {}
    keys = keys_of(T); sc = T['MLSS_v4.0'].values; s9 = np.round(sc, 9)
    com = [i for i, k in enumerate(keys) if k in REP_SCORE]
    ref = np.array([REP_SCORE[keys[i]] for i in com])
    top15 = set(keys[:15])
    genes = sorted({g for col in ('Targets_A', 'Targets_B') for v in T[col] for g in v.split('|') if g})
    c1 = T['Complementarity_C'].values >= 1 - 1e-12
    return dict(n_drugs=len(set(T['Drug_A']) | set(T['Drug_B'])), n_combinations=len(T), n_targets_involved=len(genes), targets_involved='|'.join(genes),
                n_common=len(com), spearman_vs_reported=safe_rho(sc[com], ref), top15_overlap=len(top15 & REP_TOP15),
                top15_tie_at_cutoff=bool(len(sc) > 15 and s9[14] == s9[15]), top_pair=T['Combination'].iloc[0], n_tied_at_top=int((s9 == s9[0]).sum()),
                top_score=float(sc[0]), top_pair_in_reported_top15=bool(keys[0] in REP_TOP15),
                max_abs_score_change=float(np.abs(sc[com] - ref).max()) if com else np.nan,
                r2_potency_all=safe_rho(sc, T['Mean_pIC50'].values) ** 2,
                r2_potency_within_C1=safe_rho(sc[c1], T['Mean_pIC50'].values[c1]) ** 2 if c1.sum() > 2 else np.nan,
                r2_C_plus_potency=rank_r2(sc, [T['Complementarity_C'].values, stats.rankdata(T['Mean_pIC50'].values)]))


ROWS, SCORES = [], {}
def add(sid, family, desc, T, hubs, corr, extra=None):
    e = endpoints(T)
    row = dict(scenario_id=sid, family=family, description=desc, n_hubs=len(hubs), hubs='|'.join(sorted(hubs)),
               r_key_pair=corr.get(M.KEY, np.nan) if corr is not None else np.nan, **e)
    if extra:
        row.update(extra)
    ROWS.append(row)
    SCORES[sid] = dict(zip(keys_of(T), T['MLSS_v4.0'].values))
    log('%-34s rho %s  top15 overlap %s  top %s (%s tied)  targets %s  drugs %s' % (sid, '%.4f' % e['spearman_vs_reported'] if not pd.isna(e.get('spearman_vs_reported', np.nan)) else 'nan',
        e['top15_overlap'], e['top_pair'], e['n_tied_at_top'], e['n_targets_involved'], e['n_drugs']))


DRAWS = []
COLS = ['scenario_id', 'family', 'description', 'n_hubs', 'hubs', 'n_drugs', 'n_combinations', 'n_targets_involved', 'targets_involved', 'r_key_pair', 'n_common', 'spearman_vs_reported',
        'top15_overlap', 'top15_tie_at_cutoff', 'top_pair', 'n_tied_at_top', 'top_score', 'top_pair_in_reported_top15', 'max_abs_score_change', 'r2_potency_all', 'r2_potency_within_C1', 'r2_C_plus_potency',
        'n_draws', 'spearman_p2_5', 'spearman_p97_5', 'spearman_min', 'spearman_max', 'top15_overlap_p2_5', 'top15_overlap_p97_5', 'top15_overlap_max', 'top_pair_in_reported_top15_frac',
        'r2_potency_all_p2_5', 'r2_potency_all_p97_5']
def write_outputs():
    """Checkpoint: everything finished so far is written after every family, so a run that is killed keeps its results."""
    OUT = pd.DataFrame(ROWS).reindex(columns=COLS)
    OUT.to_csv(os.path.join(HERE, 'mlss_scenarios.csv'), index=False, float_format='%.6g')
    if DRAWS:
        pd.concat(DRAWS).to_csv(os.path.join(HERE, 'mlss_permutation_draws.csv'), index=False, float_format='%.6g')
    W = R_CSV[['Rank', 'Combination']].copy(); W['reported_MLSS'] = R_CSV['MLSS_v4.0'].values
    for sid, d in SCORES.items():
        W[sid] = [d.get(fs(a, b), np.nan) for a, b in zip(R_CSV['Drug_A'], R_CSV['Drug_B'])]
    W.to_csv(os.path.join(HERE, 'mlss_scenario_scores.csv'), index=False, float_format='%.9g')
    with open(os.path.join(HERE, 'mlss_scenarios_log.txt'), 'w', encoding='utf8') as f:
        f.write('\n'.join(LOG) + '\n')
    return len(OUT)


# ------------------------------------------------------------------ 0 reported
add('S00_reported_port', 'gate', 'Port with the reported inputs (reproduces the R output; gate in gate_output.txt)', REP, HUBS5, CORR_REP)

# ------------------------------------------------------------------ 1 Gapdh variants carried through the pipeline
for w in ('C0', 'C1'):
    net = M.network_from_matrix(M.gapdh_variant(X, w)); hubs = [g.upper() for g, d in net['degree'].items() if d >= 5]
    corr = M.hes1_lookup(net); tg = M.targets_from(dc, hubs)
    add('S01_gapdh_%s' % w, 'gapdh', 'Gapdh %s: 26 Aug DAPT reference %s; consensus %d pairs, hubs (degree >= 5) %s; Hes1-normalized r lookup from the recalculated matrix' % (
        w, 'mean of the two valid wells' if w == 'C0' else 'mean of all four wells', len(net['consensus']), sorted(hubs)), M.mlss_table(tg, POT, corr), hubs, corr)

# ------------------------------------------------------------------ 2 source of r in the Synergy Bonus / Antagonism Penalty
net_r = NET_REP
sources = [('kcnj11', 'Kcnj11-normalized pooled r', M.hes1_lookup(net_r, source='Kcnj11-normalized')),
           ('ins1', 'Ins1-normalized pooled r', M.hes1_lookup(net_r, source='Ins1-normalized')),
           ('raw_pooled', 'raw pooled r', M.hes1_lookup(net_r, source='Raw')),
           ('raw_adj_linear', 'raw treatment-adjusted r (linear scale, df 8), kept only if p < 0.05 as in the R lookup', M.adjusted_lookup(X, log=False, signif_only=True)),
           ('raw_adj_linear_nofilter', 'raw treatment-adjusted r (linear scale, df 8), value used whatever its p', M.adjusted_lookup(X, log=False, signif_only=False)),
           ('raw_adj_log', 'raw treatment-adjusted r (log scale, df 8; protocol 2.2 scale), kept only if p < 0.05', M.adjusted_lookup(X, log=True, signif_only=True)),
           ('raw_adj_log_nofilter', 'raw treatment-adjusted r (log scale, df 8), value used whatever its p', M.adjusted_lookup(X, log=True, signif_only=False)),
           ('none', 'no correlation in the Synergy Bonus and Antagonism Penalty (V kept as reported)', {})]
for name, desc, lk in sources:
    T = M.mlss_table(TG_REP, POT, lk, corr_v=CORR_REP)
    add('S02_r_%s' % name, 'r_source', 'r in Synergy/Antagonism: ' + desc + '; r(Cacna1d,Glp1r) = %s' % ('absent' if M.KEY not in lk else '%.4f' % lk[M.KEY]), T, HUBS5, lk)
    T2 = M.mlss_table(TG_REP, POT, lk)
    if np.abs(T2.set_index(['Drug_A', 'Drug_B'])['MLSS_v4.0'].sort_index().values - T.set_index(['Drug_A', 'Drug_B'])['MLSS_v4.0'].sort_index().values).max() > 1e-12:
        add('S02_r_%s_all' % name, 'r_source', 'as above but the same lookup also feeds V (as in the R script; an empty lookup makes V = 0 for every pair)', T2, HUBS5, lk)

# ------------------------------------------------------------------ 3 hub rule
for k in (4, 5, 6):
    hubs, _ = M.hubs_from_centrality(min_degree=k); tg = M.targets_from(dc, hubs)
    add('S03_hub_ge%d' % k, 'hub_rule', 'Hub rule degree >= %d on the reported network: hubs %s' % (k, sorted(hubs)), M.mlss_table(tg, POT, CORR_REP), hubs, CORR_REP)

write_outputs()          # checkpoint: deterministic scenarios done

# ------------------------------------------------------------------ 4 permutations
DRUGS = sorted(TG_REP, key=M.r_sort_key); GENE_OF = {d: TG_REP[d][0] for d in DRUGS}
def perm_potency(rng):
    slots = [POT.get((d, GENE_OF[d])) for d in DRUGS]; idx = rng.permutation(len(DRUGS)); p = {}
    for i, d in enumerate(DRUGS):
        if slots[idx[i]] is not None:
            p[(d, GENE_OF[d])] = slots[idx[i]]
    return TG_REP, p
def perm_target(rng):
    genes = [GENE_OF[d] for d in DRUGS]; idx = rng.permutation(len(DRUGS)); tg, p = {}, {}
    for i, d in enumerate(DRUGS):
        tg[d] = [genes[idx[i]]]
        if (d, genes[i]) in POT:
            p[(d, genes[idx[i]])] = POT[(d, genes[i])]
    return tg, p

def run_perm(name, fn):
    rng = np.random.default_rng(SEED); out = []
    for b in range(NPERM):
        tg, p = fn(rng); T = M.mlss_table(tg, p, CORR_REP); e = endpoints(T)
        out.append(dict(family=name, draw=b, **{k: e[k] for k in ('spearman_vs_reported', 'top15_overlap', 'top_pair', 'n_tied_at_top', 'top_pair_in_reported_top15', 'r2_potency_all', 'r2_potency_within_C1', 'r2_C_plus_potency', 'n_targets_involved')}))
    D = pd.DataFrame(out); DRAWS.append(D)
    q = lambda c, p_: float(np.percentile(D[c], p_))
    top = D['top_pair'].value_counts()
    extra = dict(n_draws=NPERM, spearman_p2_5=q('spearman_vs_reported', 2.5), spearman_p97_5=q('spearman_vs_reported', 97.5), spearman_min=float(D['spearman_vs_reported'].min()),
                 spearman_max=float(D['spearman_vs_reported'].max()), top15_overlap_p2_5=q('top15_overlap', 2.5), top15_overlap_p97_5=q('top15_overlap', 97.5),
                 top15_overlap_max=int(D['top15_overlap'].max()), top_pair_in_reported_top15_frac=float(D['top_pair_in_reported_top15'].mean()),
                 r2_potency_all_p2_5=q('r2_potency_all', 2.5), r2_potency_all_p97_5=q('r2_potency_all', 97.5))
    row = dict(scenario_id='S04_perm_%s' % name, family='permutation',
               description='%s permuted across the 18 drugs, %d draws, seed %d; spearman/overlap/r2 columns are means over draws, spread in the p2_5/p97_5 columns' % (
                   'potency values (missing values travel with the permutation)' if name == 'potency' else 'drug-to-target mapping (each drug keeps its potency)', NPERM, SEED),
               n_hubs=len(HUBS5), hubs='|'.join(sorted(HUBS5)), r_key_pair=CORR_REP[M.KEY], n_drugs=18, n_combinations=153, n_targets_involved=int(D['n_targets_involved'].mean()),
               targets_involved='CACNA1D|GLP1R', n_common=153, spearman_vs_reported=float(D['spearman_vs_reported'].mean()), top15_overlap=float(D['top15_overlap'].mean()),
               top_pair='most frequent: %s (%d of %d draws)' % (top.index[0], top.iloc[0], NPERM), n_tied_at_top=float(D['n_tied_at_top'].mean()),
               top_pair_in_reported_top15=bool(D['top_pair_in_reported_top15'].mean() > 0.5), r2_potency_all=float(D['r2_potency_all'].mean()),
               r2_potency_within_C1=float(D['r2_potency_within_C1'].mean()), r2_C_plus_potency=float(D['r2_C_plus_potency'].mean()), **extra)
    ROWS.append(row)
    log('%-34s rho mean %.3f (2.5-97.5%%: %.3f to %.3f; max %.3f)  top15 overlap mean %.2f (%.0f to %.0f; max %d)  top pair in reported top 15: %.1f%%  R2 potency all %.3f, within C=1 %.3f, C+potency %.3f' % (
        row['scenario_id'], row['spearman_vs_reported'], extra['spearman_p2_5'], extra['spearman_p97_5'], extra['spearman_max'], row['top15_overlap'],
        extra['top15_overlap_p2_5'], extra['top15_overlap_p97_5'], extra['top15_overlap_max'], 100 * extra['top_pair_in_reported_top15_frac'],
        row['r2_potency_all'], row['r2_potency_within_C1'], row['r2_C_plus_potency']))
run_perm('potency', perm_potency); write_outputs()
run_perm('target', perm_target); write_outputs()

# ------------------------------------------------------------------ extras (not in the protocol)
for r in (-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0):
    lk = dict(CORR_REP); lk[M.KEY] = r
    add('X1_r_sweep_%+.2f' % r, 'r_sweep', 'EXTRA: r(Cacna1d,Glp1r) set to %+.2f in the single lookup used for V, Synergy and Antagonism (r = 0 makes V = 0 for every pair)' % r, M.mlss_table(TG_REP, POT, lk), HUBS5, lk)
net_h = NET_REP['res']['Hes1-normalized']; pair_r = sorted(net_h.items(), key=lambda kv: kv[1][0]); out = []
for pr, (r, p) in net_h.items():
    lk = dict(CORR_REP); lk[M.KEY] = r; e = endpoints(M.mlss_table(TG_REP, POT, lk)); out.append((r, e['spearman_vs_reported'], e['top15_overlap'], e['top_pair_in_reported_top15']))
O = pd.DataFrame(out, columns=['r', 'rho', 'overlap', 'top_in'])
ROWS.append(dict(scenario_id='X2_r_any_pair_120', family='r_any_pair', description='EXTRA: r(Cacna1d,Glp1r) replaced, in turn, by each of the 120 Hes1-normalized pooled r (range %.3f to %.3f); columns are means over the 120 replacements, spread in min/max' % (O.r.min(), O.r.max()),
                 n_hubs=len(HUBS5), hubs='|'.join(sorted(HUBS5)), r_key_pair=float(O.r.mean()), n_drugs=18, n_combinations=153, n_targets_involved=2, targets_involved='CACNA1D|GLP1R', n_common=153,
                 spearman_vs_reported=float(O.rho.mean()), top15_overlap=float(O.overlap.mean()), n_draws=120, spearman_min=float(O.rho.min()), spearman_max=float(O.rho.max()),
                 spearman_p2_5=float(np.percentile(O.rho, 2.5)), spearman_p97_5=float(np.percentile(O.rho, 97.5)), top15_overlap_p2_5=float(np.percentile(O.overlap, 2.5)),
                 top15_overlap_p97_5=float(np.percentile(O.overlap, 97.5)), top15_overlap_max=int(O.overlap.max()), top_pair_in_reported_top15_frac=float(O.top_in.mean())))
log('%-34s rho mean %.4f (min %.4f, max %.4f), top15 overlap mean %.1f (min %d, max %d) over 120 replacements of r(Cacna1d,Glp1r) [%.3f to %.3f]' % (
    'X2_r_any_pair_120', O.rho.mean(), O.rho.min(), O.rho.max(), O.overlap.mean(), O.overlap.min(), O.overlap.max(), O.r.min(), O.r.max()))

WS = [('MLSS (Standard)', .45, .30, .10, .15), ('C-Heavy (60:25:15)', .60, .25, .10, .05), ('C-Heavy (55:20:25)', .55, .20, .15, .10), ('B-Heavy (35:50:15)', .35, .50, .10, .05),
      ('V-Heavy (35:30:35)', .35, .30, .30, .05), ('P-Heavy (30:20:10:40)', .30, .20, .10, .40), ('Balanced (40:40:10:10)', .40, .40, .10, .10), ('Network-Focus (30:20:40:10)', .30, .20, .40, .10)]
for i, (nm, c, b, v, p) in enumerate(WS):
    add('X3_weights_%d' % i, 'weights_all153', 'EXTRA: Fig_6d weight scenario "%s" (C:B:V:P) on all 153 combinations' % nm, M.mlss_table(TG_REP, POT, CORR_REP, weights=(c, b, v, p)), HUBS5, CORR_REP)
AB = [('Without_C', (0, .30 / .55, .10 / .55, .15 / .55)), ('Without_B', (.45 / .70, 0, .10 / .70, .15 / .70)), ('Without_V', (.45 / .90, .30 / .90, 0, .15 / .90)),
      ('Without_P', (.45 / .85, .30 / .85, .10 / .85, 0)), ('Base_Only', M.W_STD)]
for nm, w in AB:
    if nm == 'Base_Only':
        T = M.mlss_table(TG_REP, POT, {}, weights=w, corr_v=CORR_REP); lk = {}
    else:
        T = M.mlss_table(TG_REP, POT, CORR_REP, weights=w); lk = CORR_REP
    add('X4_ablation_%s' % nm, 'ablation_all153', 'EXTRA: Fig_6d leave-one-out "%s" (remaining weights renormalized; Base_Only = no Synergy/Antagonism) on all 153 combinations' % nm, T, HUBS5, lk)

# (the port check of Fig_6d on the top 15 and of its Top-5 metric is in fig6d_port_check.py)

# ------------------------------------------------------------------ final write
log(''); n_rows = write_outputs()
log('written: mlss_scenarios.csv (%d rows), mlss_permutation_draws.csv (%d rows), mlss_scenario_scores.csv, mlss_scenarios_log.txt' % (n_rows, sum(len(d) for d in DRAWS)))
write_outputs()
