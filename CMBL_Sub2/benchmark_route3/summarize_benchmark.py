"""Summarize results_per_dataset.jsonl into the pre-specified endpoints (PROTOCOL section 4 and 5) and one figure. Applies the interpretive bands fixed in the protocol; nothing is tuned here.
Usage: python -I summarize_benchmark.py"""
import sys, os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels.stats.proportion import proportion_confint
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__))
recs = [json.loads(l) for l in open(os.path.join(HERE, 'results_per_dataset.jsonl'), encoding='utf8')]
status = pd.DataFrame([dict(gse=r['gse'], platform=r['platform'], order=r['order'], status=r['status'], n_de=r.get('n_de')) for r in recs]); status.to_csv(os.path.join(HERE, 'processing_log.csv'), index=False)
ok = [r for r in recs if r['status'] == 'analysed']
rows = []
for r in ok:
    for tag in ('A', 'B', 'A_spike'):
        x = r[tag]; obs = x['observed']
        rows.append(dict(gse=r['gse'], platform=r['platform'], set=tag, scale=r['scale'], n_de=r['n_de'], n_de_in_set=x.get('n_de_in_set'), observed=obs, null_mean=x['null_mean'], null_sd=x['null_sd'], null_q95=x['null_q95'], z=x['z'], perm_p=x['perm_p_plain'],
                         observed_same_sign=x['observed_same_sign'], adj_any=x['cons_adj_any'], adj_all4=x['cons_adj_all4'], adj_pairs_all4=x['adj_pairs_all4'], retained_any=(x['cons_adj_any'] / obs if obs else np.nan),
                         retained_all4=(x['cons_adj_all4'] / obs if obs else np.nan), share_between=x['share_between_raw'], pooled_raw=x['pooled_sig_per_scheme'][0], adj_q05_raw=x['adj_q05_per_scheme'][0]))
T = pd.DataFrame(rows); T.to_csv(os.path.join(HERE, 'per_dataset_results.csv'), index=False)
S = {'n_series_logged': len(recs), 'n_analysed': len(ok), 'analysed_per_platform': {k: int(v) for k, v in status[status.status == 'analysed'].groupby('platform').size().items()},
     'not_analysable': status[status.status != 'analysed'].status.str.replace(r'\d+', 'N', regex=True).value_counts().to_dict()}


def prop(mask):
    k, n = int(mask.sum()), int(mask.size); lo, hi = proportion_confint(k, n, method='wilson') if n else (np.nan, np.nan); return dict(k=k, n=n, prop=(k / n if n else np.nan), ci95=[float(lo), float(hi)])


def med(x): x = pd.Series(x).dropna(); return dict(n=int(len(x)), median=float(x.median()) if len(x) else np.nan, q1=float(x.quantile(.25)) if len(x) else np.nan, q3=float(x.quantile(.75)) if len(x) else np.nan)


for tag in ('A', 'B', 'A_spike'):
    d = T[T.set == tag]
    S[tag] = dict(E1_not_above_null_P_ge_0_05=prop(d.perm_p >= 0.05), E1_flagged_P_lt_0_05=prop(d.perm_p < 0.05), z=med(d.z), observed=med(d.observed), null_mean=med(d.null_mean),
                  observed_minus_null=med(d.observed - d.null_mean), E2_retained_any_scheme=med(d.retained_any), E2_retained_all_schemes=med(d.retained_all4), E2_pairs_sig_all4_after_adjustment=med(d.adj_pairs_all4), E3_share_between=med(d.share_between),
                  datasets_with_zero_observed=int((d.observed == 0).sum()))
    for gpl in ('GPL570', 'GPL1261'):
        dd = d[d.platform == gpl]; S[tag]['by_' + gpl] = dict(n=int(len(dd)), P_ge_0_05=prop(dd.perm_p >= 0.05) if len(dd) else None, retained_any=med(dd.retained_any))
band = S['A']['E1_not_above_null_P_ge_0_05']['prop']
S['interpretive_band_gene_set_A'] = ('pooled consensus counts are generally explained by the design (>= 80%)' if band >= 0.8 else 'often explained (50-80%)' if band >= 0.5 else 'generally not explained: within-group coupling is common (< 50%)') if band == band else 'no analysable dataset'
S['positive_control_detected_P_lt_0_05'] = S['A_spike']['E1_flagged_P_lt_0_05']; S['positive_control_ok_ge_80pct'] = bool(S['A_spike']['E1_flagged_P_lt_0_05']['prop'] >= 0.8) if S['A_spike']['E1_flagged_P_lt_0_05']['n'] else None
json.dump(S, open(os.path.join(HERE, 'summary_endpoints.json'), 'w'), indent=1, default=float); print(json.dumps(S, indent=1, default=float))
# ---- figure: (a) observed vs null per dataset, (b) permutation P distributions, (c) share retained after adjustment
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.5, 'axes.linewidth': 0.6})
A = T[T.set == 'A'].sort_values('z').reset_index(drop=True); fig = plt.figure(figsize=(6.69, 7.0)); gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.38, wspace=0.32, left=0.09, right=0.985, top=0.965, bottom=0.075)
ax = fig.add_subplot(gs[0, :]); x = np.arange(len(A))
ax.vlines(x, A.null_mean - 1.96 * A.null_sd, A.null_mean + 1.96 * A.null_sd, color='#C9C9C9', lw=3, zorder=1); ax.scatter(x, A.null_mean, s=8, color='#8A8A8A', zorder=2, label='_')
ax.scatter(x, A.observed, s=16, color=np.where(A.perm_p < 0.05, '#E64B35', '#3C5488'), zorder=3, linewidths=0); ax.set_xticks(x); ax.set_xticklabels(A.gse, rotation=90, fontsize=5.5); ax.set_ylabel('Consensus pairs (pooled analysis)')
ax.set_title('Gene set A: 17 differentially expressed genes per dataset; dot = observed, grey bar = permutation null mean ± 1.96 SD', fontsize=7.5, loc='left'); ax.set_xlim(-0.8, len(A) - 0.2)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
ax = fig.add_subplot(gs[1, 0]); bins = np.linspace(0, 1, 11)
for tag, col, off in (('A', '#3C5488', -0.012), ('A_spike', '#E64B35', 0.012)):
    h, _ = np.histogram(T[T.set == tag].perm_p, bins=bins); ax.bar(bins[:-1] + 0.05 + off * 3, h, width=0.036, color=col)
ax.set_xlabel('Permutation P (null ≥ observed)'); ax.set_ylabel('Datasets'); ax.set_title('Observed (blue) and with injected within-group coupling (red)', fontsize=7.5, loc='left')
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
ax = fig.add_subplot(gs[1, 1]); d = A.dropna(subset=['retained_any']); ax.hist(d.retained_any, bins=np.linspace(0, 1, 11), color='#4DBBD5'); ax.set_xlabel('Consensus pairs with q < 0.05 after\ntreatment adjustment (share, ≥ 1 scheme)'); ax.set_ylabel('Datasets')
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
fig.text(0.005, 0.99, '(a)', fontsize=9, fontweight='bold', va='top'); fig.text(0.005, 0.40, '(b)', fontsize=9, fontweight='bold', va='top'); fig.text(0.52, 0.40, '(c)', fontsize=9, fontweight='bold', va='top')
fig.savefig(os.path.join(HERE, 'Benchmark_figure.png'), dpi=300); fig.savefig(os.path.join(HERE, 'Benchmark_figure.pdf')); print('figure saved')
