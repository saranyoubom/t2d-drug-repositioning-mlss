"""Port check of Fig_6d_S2_Validation_Top15.r (weight scenarios and leave-one-out ablation on the TOP 15 of the reported run) and of the 'Top 5 stability'
metric as the R code defines it.  Usage: python -I fig6d_port_check.py   (milliseconds; writes fig6d_top15_port_check.csv and fig6d_port_check_log.txt here).
Compared with the R log (R_17genes/logs/Fig_6d_S2_Validation_Top15.log) and with the Additional file 4 sheets Fig_6d_ranking_correlations / Fig_6d_ablation_impact.
In the R code Top-5 stability is which(rank(x, ties.method = "first") <= 5); rank() is ascending, so this selects the FIVE LOWEST values of the 15, not the five
highest; both versions are given.  No altered input of protocol section 5 is scored here: the reported table is only re-weighted exactly as the R script does."""
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

dc = M.read_drugcentral(); hubs, _ = M.hubs_from_centrality(min_degree=5)
REP = M.mlss_table(M.targets_from(dc, hubs), M.potency_from(dc), M.corr_from_matrix_csv())
WS = [('MLSS (Standard)', .45, .30, .10, .15), ('C-Heavy (60:25:15)', .60, .25, .10, .05), ('C-Heavy (55:20:25)', .55, .20, .15, .10), ('B-Heavy (35:50:15)', .35, .50, .10, .05),
      ('V-Heavy (35:30:35)', .35, .30, .30, .05), ('P-Heavy (30:20:10:40)', .30, .20, .10, .40), ('Balanced (40:40:10:10)', .40, .40, .10, .10), ('Network-Focus (30:20:40:10)', .30, .20, .40, .10)]
log('Fig_6d_S2 on the top 15 of the reported run (port) against Fig_6d_S2_Validation_Top15.log and the Additional file 4 sheets')
top = REP.head(15)
log('  components inside the top 15: C %s, B %s, V %s (distinct values); P from %.4f to %.4f; Synergy %s' % (sorted(set(top['Complementarity_C'].round(9))), sorted(set(top['Balance_B'].round(9))),
    sorted(set(top['Network_Coverage_V'].round(9))), top['Potency_P'].min(), top['Potency_P'].max(), sorted(set(top['Synergy_Bonus'].round(6)))))
S = {'Original_MLSS': top['MLSS_v4.0'].values}
for nm, c, b, v, p in WS:
    S[nm] = (c * top['Complementarity_C'] + b * top['Balance_B'] + v * top['Network_Coverage_V'] + p * top['Potency_P']).values * 9 + top['Synergy_Bonus'].values - top['Antagonism_Penalty'].values
RK = pd.DataFrame({k: stats.rankdata(v) for k, v in S.items()}); CM = RK.corr(method='pearson')
sc8 = [k for k in S if k != 'Original_MLSS']; pairs = [(a, b) for i, a in enumerate(sc8) for b in sc8[i + 1:]]
log('  pairwise Spearman among the 8 weight scenarios (%d pairs): min %.6f, max %.6f; Original_MLSS vs each: min %.6f' % (
    len(pairs), min(CM.loc[a, b] for a, b in pairs), max(CM.loc[a, b] for a, b in pairs), CM['Original_MLSS'].min()))
log('  (R log: Spearman with the original 1.0000 for all eight; Additional file 4 matrix all 1)')
log('  mean MLSS of the top 15 per scenario (R log: 9.943, 9.955, 9.949, 9.955, 9.955, 9.915, 9.949, 9.949): ' + ', '.join('%.3f' % v.mean() for k, v in S.items() if k != 'Original_MLSS'))

def low5(x):    # R: which(rank(x, ties.method = "first") <= 5): rank() is ascending -> the five LOWEST values (ties by position)
    return set(np.argsort(x, kind='stable')[:5])
def high5(x):   # the five largest values (ties by position)
    return set(np.argsort(-x, kind='stable')[:5])
o = top
mk = lambda w: (w[0] * o['Complementarity_C'] + w[1] * o['Balance_B'] + w[2] * o['Network_Coverage_V'] + w[3] * o['Potency_P']).values * 9 + o['Synergy_Bonus'].values - o['Antagonism_Penalty'].values
ab = {'Without_C': mk((0, .30 / .55, .10 / .55, .15 / .55)), 'Without_B': mk((.45 / .70, 0, .10 / .70, .15 / .70)), 'Without_V': mk((.45 / .90, .30 / .90, 0, .15 / .90)),
      'Without_P': mk((.45 / .85, .30 / .85, .10 / .85, 0)),
      'Base_Only': ((.45 * o['Complementarity_C'] + .30 * o['Balance_B'] + .10 * o['Network_Coverage_V'] + .15 * o['Potency_P']) * 9).values}
orig = top['MLSS_v4.0'].values; AR = []
for nm, x in ab.items():
    shift = float(np.abs(stats.rankdata(orig) - stats.rankdata(x)).mean()); drop = float((orig - x).mean())
    t5c = len(low5(orig) & low5(x)) / 5; t5i = len(high5(orig) & high5(x)) / 5
    AR.append(dict(component=nm, avg_rank_shift=shift, avg_score_drop=drop, top5_as_coded_five_lowest=t5c, top5_five_largest=t5i, n_distinct_scores_after_ablation=len(set(np.round(x, 9)))))
    log('  %-10s avg rank shift %.1f, avg score drop %+.9f, Top5_Stability as coded (five LOWEST) %.1f, with the five largest %.1f; distinct scores among the 15 after ablation: %d' % (
        nm, shift, drop, t5c, t5i, len(set(np.round(x, 9)))))
log('  R log / Additional file 4: rank shift 0, 0, 0, 3.6, 0; score drop 0.013998273, 0.007332429, 0.001901, -0.017109, 0.960356949; Top5_Stability 1, 1, 1, 0, 1')
pd.DataFrame(AR).to_csv(os.path.join(HERE, 'fig6d_top15_port_check.csv'), index=False)
with open(os.path.join(HERE, 'fig6d_port_check_log.txt'), 'w', encoding='utf8') as f:
    f.write('\n'.join(LOG) + '\n')
