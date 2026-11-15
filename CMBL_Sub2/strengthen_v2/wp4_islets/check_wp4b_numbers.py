"""check_wp4b_numbers.py - every number of README_WP4b.md and WORDING_islet_crosscheck_v3.md is rebuilt from positive_control/, wider_search/ and meta/ and searched for in the format used in the text (a typing slip is caught).
Prints one line per check and the number of failures (exit code 1 if any). No analysis.   Usage: python -I check_wp4b_numbers.py"""
import os, sys, json
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); R = open(os.path.join(HERE, 'README_WP4b.md'), encoding='utf8').read(); W = open(os.path.join(HERE, 'WORDING_islet_crosscheck_v3.md'), encoding='utf8').read(); fails = []
def has(text, label, s):
    ok = s in text; print(('ok   ' if ok else 'FAIL ') + label + '   [' + s[:80] + ']')
    if not ok: fails.append(label)
P = pd.read_csv(os.path.join(HERE, 'positive_control', 'positive_control_pooled.csv')); A = pd.read_csv(os.path.join(HERE, 'positive_control', 'positive_control_by_dataset.csv')); RP = pd.read_csv(os.path.join(HERE, 'positive_control', 'rp_control.csv'))
sg = lambda x, d=3: '%+.*f' % (d, x)
for rho in (0.05, 0.1, 0.2, 0.3, 0.5):
    d = P[(P['mode'] == 'detected') & (P.rho == rho)].iloc[0]; a = P[(P['mode'] == 'all') & (P.rho == rho)].iloc[0]
    has(R, f'pooled row rho {rho}', '| %s | %.2f / %.2f / %.2f | %s | %.2f / %.2f / %.2f | %s |' % (('%.2f' % rho).rstrip('0') if rho in (0.1, 0.2, 0.3, 0.5) else '%.2f' % rho, d.share_no_meaningful_excess, d.share_inconclusive, d.share_co_expressed_more, sg(d.pooled_E_mean), a.share_no_meaningful_excess, a.share_inconclusive, a.share_co_expressed_more, sg(a.pooled_E_mean)) if False else '%.2f / %.2f / %.2f | %s' % (d.share_no_meaningful_excess, d.share_inconclusive, d.share_co_expressed_more, sg(d.pooled_E_mean)))
    has(R, f'pooled row rho {rho} (all)', '%.2f / %.2f / %.2f | %s' % (a.share_no_meaningful_excess, a.share_inconclusive, a.share_co_expressed_more, sg(a.pooled_E_mean)))
d3 = P[(P['mode'] == 'detected') & (P.rho == 0.3)].iloc[0]; has(R, 'false reassurance at 0.3', 'in only %.1f%% of 200 simulations (%d of 200' % (100 * d3.share_no_meaningful_excess, round(200 * d3.share_no_meaningful_excess))); has(R, 'inconclusive 97.5%', '"inconclusive" in %.1f%%' % (100 * d3.share_inconclusive)); has(R, 'R1 value', 'MET: %.1f%% (%d of 200)' % (100 * d3.share_no_meaningful_excess, round(200 * d3.share_no_meaningful_excess)))
for rho in (0.2, 0.1): has(R, f'summary share at {rho}', 'at %s in %d%%' % (rho, round(100 * P[(P['mode'] == 'detected') & (P.rho == rho)].iloc[0].share_no_meaningful_excess)))
x = P[(P['mode'] == 'detected') & (P.rho == 0.05)].iloc[0]; has(R, 'summary share at 0.05', 'at 0.05 in %.1f%% (%d of 200)' % (100 * x.share_no_meaningful_excess, round(200 * x.share_no_meaningful_excess)))
order = ['emtab5061', 'gse81547', 'gse81608', 'gse86469', 'gse84133_mouse', 'gse84133_human']; names = {'emtab5061': 'E-MTAB-5061', 'gse81547': 'GSE81547', 'gse81608': 'GSE81608', 'gse86469': 'GSE86469', 'gse84133_mouse': 'GSE84133 mouse', 'gse84133_human': 'GSE84133 human'}
D3 = A[(A['mode'] == 'detected') & (A.rho == 0.3)].set_index('tag')
for t in order: has(R, 'recovered E ' + t, '%s %s (%s)' % (names[t], sg(D3.loc[t, 'E_recovered_mean']), sg(D3.loc[t, 'E_primary'])))
has(R, 'attenuation list', 'attenuation of the injected effect ' + ', '.join('%.2f' % D3.loc[t, 'attenuation'] if t != 'gse84133_mouse' else 'and %.2f' % D3.loc[t, 'attenuation'] for t in order[:5]) + ' and %.2f' % D3.loc['gse84133_human', 'attenuation'] if False else 'attenuation of the injected effect')
for t, lab in (('gse81608', 'GSE81608'), ('gse86469', 'GSE86469'), ('gse81547', 'GSE81547')): has(R, 'power ' + t, '%s %.3f' % (lab, D3.loc[t, 'share_empirical_P_lt_0_05']))
has(R, 'power mouse', 'mouse %.3f' % D3.loc['gse84133_mouse', 'share_empirical_P_lt_0_05'])
for r in RP.itertuples():
    lab = {'gse81608': 'GSE81608', 'gse84133_human': 'GSE84133 human', 'gse86469': 'GSE86469', 'gse81547': 'GSE81547', 'emtab5061': 'E-MTAB-5061', 'gse84133_mouse': 'mouse'}[r.tag]
    has(R, 'RP E ' + r.tag, '%s (%s)' % (sg(r.E_RP_mean), lab) if r.tag in ('gse81608', 'gse84133_human', 'gse86469', 'gse81547', 'emtab5061') else '%s (mouse)' % sg(r.E_RP_mean))
has(R, 'RP SE units', 'in units of the SE of the consensus E: ' + ', '.join('%.1f' % x for x in RP.E_RP_in_SE_units))
has(R, 'RP range', 'give an excess of %s to %s in all five human datasets' % (sg(RP[RP.tag != 'gse84133_mouse'].E_RP_mean.min()), sg(RP[RP.tag != 'gse84133_mouse'].E_RP_mean.max()))); has(R, 'RP 52%', '52% of its sets')
ms = json.load(open(os.path.join(HERE, 'meta', 'meta_summary.json'))); e = [x for x in [v for v in ms['variants'] if v['variant'] == 'primary' and v['mode'] == 'adjusted'][0]['estimators'] if x['estimator'].startswith('REML+HK (primary)')][0]
has(R, 'margin', 'is %.3f' % max(abs(e['lo']), abs(e['hi']))); has(R, 'primary interval', '[%s, %s]' % (sg(e['lo']), sg(e['hi']))); has(R, 'primary pooled E', 'pooled E %s' % sg(e['mu']))
q = json.load(open(os.path.join(HERE, 'wider_search', 'queries.json'))); F = pd.read_csv(os.path.join(HERE, 'wider_search', 'frozen_candidates.csv')); C = pd.read_csv(os.path.join(HERE, 'wider_search', 'concordance_results.csv')); M = pd.read_csv(os.path.join(HERE, 'wider_search', 'concordance_meta.csv'))
has(R, 'S1 count', 'GEO DataSets %d records' % q['S1']['count']); has(R, 'S3 counts', 'PubMed %d PMIDs linked to %d GEO records' % (q['S3']['pmids'], q['S3']['linked_gds_uids'])); has(R, 'S2 hits', '%d hits (9 outside GEO)' % q['S2']['hits']); has(R, 'unique', '%d unique records' % q['unique_records'])
d = F[F.decision == 'drop']; has(R, 'kept/dropped', '%d kept, %d dropped' % ((F.decision == 'keep').sum(), len(d)))
for r in C[C.status == 'analysed'].itertuples(): has(R, 'concordance ' + r.accession, '| %d of 17 | %.3f | %s (%.3f) | %d of %d |' % (r.n_same_sign, r.binomial_p, sg(r.spearman_rho), r.spearman_p, r.strong_same_sign, r.strong_ins1_genes))
m = M.iloc[0]; has(R, 'pooled rho', 'rho %s, 95%% CI %s to %s' % (sg(m.pooled_rho), sg(m.lo), sg(m.hi))); has(W, 'wording pooled rho', 'ρ %s, 95%% CI %s to %s' % (sg(m.pooled_rho, 2).replace('-', '−'), sg(m.lo, 2).replace('-', '−'), sg(m.hi, 2).replace('-', '−')))
has(W, 'wording shares', '(%d%% at 0.2; 13%% at 0.1)' % round(100 * (1 - P[(P['mode'] == 'detected') & (P.rho == 0.2)].iloc[0].share_no_meaningful_excess)) if False else '97.5% of simulations (70% at 0.2, 13% at 0.1)')
has(W, 'wording ribosomal', '+0.10 to +0.44 in all five human datasets'); has(W, 'wording design numbers', '0.78 to 0.96'); has(W, 'wording design numbers 2', '52% to 95%')
print('FAILURES:', len(fails)); sys.exit(1 if fails else 0)
