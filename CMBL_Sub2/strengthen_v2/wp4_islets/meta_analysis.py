"""meta_analysis.py - protocol 4.4 (primary endpoint of WP4): random-effects pooling of the excess E of the consensus pairs over matched random pairs across the human public islet datasets, the pre-specified decision rule
(margin +/-0.10 on the E scale), leave-one-out and estimator sensitivities (PRESPEC_WP4_choices.md section B), and pair-level pooled rho for consensus pairs evaluable in at least 3 human datasets.
Estimators (fixed in advance): tau-squared by REML; Hartung-Knapp standard error with t(k-1) quantile and no floor = PRIMARY interval; sensitivity (never used for the rule): Knapp-Hartung factor floored at 1,
DerSimonian-Laird tau-squared with the Hartung-Knapp interval, DerSimonian-Laird with the normal interval, fixed effect. Datasets qualify with >= 3 evaluable consensus pairs; with < 2 qualifying datasets the rule is not applied.
Input: results/<tag>__<variant>__summary.json and results/<tag>__<variant>__adjusted__cons.csv from run_datasets.py. Output: meta/*.csv, meta/meta_summary.json, meta/meta_report.md.
Usage: python -I meta_analysis.py"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results'); OUT = os.path.join(HERE, 'meta'); os.makedirs(OUT, exist_ok=True)
HUMAN = ['gse81608', 'gse84133_human', 'gse86469', 'gse81547', 'emtab5061']
LABEL = {'gse81608': 'GSE81608 (Xin 2016)', 'gse84133_human': 'GSE84133 human (Baron 2016)', 'gse86469': 'GSE86469 (Lawlor 2017)', 'gse81547': 'GSE81547 (Enge 2017)', 'emtab5061': 'E-MTAB-5061 (Segerstolpe 2016)'}
MARGIN = 0.10


def reml_tau2(y, v, tol=1e-12, maxit=500):
    tau2 = max(0.0, float(np.var(y, ddof=1) - np.mean(v))) if len(y) > 1 else 0.0
    for _ in range(maxit):
        w = 1 / (v + tau2); mu = np.sum(w * y) / np.sum(w)
        new = max(0.0, float(np.sum(w ** 2 * ((y - mu) ** 2 - v)) / np.sum(w ** 2) + 1 / np.sum(w)))
        if abs(new - tau2) < tol: return new
        tau2 = new
    return tau2


def dl_tau2(y, v):
    w = 1 / v; mu = np.sum(w * y) / np.sum(w); Q = float(np.sum(w * (y - mu) ** 2)); k = len(y)
    return max(0.0, (Q - (k - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w)))


def pool(y, se, tau='reml', ci='hk', floor=False):
    """returns dict(mu, lo, hi, se, tau2, I2, Q, k); ci in {'hk', 'normal', 'fixed'}"""
    y, v = np.asarray(y, float), np.asarray(se, float) ** 2; k = len(y)
    w0 = 1 / v; mu0 = np.sum(w0 * y) / np.sum(w0); Q = float(np.sum(w0 * (y - mu0) ** 2)); I2 = max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0
    if ci == 'fixed': return dict(mu=float(mu0), se=float(np.sqrt(1 / np.sum(w0))), lo=float(mu0 - 1.959964 / np.sqrt(np.sum(w0))), hi=float(mu0 + 1.959964 / np.sqrt(np.sum(w0))), tau2=0.0, I2=I2, Q=Q, k=k)
    tau2 = reml_tau2(y, v) if tau == 'reml' else dl_tau2(y, v); w = 1 / (v + tau2); mu = float(np.sum(w * y) / np.sum(w))
    if ci == 'normal': se = float(np.sqrt(1 / np.sum(w))); q = 1.959964
    else:
        fac = float(np.sum(w * (y - mu) ** 2) / (k - 1)); fac = max(1.0, fac) if floor else fac; se = float(np.sqrt(fac / np.sum(w))); q = float(stats.t.ppf(0.975, k - 1))
    return dict(mu=mu, se=se, lo=mu - q * se, hi=mu + q * se, tau2=float(tau2), I2=I2, Q=Q, k=k)


def verdict(lo, hi, m=MARGIN):
    if lo > m: return 'consensus pairs are co-expressed more than random pairs'
    if lo >= -m and hi <= m: return 'no meaningful excess co-expression'
    return 'inconclusive'


def study_table(variant, mode):
    rows = []
    for tag in HUMAN:
        p = os.path.join(RES, f'{tag}__{variant}__summary.json')
        if not os.path.exists(p): continue
        s = json.load(open(p))[mode]
        rows.append(dict(tag=tag, label=LABEL[tag], variant=variant, mode=mode, beta_cells=s['beta_cells'], donors=s['donors'], genes_evaluable=s['genes_evaluable'], pairs=s['consensus_pairs_evaluable'], E=s.get('E'), SE=s.get('SE_E'),
                         median_rho=s.get('median_rho_consensus'), null_rho=s.get('null_median_rho_mean')))
    T = pd.DataFrame(rows); T['qualifies'] = T.pairs >= 3; return T


def meta_for(T):
    Q = T[T.qualifies]
    if len(Q) < 2: return None
    return Q


def run():
    out = {}; variants = []
    # primary: adjusted E, all frozen beta-cells, REML + HK
    for variant in ('primary', 'nondiabetic'):
        for mode in ('adjusted', 'unadjusted'):
            T = study_table(variant, mode)
            if variant == 'nondiabetic':      # datasets without T2D donors keep their primary result in the non-diabetic-only analysis
                P = study_table('primary', mode); T = pd.concat([T, P[~P.tag.isin(T.tag)].assign(variant='nondiabetic')], ignore_index=True)
                T = T.set_index('tag').loc[[t for t in HUMAN if t in set(T.tag)]].reset_index()
            T.to_csv(os.path.join(OUT, f'study_effects__{variant}__{mode}.csv'), index=False); Q = meta_for(T)
            rec = dict(variant=variant, mode=mode, qualifying=[] if Q is None else Q.tag.tolist(), n_qualifying=0 if Q is None else len(Q))
            if Q is not None:
                y, se = Q.E.to_numpy(float), Q.SE.to_numpy(float)
                for name, kw in (('REML+HK (primary)', dict(tau='reml', ci='hk', floor=False)), ('REML+HK, factor floored at 1', dict(tau='reml', ci='hk', floor=True)), ('DL+HK', dict(tau='dl', ci='hk')),
                                 ('DL, normal CI', dict(tau='dl', ci='normal')), ('fixed effect', dict(ci='fixed'))):
                    r = pool(y, se, **kw); r.update(estimator=name, verdict=verdict(r['lo'], r['hi'])); rec.setdefault('estimators', []).append(r)
                if len(Q) >= 3:                                                      # leave-one-out (needs >= 2 left)
                    rec['leave_one_out'] = []
                    for t in Q.tag:
                        R = Q[Q.tag != t]; r = pool(R.E.to_numpy(float), R.SE.to_numpy(float)); r.update(left_out=t, verdict=verdict(r['lo'], r['hi'])); rec['leave_one_out'].append(r)
            variants.append(rec)
    out['variants'] = variants
    # pair-level pooled rho (adjusted, primary, human datasets)
    pairs = {}
    for tag in HUMAN:
        p = os.path.join(RES, f'{tag}__primary__adjusted__cons.csv'); sp = os.path.join(RES, f'{tag}__primary__summary.json')
        if not (os.path.exists(p) and os.path.exists(sp)): continue
        df = json.load(open(sp))['adjusted']['df_for_p']; C = pd.read_csv(p)
        for r in C.itertuples(): pairs.setdefault((r.gene1, r.gene2), []).append((tag, r.rho, df))
    rows = []
    for (g1, g2), L in pairs.items():
        if len(L) < 3: continue
        z = np.array([np.arctanh(x[1]) for x in L]); se = np.array([np.sqrt(1.06 / (x[2] - 1)) for x in L]); r = pool(z, se)
        rows.append(dict(gene1=g1, gene2=g2, datasets=len(L), tags=';'.join(x[0] for x in L), rho_by_dataset=';'.join('%.3f' % x[1] for x in L), pooled_rho=np.tanh(r['mu']), lo=np.tanh(r['lo']), hi=np.tanh(r['hi']), tau2=r['tau2'], I2=r['I2']))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'pair_level_pooled_rho.csv'), index=False); out['pairs_evaluable_in_ge3'] = len(rows)
    json.dump(out, open(os.path.join(OUT, 'meta_summary.json'), 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    L = ['# WP4 meta-analysis (protocol 4.4)', '', 'Primary = depth-adjusted E, all frozen beta-cells, REML tau-squared + Hartung-Knapp 95%% CI (t, k - 1 df), margin +/-%.2f on the E (Fisher z) scale.' % MARGIN, '']
    for rec in variants:
        L.append('## %s, %s mode: %d qualifying dataset(s): %s' % (rec['variant'], rec['mode'], rec['n_qualifying'], ', '.join(rec['qualifying']) or '-'))
        if 'estimators' in rec:
            for r in rec['estimators']: L.append('* %-30s pooled E %+.3f  95%% CI [%+.3f, %+.3f]  tau2 %.4f  I2 %.0f%%  -> %s' % (r['estimator'], r['mu'], r['lo'], r['hi'], r['tau2'], 100 * r['I2'], r['verdict']))
            for r in rec.get('leave_one_out', []): L.append('  leave out %-15s pooled E %+.3f  [%+.3f, %+.3f] -> %s' % (r['left_out'], r['mu'], r['lo'], r['hi'], r['verdict']))
        else: L.append('* decision rule not applied (fewer than 2 qualifying datasets)')
        L.append('')
    open(os.path.join(OUT, 'meta_report.md'), 'w', encoding='utf8').write('\n'.join(L)); print('\n'.join(L))


if __name__ == '__main__': run()
