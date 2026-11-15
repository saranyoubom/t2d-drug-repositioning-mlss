"""make_wp2_figure.py - DRAFT figure for benchmark part 2 (WP2 deliverable F); WP6a refines it for the manuscript. Plots saved results only; no analysis.
(a) Power of the coupling tests and of the consensus-count test when one latent factor is injected into 6 genes of decoupled public data (3 groups x 4 replicates), by within-group correlation r;
    the leftmost points are the type I error without injection. (b) Power of T_r2 by the number of coupled genes.
(c) Share of public data sets in which T_r2 detects within-group coupling (P < 0.05), by arm and gene set, with Wilson 95% intervals; the dashed line is the 5% expected without coupling.
(d) Observed adjusted q < 0.05 count divided by its own null mean, by normalization scheme (gene set A, all analysed series; the line marks a count equal to its null).
(e, f) Design guidance by simulation (G = 17, one module of 6 genes, 3 groups): power of T_r2 and share of true module pairs detected, by number of independent experiments per group.
No key strip is drawn: curves are labelled at their ends, and the symbols are explained in the legend text. Usage: python -I make_wp2_figure.py"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'mathtext.default': 'regular'})
NAVY, LBLUE, SALMON, SLATE, GREY = '#3C5488', '#4DBBD5', '#F39B7F', '#8491B4', '#9A9A9A'
STAT_COL = {'T_r2': NAVY, 'T_p': LBLUE, 'T_eig': SLATE, 'count': SALMON}
STAT_LAB = {'T_r2': r'$T_{r2}$', 'T_p': r'$T_{p}$', 'T_eig': r'$T_{eig}$', 'count': 'count'}
STAT_DY = {'T_r2': 0.0, 'T_p': -0.075, 'T_eig': 0.075, 'count': 0.0}


def have(name): return os.path.exists(os.path.join(HERE, name))


def clean(ax):
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)


fig = plt.figure(figsize=(6.69, 4.9))
gs = fig.add_gridspec(2, 3, left=0.075, right=0.925, top=0.925, bottom=0.125, wspace=0.62, hspace=0.70, width_ratios=[1, 1, 1.3])


def tag(ax, s, dx=0.062):
    p = ax.get_position(); fig.text(p.x0 - dx, p.y1 + 0.062, s, fontsize=9, fontweight='bold', va='top')


XP = {0.0: 0.0, 0.3: 1.0, 0.5: 2.0, 0.7: 3.0, 0.9: 4.0}
XT = ['none', '0.3', '0.5', '0.7', '0.9']
# ---- (a) and (b) calibration
if have('calibration_tables.csv'):
    T = pd.read_csv(os.path.join(HERE, 'calibration_tables.csv'))
    ax = fig.add_subplot(gs[0, 0]); tag(ax, '(a)')
    for s in ('T_r2', 'T_p', 'T_eig', 'count'):
        d = pd.concat([T[(T.statistic == s) & (T.m == 0)], T[(T.statistic == s) & (T.m == 6)]]).sort_values('r'); x = np.array([XP[r] for r in d.r]); y = d.rate.values
        ax.errorbar(x, y, yerr=[y - d.ci_lo.values, d.ci_hi.values - y], color=STAT_COL[s], lw=1.1, marker='o', ms=3, capsize=1.5, elinewidth=0.6, zorder=3)
        ax.text(4.3, y[-1] + STAT_DY[s], STAT_LAB[s], color=STAT_COL[s], fontsize=7, va='center')
    ax.axhline(0.8, color='black', lw=0.5, ls='--', zorder=1); ax.axhline(0.05, color='black', lw=0.5, ls=':', zorder=1)
    ax.set_xticks(list(XP.values())); ax.set_xticklabels(XT); ax.set_xlim(-0.3, 5.5); ax.set_ylim(-0.02, 1.1); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel('Injected within-group r (6 genes)'); ax.set_ylabel('Share of tests with P < 0.05'); ax.set_title('Power and type I error', fontsize=7.5, fontweight='bold'); clean(ax)
    ax = fig.add_subplot(gs[0, 1]); tag(ax, '(b)'); null = T[(T.statistic == 'T_r2') & (T.m == 0)].iloc[0]
    for m, col, dy in ((3, LBLUE, 0.0), (6, NAVY, -0.06), (10, SALMON, 0.06)):
        d = T[(T.statistic == 'T_r2') & (T.m == m)].sort_values('r'); x = np.array([0.0] + [XP[r] for r in d.r]); y = np.r_[null.rate, d.rate.values]; lo = np.r_[null.ci_lo, d.ci_lo.values]; hi = np.r_[null.ci_hi, d.ci_hi.values]
        ax.errorbar(x, y, yerr=[y - lo, hi - y], color=col, lw=1.1, marker='o', ms=3, capsize=1.5, elinewidth=0.6, zorder=3); ax.text(4.3, y[-1] + dy, '%d genes' % m, color=col, fontsize=7, va='center')
    ax.axhline(0.8, color='black', lw=0.5, ls='--', zorder=1); ax.axhline(0.05, color='black', lw=0.5, ls=':', zorder=1)
    ax.set_xticks(list(XP.values())); ax.set_xticklabels(XT); ax.set_xlim(-0.3, 5.7); ax.set_ylim(-0.02, 1.1); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel('Injected within-group r'); ax.set_ylabel('Share of tests with P < 0.05'); ax.set_title(r'$T_{r2}$ by number of coupled genes', fontsize=7.5, fontweight='bold'); clean(ax)
else: print('missing calibration_tables.csv: panels (a), (b) skipped')
# ---- (c) prevalence of within-group coupling
lf = 'validation_endpoints_long.csv' if have('validation_endpoints_long.csv') else 'dev_endpoints_long.csv'
if have(lf):
    L = pd.read_csv(os.path.join(HERE, lf)); L = L[(L.endpoint == 'coupling_T_r2')]
    order = [('Development arrays (37)', 'Development\narrays'), ('Held-out GPL570', 'Held-out\nhuman arrays'), ('Held-out GPL1261', 'Held-out\nmouse arrays'), ('RNA-seq human', 'RNA-seq\nhuman')]
    if lf.startswith('dev_'): order = [('development, human GPL570', 'Development\nhuman arrays'), ('development, mouse GPL1261', 'Development\nmouse arrays')]
    ax = fig.add_subplot(gs[0, 2]); tag(ax, '(c)', dx=0.125); rows = []
    for key, lab in order:
        for gset, dy, mk in (('A', 0.17, 'o'), ('B', -0.17, 's')):
            r = L[(L.arm == key) & (L.gene_set == gset)]
            if len(r): rows.append((lab, gset, dy, mk, r.iloc[0]))
    labs = [lab for _, lab in order if any(x[0] == lab for x in rows)]; ypos = {lab: len(labs) - 1 - i for i, lab in enumerate(labs)}
    for lab, gset, dy, mk, r in rows:
        y = ypos[lab] + dy; col = NAVY if gset == 'A' else SALMON
        ax.plot([r.lo, r.hi], [y, y], color=col, lw=1.0, solid_capstyle='butt', zorder=2); ax.scatter([r.estimate], [y], s=14, marker=mk, color=col, zorder=3, linewidths=0)
        ax.text(1.04, y, '%d/%d' % (r.k, r.n), fontsize=6, va='center', color=col)
    ax.axvline(0.05, color='black', lw=0.5, ls='--', zorder=1); ax.set_xlim(0, 1.0); ax.set_ylim(-0.6, len(labs) - 0.4); ax.set_yticks([ypos[l] for l in labs]); ax.set_yticklabels(labs, fontsize=6.8, linespacing=1.1)
    ax.set_xlabel('Data sets with within-group coupling'); ax.set_title('Prevalence in public data', fontsize=7.5, fontweight='bold'); ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0]); ax.set_xticklabels(['0', '25%', '50%', '75%', '100%']); clean(ax)
else: print('missing endpoint tables: panel (c) skipped')
# ---- (d) adjusted counts against their own null
pf = 'validation_per_series.csv'
if have(pf):
    P = pd.read_csv(os.path.join(HERE, pf)); P = P[P.gene_set == 'A']; ax = fig.add_subplot(gs[1, 0]); tag(ax, '(d)'); rng = np.random.default_rng(1)
    for i, s in enumerate(('raw', 'den1', 'den2', 'den3')):
        o, nm = P['adjq_%s_obs' % s], P['adjq_%s_null_mean' % s]; r = (o / nm).where(nm > 0).dropna(); r = r.clip(lower=0.02, upper=50)
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(r)), r, s=5, color=NAVY if s != 'raw' else SALMON, alpha=0.55, linewidths=0, zorder=2)
        q1, md, q3 = r.quantile([.25, .5, .75]); ax.plot([i, i], [q1, q3], color='black', lw=1.0, zorder=4); ax.plot([i - 0.2, i + 0.2], [md, md], color='black', lw=1.2, zorder=4)
    ax.axhline(1, color='black', lw=0.5, ls='--', zorder=1); ax.set_yscale('log'); ax.set_ylim(0.015, 70); ax.set_xticks(range(4)); ax.set_xticklabels(['raw', 'ratio\nto 1', 'ratio\nto 2', 'ratio\nto 3'], fontsize=6.8)
    ax.set_ylabel('Adjusted count / null mean'); ax.set_title('Adjusted q < 0.05 counts', fontsize=7.5, fontweight='bold'); clean(ax)
else: print('missing validation_per_series.csv: panel (d) skipped')
# ---- (e, f) design table
if have('design_table.csv'):
    D = pd.read_csv(os.path.join(HERE, 'design_table.csv')); cols = {0.3: LBLUE, 0.5: SLATE, 0.7: NAVY, 0.9: SALMON}
    for k, (col, ttl, yl) in enumerate((('power_T_r2', r'Power of $T_{r2}$', 'Share of experiments with P < 0.05'), ('pair_power', 'True module pairs detected', 'Share detected at BH q < 0.05'))):
        ax = fig.add_subplot(gs[1, 1 + k]); tag(ax, '(e)' if k == 0 else '(f)', dx=0.062 if k == 0 else 0.125)
        for r, c in cols.items():
            d = D[D.r == r].sort_values('n_per_group'); ax.plot(np.log2(d.n_per_group), d[col], color=c, lw=1.1, marker='o', ms=3, zorder=3); ax.text(np.log2(30) + 0.15, d[col].iloc[-1] + {0.3: 0.0, 0.5: 0.0, 0.7: 0.0, 0.9: 0.0}[r] + (0.0 if d[col].iloc[-1] < 0.95 else {0.5: -0.14, 0.7: -0.05, 0.9: 0.04, 0.3: 0.13}[r]), 'r = %.1f' % r, color=c, fontsize=6.5, va='center')
        ax.axvline(np.log2(4), color='black', lw=0.5, ls=':', zorder=1)
        if k == 0: ax.axhline(0.8, color='black', lw=0.5, ls='--', zorder=1)
        ax.set_xticks(np.log2([4, 6, 8, 12, 20, 30])); ax.set_xticklabels(['4', '6', '8', '12', '20', '30'], fontsize=6.8); ax.set_xlim(np.log2(4) - 0.2, np.log2(30) + 1.2); ax.set_ylim(-0.02, 1.03)
        ax.set_xlabel('Independent experiments per group'); ax.set_ylabel(yl, fontsize=7); ax.set_title(ttl, fontsize=7.5, fontweight='bold'); clean(ax)
else: print('missing design_table.csv: panels (e), (f) skipped')
out = os.path.join(HERE, 'WP2_benchmark2_figure_draft'); fig.savefig(out + '.png', dpi=300); fig.savefig(out + '.pdf'); print('saved', out)
