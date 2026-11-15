"""make_fig7_v2.py - DRAFT figures for the version-2 public islet cross-check (WP4). Plots saved results only; no new analysis. Written beside the earlier Fig. 7 script (islet_figures/make_islet_figure.py), same layout
and palette, NO key strips (symbols are explained in the legend text), panel tags (a), (b), ...
Figure 7 (updated, 170 x 190 mm): (a) detection of the 17 genes in the author-labelled beta-cells of six datasets (numbers = % of beta-cells; grey = below the 10% rule, not evaluable);
 (b) the 21 consensus pairs in the same six datasets: depth-adjusted partial Spearman rho (dot) and expression-matched random pairs (grey bar = null mean to 95th percentile); NE = not evaluable.
Figure 8 (draft, 170 x 150 mm): (a) excess E of the consensus pairs over matched random pairs per human dataset with the pooled estimates (diamonds) and the +/-0.10 margin; (b) pair-level pooled rho for the pairs
 evaluable in at least 3 human datasets (dots = datasets, diamond = pooled); (c) log2 fold changes of the 17 genes in GSE60803 (beta-catenin shRNA) against the INS-1 DKK-1 arm.
Output: figures/Figure_7_v2_islet_crosscheck.{png,pdf}, figures/Figure_8_draft_islet_meta.{png,pdf}.   Usage: python -I make_fig7_v2.py"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results'); OUT = os.path.join(HERE, 'figures'); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5})
GENES = ['Hes1', 'Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r']
NAVY, LBLUE, SALMON, SLATE, GREY, GREEN, RED = '#3C5488', '#4DBBD5', '#F39B7F', '#8491B4', '#9A9A9A', '#00A087', '#E64B35'
DS = [('gse81608', 'GSE81608\nhuman', NAVY), ('gse84133_human', 'GSE84133\nhuman', LBLUE), ('gse86469', 'GSE86469\nhuman', GREEN), ('gse81547', 'GSE81547\nhuman', RED), ('emtab5061', 'E-MTAB-5061\nhuman', SLATE), ('gse84133_mouse', 'GSE84133\nmouse', SALMON)]


def rd(tag, mode, part, variant='primary'): return pd.read_csv(os.path.join(RES, f'{tag}__{variant}__{mode}__{part}.csv'))


def panel_tag(fig, x, y, s): fig.text(x, y, s, fontsize=9, fontweight='bold', va='top')


# ------------------------------------------------------------------------------------------------------------------------------ Figure 7
summ = {t: json.load(open(os.path.join(RES, f'{t}__primary__summary.json')))['adjusted'] for t, _, _ in DS}
cons = rd('gse81608', 'unadjusted', 'pairs'); cons = cons[cons.consensus_invitro][['gene1', 'gene2']].values.tolist(); assert len(cons) == 21
fig = plt.figure(figsize=(6.69, 7.5))
top = fig.add_gridspec(1, 1, left=0.105, right=0.985, top=0.945, bottom=0.60)
bot = fig.add_gridspec(1, 6, left=0.185, right=0.985, top=0.50, bottom=0.055, wspace=0.12)
panel_tag(fig, 0.005, 0.985, '(a)'); panel_tag(fig, 0.005, 0.555, '(b)')
ax = fig.add_subplot(top[0, 0]); M = np.full((len(GENES), len(DS)), np.nan)
for j, (t, _, _) in enumerate(DS):
    d = rd(t, 'unadjusted', 'detection').set_index('gene_invitro').detected_pct; M[:, j] = d.reindex(GENES).values
for i in range(len(GENES)):
    for j in range(len(DS)):
        v = M[i, j]
        if np.isnan(v): ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor='white', edgecolor='#D0D0D0', lw=0.4)); ax.text(j, i, 'absent', fontsize=5.5, color=GREY, ha='center', va='center'); continue
        col = DS[j][2] if v >= 10 else '#E9E9E9'
        ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor=col, alpha=(0.25 + 0.75 * min(v, 100) / 100) if v >= 10 else 1.0, edgecolor='white', lw=0.6))
        ax.text(j, i, ('%.0f' % v) if v >= 1 else ('%.1f' % v), fontsize=6.2, ha='center', va='center', color='black' if (v < 55 or v < 10) else 'white')
ax.set_xlim(-0.5, len(DS) - 0.5); ax.set_ylim(len(GENES) - 0.5, -0.5); ax.set_yticks(range(len(GENES))); ax.set_yticklabels(GENES, fontstyle='italic'); ax.set_xticks(range(len(DS)))
ax.set_xticklabels(['%s\n(%d β-cells)' % (lab, summ[t]['beta_cells']) for t, lab, _ in DS], fontsize=6.8, fontweight='bold'); ax.xaxis.tick_top(); ax.tick_params(length=0)
for sp in ax.spines.values(): sp.set_visible(False)
for k, (t, lab, col) in enumerate(DS):
    ax = fig.add_subplot(bot[0, k]); c = rd(t, 'adjusted', 'cons'); got = {(r.gene1, r.gene2): r for r in c.itertuples()}; yp = np.arange(len(cons))[::-1]
    ax.axvline(0, color='black', lw=0.5, zorder=1)
    for (a, b), y in zip(cons, yp):
        if (a, b) in got:
            r = got[(a, b)]; ax.plot([r.null_mean, r.null_95pct], [y, y], color='#C9C9C9', lw=2.8, solid_capstyle='butt', zorder=2); ax.scatter([r.rho], [y], s=16, color=col, zorder=3, linewidths=0)
        else: ax.text(0.3, y, 'NE', fontsize=6, color='#B5B5B5', ha='center', va='center')
    ax.set_ylim(-0.8, len(cons) - 0.2); ax.set_xlim(-0.25, 0.72); ax.set_xticks([0, 0.3, 0.6]); ax.set_yticks(yp)
    ax.set_yticklabels([f'{a}–{b}' for a, b in cons] if k == 0 else [], fontstyle='italic', fontsize=6.8); ax.grid(axis='y', color='#EFEFEF', lw=0.5, zorder=0)
    ax.set_title('%s\n%d evaluable' % (lab, len(got)), fontsize=7, fontweight='bold', linespacing=1.15); ax.set_xlabel('Partial Spearman ρ', fontsize=7)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
fig.savefig(os.path.join(OUT, 'Figure_7_v2_islet_crosscheck.png'), dpi=300); fig.savefig(os.path.join(OUT, 'Figure_7_v2_islet_crosscheck.pdf')); plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------------------ Figure 8 (draft)
ms = json.load(open(os.path.join(HERE, 'meta', 'meta_summary.json'))); V = {(v['variant'], v['mode']): v for v in ms['variants']}
T = pd.read_csv(os.path.join(HERE, 'meta', 'study_effects__primary__adjusted.csv')); PP = pd.read_csv(os.path.join(HERE, 'meta', 'pair_level_pooled_rho.csv')); GP = pd.read_csv(os.path.join(HERE, 'perturbation', 'gse60803_fold_changes.csv'))
fig = plt.figure(figsize=(6.69, 6.0))
g1 = fig.add_gridspec(1, 1, left=0.30, right=0.78, top=0.945, bottom=0.67); g2 = fig.add_gridspec(1, 1, left=0.215, right=0.62, top=0.54, bottom=0.075); g3 = fig.add_gridspec(1, 1, left=0.745, right=0.985, top=0.54, bottom=0.075)
panel_tag(fig, 0.005, 0.985, '(a)'); panel_tag(fig, 0.005, 0.585, '(b)'); panel_tag(fig, 0.665, 0.585, '(c)')
ax = fig.add_subplot(g1[0, 0]); cols = {t: c for t, _, c in DS}
ax.axvspan(-0.10, 0.10, color='#EDEDED', zorder=0); ax.axvline(0, color='black', lw=0.5, zorder=1)
ylab, yy, txt = [], [], []                                  # row labels, row positions, numeric column on the right
for y, r in zip([8, 7, 6, 5, 4], T.itertuples()):
    ax.plot([r.E - 1.96 * r.SE, r.E + 1.96 * r.SE], [y, y], color=cols[r.tag], lw=1.1, zorder=2); ax.scatter([r.E], [y], marker='s', s=16 + 22 * np.log10(r.beta_cells) / 3.4, color=cols[r.tag], zorder=3, linewidths=0)
    ylab.append('%s (%d β-cells, %d pairs)' % (r.label.split(' (')[0], r.beta_cells, r.pairs)); yy.append(y); txt.append('%+.3f [%+.3f, %+.3f]' % (r.E, r.E - 1.96 * r.SE, r.E + 1.96 * r.SE))
def pooled(y, rec, name, col):
    r = [e for e in rec['estimators'] if e['estimator'].startswith('REML+HK (primary)')][0]
    ax.fill([r['lo'], r['mu'], r['hi'], r['mu']], [y, y + 0.28, y, y - 0.28], color=col, zorder=3, linewidth=0); ylab.append(name); yy.append(y); txt.append('%+.3f [%+.3f, %+.3f]' % (r['mu'], r['lo'], r['hi']))
pooled(2.4, V[('primary', 'adjusted')], 'Pooled, depth-adjusted (primary)', 'black'); pooled(1.4, V[('nondiabetic', 'adjusted')], 'Pooled, non-diabetic donors only', '#6B6B6B'); pooled(0.4, V[('primary', 'unadjusted')], 'Pooled, not depth-adjusted', '#6B6B6B')
ax.axhline(3.2, color='#BDBDBD', lw=0.5, zorder=1)
for y, t_ in zip(yy, txt): ax.text(1.03, y, t_, transform=ax.get_yaxis_transform(), fontsize=6.3, va='center', ha='left')
ax.set_yticks(yy); ax.set_yticklabels(ylab, fontsize=6.5); ax.set_ylim(-0.5, 8.7); ax.set_xlim(-0.22, 0.42); ax.set_xlabel('Excess co-expression of the consensus pairs over matched\nrandom pairs, E (Fisher z)', fontsize=7)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
ax = fig.add_subplot(g2[0, 0]); P = PP.sort_values('pooled_rho').reset_index(drop=True); yp = np.arange(len(P)); ax.axvline(0, color='black', lw=0.5, zorder=1)
for y, r in zip(yp, P.itertuples()):
    ax.plot([r.lo, r.hi], [y, y], color='#8C8C8C', lw=1.1, zorder=2)
    for t, v in zip(r.tags.split(';'), r.rho_by_dataset.split(';')): ax.scatter([float(v)], [y], s=9, color=cols[t], zorder=3, linewidths=0)
    ax.fill([r.pooled_rho - 0.012, r.pooled_rho, r.pooled_rho + 0.012, r.pooled_rho], [y, y + 0.3, y, y - 0.3], color='black', zorder=4, linewidth=0)
ax.set_yticks(yp); ax.set_yticklabels(['%s–%s (%d)' % (r.gene1, r.gene2, r.datasets) for r in P.itertuples()], fontstyle='italic', fontsize=6.5); ax.set_xlabel('Pooled Spearman ρ (95% CI)', fontsize=7); ax.set_ylim(-0.7, len(P) - 0.3)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
ax = fig.add_subplot(g3[0, 0]); D = GP.dropna(subset=['log2fc_all_donors']); ax.axhline(0, color='black', lw=0.5, zorder=1); ax.axvline(0, color='black', lw=0.5, zorder=1)
ax.scatter(D.ins1_log2fc_DKK1, D.log2fc_all_donors, s=14, color=NAVY, zorder=3, linewidths=0)
for r in D.itertuples():
    if abs(r.ins1_log2fc_DKK1) > 0.4 or abs(r.log2fc_all_donors) > 0.2: ax.annotate(r.gene, (r.ins1_log2fc_DKK1, r.log2fc_all_donors), fontsize=5.5, fontstyle='italic', xytext=(2, 2), textcoords='offset points')
ax.set_xlabel('INS-1, DKK-1 vs control\n(log2 fold change)', fontsize=7); ax.set_ylabel('GSE60803, β-catenin shRNA\nvs control (log2 fold change)', fontsize=7)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
fig.savefig(os.path.join(OUT, 'Figure_8_draft_islet_meta.png'), dpi=300); fig.savefig(os.path.join(OUT, 'Figure_8_draft_islet_meta.pdf')); plt.close(fig)
print('saved', OUT)
