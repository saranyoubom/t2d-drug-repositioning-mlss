"""make_fig10.py - Fig. 10 of the JTM manuscript (20 Nov 2026, main session): the candidate checks for pooled co-expression networks from small experiments, with the evidence behind each step.
A flow of six steps; for each: the check, what the benchmark and the INS-1 data show (every number carries its status: simulation, benchmark, post hoc, INS-1) and how to read it.
Plots no new analysis. Every number is read from a saved result file of the study (calibration_tables.csv, calibration_sensitivity_no_degenerate.csv, design_table.csv, validation_endpoints_long.csv,
post_hoc/results_*/summary.json) or is an INS-1 number of the manuscript text; fig10_numbers.csv lists each number with its source, and the script checks that the INS-1 numbers occur in the manuscript text.
Style of Figs 7-9 (figstyle.py of wp6a_assets, imported read-only; output goes to THIS folder, never to wp6a_assets): 170 mm wide, Arial >= 6 pt, PNG 300 dpi and PDF, no key strip,
statistical parameters (P, q, r, D, T) in italics.
Usage (python -I): python -I make_fig10.py validation      final figure (needs post_hoc/results_validation)
                   python -I make_fig10.py dev              layout test with development-series numbers; the file name carries _TEST and must not be used"""
import os, sys, json, csv, re, zipfile, html
HERE = os.path.dirname(os.path.abspath(__file__)); os.environ['WP6A_OUT'] = HERE
V2 = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(V2, 'wp6a_assets'))
import numpy as np, pandas as pd
import figstyle as fs
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle

SOURCE = sys.argv[1] if len(sys.argv) > 1 else 'validation'
assert SOURCE in ('validation', 'dev')
fs.setup()
W2 = os.path.join(V2, 'wp2_benchmark2'); PH = os.path.join(V2, 'ars_review_191126', 'post_hoc')
MS_PATH = os.environ.get('FIG10_MANUSCRIPT', os.path.abspath(os.path.join(V2, '..', '..', 'Manuscript', 'Manuscript_JTM_v2_SUBMISSION_191126c.docx')))

# ------------------------------------------------------------------ numbers and their sources
trace = []
def put(key, value, source): trace.append((key, value, source)); return value
cal = pd.read_csv(os.path.join(W2, 'calibration_tables.csv'))
def cal_row(r, m, stat): return cal[(cal.r == r) & (cal.m == m) & (cal.statistic == stat)].iloc[0]
t1 = cal_row(0.0, 0, 'T_r2'); p_tr2 = cal_row(0.7, 6, 'T_r2'); c07 = cal_row(0.7, 6, 'count'); c09 = cal_row(0.9, 6, 'count'); c05 = cal_row(0.5, 6, 'count')
pct = lambda x: '%d%%' % round(100 * x)
N = {}
N['typeI'] = put('T_r2 type I error', '%.1f%%' % (100 * t1.rate), 'calibration_tables.csv T_r2 r=0 m=0: %d of %d' % (t1.n_reject, t1.n_tests))
N['typeI_k'] = (int(t1.n_reject), int(t1.n_tests))
nd = pd.read_csv(os.path.join(W2, 'calibration_sensitivity_no_degenerate.csv')).set_index('statistic').loc['T_r2']
N['typeI_nd'] = put('T_r2 type I error without 2 series with constant genes (post hoc)', '%.1f%%' % (100 * nd.type1), 'calibration_sensitivity_no_degenerate.csv: %d of %d' % (nd.type1_k, nd.type1_n))
N['power_tr2'] = put('T_r2 power r 0.7, 6 genes', pct(p_tr2.rate), 'calibration_tables.csv: %d of %d' % (p_tr2.n_reject, p_tr2.n_tests))
N['count07'] = put('count test power r 0.7', pct(c07.rate), 'calibration_tables.csv: %d of %d' % (c07.n_reject, c07.n_tests))
N['count09'] = put('count test power r 0.9', pct(c09.rate), 'calibration_tables.csv: %d of %d' % (c09.n_reject, c09.n_tests))
N['count05'] = put('count test power r 0.5', pct(c05.rate), 'calibration_tables.csv: %d of %d' % (c05.n_reject, c05.n_tests))
dt = pd.read_csv(os.path.join(W2, 'design_table.csv')); d4 = dt[(dt.n_per_group == 4) & (dt.r == 0.7)].iloc[0]
N['d4_power'] = put('design 4 per group: T_r2 power r 0.7', pct(d4.power_T_r2), 'design_table.csv n=4 r=0.7')
N['d4_pairs'] = put('design 4 per group: share of true module pairs detected', pct(d4.pair_power), 'design_table.csv n=4 r=0.7')
N['d4_rbest'] = put('design 4 per group: |r| the best single pair needs', '%.2f' % d4.r_best_pair, 'design_table.csv n=4')
d6 = dt[(dt.n_per_group == 6) & (dt.r == 0.5)].iloc[0]; d8 = dt[(dt.n_per_group == 8) & (dt.r == 0.5)].iloc[0]
N['d6_8'] = put('design 6 and 8 per group: T_r2 power r 0.5', '%.2f to %.2f' % (d6.power_T_r2, d8.power_T_r2), 'design_table.csv n=6, n=8 r=0.5')
vl = pd.read_csv(os.path.join(W2, 'validation_endpoints_long.csv'))
va = vl[(vl.arm.str.startswith('All validation sets pooled')) & (vl.endpoint == 'coupling_T_r2') & (vl.gene_set == 'A')].iloc[0]
N['flag_k'] = (int(va.k), int(va.n)); N['flag_pct'] = put('T_r2 flagged, validation sets, gene set A', pct(va.estimate), 'validation_endpoints_long.csv: %d of %d' % (va.k, va.n))
ph_path = os.path.join(PH, 'results_validation' if SOURCE == 'validation' else 'results_dev', 'summary.json')
S = json.load(open(ph_path, encoding='utf-8')); SA, SB, SAB = S['sets']['A'], S['sets']['B'], S['sets']['A+B']
q2 = lambda v: '%.2f' % v
N['A_obs'] = put('share carried by group means, gene set A (median of series)', q2(SA['obs_share_quartiles'][1]), os.path.relpath(ph_path, V2) + ' sets.A.obs_share_quartiles')
N['A_exp'] = put('expected share without coupling, gene set A', q2(SA['null_share_quartiles'][1]), os.path.relpath(ph_path, V2) + ' sets.A.null_share_quartiles')
N['B_obs'] = put('share carried by group means, gene set B', q2(SB['obs_share_quartiles'][1]), os.path.relpath(ph_path, V2) + ' sets.B.obs_share_quartiles')
N['B_exp'] = put('expected share without coupling, gene set B', q2(SB['null_share_quartiles'][1]), os.path.relpath(ph_path, V2) + ' sets.B.null_share_quartiles')
N['above'] = put('matrices above the 97.5th null percentile', '%d of %d' % (SAB['n_above_null_q975'], SAB['n_share_defined']), os.path.relpath(ph_path, V2) + ' sets.A+B')
N['pc1_k'] = put('T_r2 flagged after removal of PC1, gene set A', '%d%%' % round(100 * SA['flagged_pc1'] / SA['n_matrices']), os.path.relpath(ph_path, V2) + ' sets.A: %d of %d' % (SA['flagged_pc1'], SA['n_matrices']))
N['orig_k'] = put('T_r2 flagged, original, same sets (gene set A)', '%d%%' % round(100 * SA['flagged_original'] / SA['n_matrices']), os.path.relpath(ph_path, V2) + ' sets.A: %d of %d' % (SA['flagged_original'], SA['n_matrices']))
N['n_sets'] = SA['n_matrices']
if SOURCE == 'validation':                                   # the original flagged fraction recomputed by the post hoc script must equal the benchmark's (gate)
    assert N['orig_k'] == N['flag_pct'] and SA['flagged_original'] == N['flag_k'][0] and SA['n_matrices'] == N['flag_k'][1] == 59, 'post hoc validation results do not reproduce the benchmark prevalence'
I = json.load(open(os.path.join(PH, 'results_ins1', 'summary.json'), encoding='utf-8'))
N['I_obs'] = put('INS-1 share carried by group means', q2(I['share']['obs_median']), 'results_ins1/summary.json'); N['I_exp'] = put('INS-1 expected share', q2(I['share']['null_mean']), 'results_ins1/summary.json')
# INS-1 numbers of the manuscript text (checked against the text below)
MS_NUM = {'null mean 22.1': '22.1', 'P = 0.68': '0.68', 'Tr2 P = 0.061': '0.061', 'same-run r': '+0.149', 'between-run r': '−0.019', 'D': '+0.168', 'D P': '0.020', 'row-layout simulation': '+0.520', 'observed same-run, 26 August': '+0.536'}
for k, v in MS_NUM.items(): trace.append(('INS-1 ' + k, v, 'manuscript Results, paragraphs on the permutation null, the coupling test and the run structure'))


def manuscript_text(path):
    z = zipfile.ZipFile(path); x = z.read('word/document.xml').decode('utf8'); body = x[x.index('<w:body>'):]
    return '\n'.join(''.join(html.unescape(t) for t in re.findall(r'<w:t(?: [^>]*)?>([^<]*)</w:t>', p)) for p in re.findall(r'<w:p[ >].*?</w:p>', body, flags=re.S))
MST = manuscript_text(MS_PATH)
missing = [v for v in MS_NUM.values() if v not in MST]
assert not missing, 'numbers not found in the manuscript text: %s' % missing

# ------------------------------------------------------------------ text of the figure
it = fs.it
TR2 = r'$\mathit{T}_{\mathit{r}\mathrm{2}}$'
rr, PP, QQ, DD = it('r'), it('P'), it('q'), it('D')
ROWS = [
 ('Design and genes',
  'With 4 replicates per group, %s detects one module of 6 genes at %s = 0.7 with power %s, but only %s of the module’s pairs reach %s < 0.05, and the best single pair would need |%s| ≥ %s. Six to eight independent experiments per group give power %s at %s = 0.5 (simulation).' % (TR2, rr, N['d4_power'], N['d4_pairs'], QQ, rr, N['d4_rbest'], N['d6_8'], rr),
  'Wells from one preparation are not independent experiments. The design limits what any pair-level result can show.'),
 ('Permutation null of the consensus count',
  'The count detects injected within-group coupling in %s of tests (%s = 0.7, 6 genes; %s at %s = 0.9, %s at %s = 0.5) (benchmark, simulation). INS-1: 21 pairs, null mean 22.1, %s = 0.68.' % (N['count07'], rr, N['count09'], rr, N['count05'], rr, PP),
  'A count inside its null does not show that coupling is absent.'),
 ('Share carried by group means, against its expectation',
  'Expected without within-group coupling (post hoc, %d public data sets): genes selected for a group effect %s observed, %s expected; random expressed genes %s observed, %s expected. INS-1: %s observed, %s expected. The observed share was above its null in %s matrices.' % (N['n_sets'], N['A_obs'], N['A_exp'], N['B_obs'], N['B_exp'], N['I_obs'], N['I_exp'], N['above']),
  'A share near 1 reflects the selected genes and the three groups of four, not co-regulation.'),
 ('Within-group coupling test, %s' % TR2,
  'Power %s (%s = 0.7, 6 genes). Type I error %s, above the limit of 6%% set in the protocol (%s after a post hoc exclusion of 2 series with constant genes). Flagged in %s of 59 validation sets (gene set A); %s after removal of one shared axis (post hoc). INS-1: %s = 0.061.' % (N['power_tr2'], rr, N['typeI'], N['typeI_nd'], N['flag_pct'], N['pc1_k'], PP),
  'A flag shows dependence within groups, not its source: shared sample-level variation or shared biology.'),
 ('Treatment adjustment and ratio schemes',
  'After adjustment for treatment no pair was significant in all four schemes, in INS-1 or (median) in the public sets. Counts of the ratio schemes sat at their own null (median observed to expected 1.05 to 1.19).',
  'Ratio schemes share a denominator: their counts are not findings.'),
 ('Run and reference structure',
  'INS-1: within-group %s was +0.149 among same-run pairs and −0.019 between runs (%s = +0.168, %s = 0.020); the row layout of one run alone gave +0.520 (observed +0.536).' % (rr, DD, PP),
  'Compatible with a shared reference, plate or layout, not with co-regulation.'),
]
NOTE = 'Candidate checks. Whether they improve the conclusions drawn from small experiments, and whether they apply to toxicogenomic data, was not tested.'

# ------------------------------------------------------------------ layout in millimeters
FS = 6.5; LH = 2.95; PAD = 1.6; GAP = 2.0; FIGW = fs.FIGW
CX, CW = 9.5, 38.0; EX, EW = 49.5, 72.0; RX, RW = 123.5, 46.0
probe = fs.new_fig(50.0); rend = probe.canvas.get_renderer(); dpi = probe.dpi


def wmm(s, weight='normal', size=FS):
    t = probe.text(0, 0, s, fontsize=size, fontweight=weight); w = t.get_window_extent(rend).width / dpi * 25.4; t.remove(); return w


def nb(s):
    """non-breaking spaces inside 'x = value' and 'x < value' so that a statistic and its value are never split across two lines"""
    return re.sub(r'(\S) ([=<≥]) (?=[−+\-]?\d)', lambda m: m.group(1) + ' ' + m.group(2) + ' ', s)


def wrap(s, width, weight='normal'):
    lines, cur = [], ''
    for w in nb(s).split(' '):
        trial = (cur + ' ' + w) if cur else w
        if wmm(trial, weight) <= width - 2 * 1.4 or not cur: cur = trial
        else: lines.append(cur); cur = w
    lines.append(cur); return lines


rows = []; y = 9.0
for i, (title, ev, rd) in enumerate(ROWS, 1):
    tl, el, rl = wrap(title, CW, 'bold'), wrap(ev, EW), wrap(rd, RW); n = max(len(tl), len(el), len(rl)); h = n * LH + 2 * PAD
    rows.append((i, y, h, tl, el, rl)); y += h + GAP
nl = wrap(NOTE, FIGW - 2 * 1.0); nh = len(nl) * LH + 2 * PAD; H = y + 1.0 + nh + 1.5
plt.close(probe)
fig = fs.new_fig(H); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, FIGW); ax.set_ylim(H, 0); ax.axis('off')
for x, w, lab in ((CX, CW, 'Check'), (EX, EW, 'What the benchmark and the INS-1 data show'), (RX, RW, 'How to read it')): ax.text(x + 1.4, 4.8, lab, fontsize=7, fontweight='bold', va='center', ha='left')
prev = None
for i, y0, h, tl, el, rl in rows:
    for x, w, lines, fill, bold in ((CX, CW, tl, '#e8f1fc', 'bold'), (EX, EW, el, 'white', 'normal'), (RX, RW, rl, '#f4f3f0', 'normal')):
        ax.add_patch(FancyBboxPatch((x, y0), w, h, boxstyle='round,pad=0,rounding_size=1.2', fc=fill, ec=fs.GREYLINE, lw=0.5))
        for k, ln in enumerate(lines): ax.text(x + 1.4, y0 + PAD + (k + 0.5) * LH, ln, fontsize=FS, fontweight=bold, va='center', ha='left')
    cy = y0 + h / 2; ax.add_patch(Circle((4.6, cy), 2.6, fc=fs.BLUE, ec='none')); ax.text(4.6, cy + 0.1, str(i), fontsize=7.5, fontweight='bold', color='white', ha='center', va='center')
    if prev is not None: ax.plot([4.6, 4.6], [prev + 2.6, cy - 2.6], color=fs.GREYLINE, lw=0.8, solid_capstyle='butt')
    prev = cy
ny = H - 1.0 - nh
ax.add_patch(FancyBboxPatch((1.0, ny), FIGW - 2.0, nh, boxstyle='round,pad=0,rounding_size=1.2', fc='#f4f3f0', ec=fs.GREYLINE, lw=0.5))
for k, ln in enumerate(nl): ax.text(2.4, ny + PAD + (k + 0.5) * LH, ln, fontsize=FS, va='center', ha='left')
stem = 'Figure_10_candidate_checks' + ('' if SOURCE == 'validation' else '_TEST')
fs.save(fig, stem)
with open(os.path.join(HERE, 'fig10_numbers%s.csv' % ('' if SOURCE == 'validation' else '_TEST')), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f); w.writerow(['number', 'value', 'source']); w.writerows(trace)
print('figure height %.1f mm, %d rows; numbers traced: %d' % (H, len(rows), len(trace)))
