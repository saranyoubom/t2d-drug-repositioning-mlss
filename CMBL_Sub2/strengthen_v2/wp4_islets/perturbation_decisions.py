"""perturbation_decisions.py - protocol 4.6, step 2: the keep / drop decision with its reason for EVERY candidate of perturbation/candidates.csv, made from title, organism, assay type and summary text only (no expression
value was read). Rules (PRESPEC_WP4_choices.md section D): keep = expression profiling of islets, beta-cell lines or beta-cell-derived cells with an INHIBITION of the Notch (gamma-secretase) or Wnt (DKK1, Wnt inhibitor,
beta-catenin or TCF7L2 loss) pathway against a control, at least 3 samples per arm, on a platform that covers the 17 genes. Everything else is dropped with the reason below. Series read in full before deciding are marked 'read'.
Output: perturbation/candidates_decisions.csv.   Usage: python -I perturbation_decisions.py"""
import os, sys
sys.stdout.reconfigure(encoding='utf8')
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); D = pd.read_csv(os.path.join(HERE, 'perturbation', 'candidates.csv'))
MANUAL = {
    'GSE60803': ('keep', 'read: adult human beta-cell-derived cells, beta-catenin shRNA against control (WNT inhibition), several donors; expanded adult cells, not a stem-cell differentiation design (the regex flag fires on the word redifferentiation); sample arms to be checked in the series matrix before analysis'),
    'GSE44045': ('drop', 'CORRECTED 6 Nov (addendum H): the sample records show that the arrays compare islets of rats fed a standard diet with a cafeteria diet (the series title refers to Sfrp5); no Notch or Wnt perturbation (K2). The first version of this reason described the series as Sfrp5 silencing (Wnt activation) from the title and summary; that description is withdrawn'),
    'GSE44046': ('drop', 'CORRECTED 6 Nov (addendum H): same study as GSE44045 (30-day arm), a diet comparison in the sample records; no Notch or Wnt perturbation (K2)'),
    'GSE30526': ('drop', 'read: COUP-TFII loss or gain of function; beta-catenin is a downstream target, the perturbation is of a transcription factor, not of the Wnt pathway'),
    'GSE12802': ('drop', 'read: small-molecule screen for beta-cell proliferation in an immortalised mouse line; compounds ACTIVATE Wnt (thiophene-pyrimidines) or L-type calcium channels; no pathway inhibition arm'),
    'GSE193888': ('drop', 'read: doxycycline-inducible beta-cell Notch GAIN of function in mouse islets; direction opposite to DAPT'),
    'GSE249170': ('drop', 'read: beta-cell Jagged1 gain and loss of function in mice (ligand level, in vivo); not a gamma-secretase or pathway inhibition'),
    'GSE139147': ('drop', 'read: HDAC inhibitor MS-275 in BRIN-BD11 cells; not a Notch or Wnt perturbation'),
    'GSE224490': ('drop', 'read: BCAS2 splicing factor in mouse islets; not a Notch or Wnt perturbation'),
    'GSE136887': ('drop', 'read: single-cell gain-of-function screen of T2D regulators in human islets; not a Notch or Wnt inhibition'),
    'GSE195647': ('drop', 'read: mammary stem cell ageing; not islets'),
    'GSE41284': ('drop', 'title: Tcf7l2 in liver (perinatal hepatic response); not islets or beta-cells'),
    'GSE53352': ('drop', 'title: haematopoietic stem cells deficient of glycogen synthase kinase; not islets'),
    'GSE50103': ('drop', 'title: bi-potent pancreas progenitor expansion (Lgr5 / R-spondin); progenitor design'),
    'GSE14': ('drop', 'title: Cancer Genome Anatomy Project SAGE libraries; not a perturbation design'),
    'GSE118725': ('drop', 'title: transcription-factor binding to non-coding variants (reporter and binding assays); not an expression perturbation'),
    'GSE126550': ('drop', 'title: saturation mutagenesis of regulatory elements; reporter assay'),
    'GSE308943': ('drop', 'title: microRNA profiling of rat islet cells; non-coding RNA assay, no Notch or Wnt perturbation'),
    'GSE323218': ('drop', 'title: MPRA control from MIN6; reporter assay'),
    'GSE323299': ('drop', 'title: MPRA from MIN6; reporter assay'),
    'GSE228665': ('drop', 'title: isogenic hESC knock-outs of T2D genes and their differentiation; stem-cell design'),
}
rows = []
for _, r in D.iterrows():
    if r.accession in MANUAL: dec, why = MANUAL[r.accession]
    elif 'Expression profiling' not in str(r.assay): dec, why = 'drop', 'not an expression profile (assay: %s)' % r.assay
    elif r.stem_cell_flag: dec, why = 'drop', 'title / summary describe a stem-cell, progenitor, reprogramming or lineage-differentiation design (flag from the regex, checked against the title)'
    else: dec, why = 'drop', 'no Notch or Wnt inhibition arm in the title or summary'
    rows.append(dict(accession=r.accession, queries=r.queries, organism=r.organism, assay=r.assay, n_samples=r.n_samples, platform=r.platform, decision=dec, reason=why, title=r.title))
O = pd.DataFrame(rows); O.to_csv(os.path.join(HERE, 'perturbation', 'candidates_decisions.csv'), index=False)
print(O.decision.value_counts().to_dict()); print(O[O.decision == 'keep'][['accession', 'organism', 'n_samples', 'platform', 'title']].to_string(index=False))
print(O[O.reason.str.startswith('title / summary describe')].accession.tolist())
