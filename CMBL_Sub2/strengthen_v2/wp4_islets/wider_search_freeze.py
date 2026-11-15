"""wider_search_freeze.py - B2/B3 of PRESPEC_WP4b_positive_control_and_wider_search.md: the keep / drop decision with its reason for EVERY one of the candidates of wider_search/candidates_all.csv, made from title,
summary, design text and sample records only (no expression value was read), and the contrast of every kept series. Generic rules (organism, assay, sample number, stem-cell flag) give the mechanical reasons; the
series read in full carry a specific reason. The file is fingerprinted (sha256 in ../PROTOCOL_v2_SHA256.txt) BEFORE any expression value of a candidate is read.
Output: wider_search/frozen_candidates.csv.   Usage: python -I wider_search_freeze.py"""
import os, sys, re
sys.stdout.reconfigure(encoding='utf8')
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'wider_search')
D = pd.read_csv(os.path.join(OUT, 'candidates_all.csv'))
STEM = re.compile(r'differentiat|iPSC|\bESC\b|embryonic stem|pluripoten|stem cell|stem-cell|progenitor|reprogram', re.I)
# kept series: expected sign (+1 inhibition, -1 mirror), arm, contrast (patterns on sample titles), pairing, data route, K4 status
KEEP = {
    'GSE60803': dict(sign=+1, mirror=0, arm='DKK-1', weak='confounded: controls are not expanded cells (addendum G)', contrast='as fixed in addendum G of PRESPEC_WP4_choices.md: per donor, control = untreated or nontarget array, perturbed = mean of the beta-catenin shRNA arrays',
                     pairing='by donor', route='series matrix (GPL11532); already analysed', k4='obtainable',
                     why='read: adult human beta-cell-derived cells, beta-catenin shRNA against control (Wnt inhibition), 4 donors (K1, K2, K3 met)'),
    'GSE167073': dict(sign=+1, mirror=0, arm='DKK-1', weak='receptor-level loss (Fzd9) in isolated islets of an insulinoma model under MycER activation', contrast='perturbed = samples titled "Fzd9 ko" (4); control = samples titled "Fzd9 wt" (4)',
                      pairing='none (unpaired)', route='series matrix values table (41346 rows, GPL17400); platform annotation from the GEO platform text export', k4='obtainable (annotation to be fetched at the analysis step)',
                      why='read: Fzd9 knock-out against wild type in islets of MycER-BclXL mice, 4 vs 4 (K1 islets of an insulinoma model, K2 Wnt receptor loss, K3 met)'),
    'GSE99503': dict(sign=-1, mirror=1, arm='DKK-1', weak='GSK3 inhibitor 1-azakenpaullone is not specific to the Wnt pathway', contrast='perturbed = "Neg siRNA + 1-AKP #1-3" (scramble siRNA plus GSK3 inhibitor); control = "Neg siRNA #1-3" (scramble siRNA, DMSO); the TUNAR siRNA arms are not used',
                     pairing='none (passages are not stated as pairs)', route='NCBI-generated raw counts (GRCh38.p13, Entrez ids; HEAD 200, 566 KB)', k4='obtainable',
                     why='read: EndoC-betaH1 cells, GSK3 inhibition (Wnt activation, a MIRROR design) against vehicle in scramble-siRNA cells, 3 vs 3 (K1, K2, K3 met)'),
    'GSE193888': dict(sign=-1, mirror=1, arm='DAPT', weak='', contrast='perturbed = "NICD-1..4" (tetO_NICD+, doxycycline 8 weeks); control = "Control-1..4" (tetO_NICD-)',
                      pairing='none (unpaired)', route='none by plain download: no supplementary file list (HTTP 404) and no NCBI-generated counts for mouse (HTTP 404)', k4='NOT obtainable',
                      why='read: beta-cell Notch gain of function in mouse islets (a MIRROR design for the DAPT arm), 4 vs 4 (K1, K2, K3 met); data not obtainable'),
    'GSE249170': dict(sign=-1, mirror=1, arm='DAPT', weak='', contrast='perturbed = "Jag1-1..3" (tetO_Jag1+, doxycycline 8 weeks); control = "Control-1..3" (tetO_Jag1-)',
                      pairing='none (unpaired)', route='none by plain download: no supplementary file list (HTTP 404) and no NCBI-generated counts for mouse (HTTP 404)', k4='NOT obtainable',
                      why='read: beta-cell Jagged1 gain of function in mouse islets (a MIRROR design for the DAPT arm), 3 vs 3 (K1, K2, K3 met); data not obtainable'),
    'E-MTAB-2034': dict(sign=+1, mirror=0, arm='DKK-1', weak='', contrast='perturbed = "TCF7L2 knockdown_Replicate 1-4"; control = "Control sample-Replicate 1-4"; paired by replicate (passage) as named in the sample records',
                        pairing='by replicate (passages 62, 63, 81, 82)', route='none: the ArrayExpress record holds raw reads only (ENA); no processed expression table; reads are not aligned here', k4='NOT obtainable',
                        why='read: TCF7L2 knockdown against control in rat INS-1 832/13 cells, 4 vs 4 (K1, K2, K3 met; the closest model to the INS-1 experiment); no processed data'),
}
DROP = {   # specific reasons for candidates read in full or from the sample records
    'GSE44045': 'read, sample records: the arrays compare islets of rats fed a standard diet with a cafeteria diet (the series title refers to Sfrp5, the array experiment is a diet comparison); no Notch or Wnt perturbation (K2)',
    'GSE44046': 'read, sample records: same study as GSE44045 (30-day arm), a diet comparison; no Notch or Wnt perturbation (K2)',
    'GSE44047': 'read: SuperSeries of GSE44045 and GSE44046 (diet comparison); no Notch or Wnt perturbation (K2)',
    'GSE12802': 'read, sample records: compounds (BayK8644, exendin-4, their combination) against DMSO over a time course; no Wnt or Notch arm (K2)',
    'GSE30526': 'read: COUP-TFII loss or gain of function; beta-catenin is a downstream target, the perturbation is of a transcription factor (K2)',
    'GSE64444': 'read: MIN6 culture conditions that change Hes3; Hes3 is not a Notch pathway perturbation (K2)', 'GSE64449': 'read: Hes3 knockdown in MIN6 cells; Hes3 is not a Notch pathway perturbation (K2)',
    'GSE7430': 'read: conditional beta-catenin knockout in the developing exocrine pancreas; embryonic pancreas tissue, not islets or beta-cells (K1)', 'GSE76189': 'read: Tspan2 and high glucose in beta-cells; no Notch or Wnt perturbation (K2)',
    'GSE186944': 'read: Cdh13 null islets; no Notch or Wnt perturbation (K2)', 'E-MEXP-3836': 'read: insulin receptor overexpression in INS-1 cells; Wnt activation is a consequence, not the perturbation (K2)',
    'E-MTAB-5325': 'title: pericytes of Tcf7l2-mutated mice; not islets or beta-cells (K1)', 'GSE24296': 'title: anti-DLL4 treatment, immune cell populations; not islets or beta-cells (K1)',
    'GSE41284': 'title: Tcf7l2 in liver (K1)', 'GSE193888_dummy': '',
    'GSE139147': 'read: HDAC inhibitor in BRIN-BD11 cells; no Notch or Wnt perturbation (K2)', 'GSE224490': 'read: BCAS2 splicing factor; no Notch or Wnt perturbation (K2)',
    'GSE195647': 'read: mammary stem cell ageing (K1)', 'GSE136887': 'read: single-cell gain-of-function screen of T2D regulators in islets; no Notch or Wnt perturbation (K2)',
    'GSE125151': 'title: serpin B13 in islets; no Notch or Wnt perturbation (K2)', 'GSE134272': 'title: glucolipotoxicity and a protective compound in rat islets; no Notch or Wnt perturbation (K2)',
    'GSE17647': 'title: 4E-BP1 and HDL in beta-cells (K2)', 'GSE19877': 'title: paternal obesity and offspring islets (K2)', 'GSE44373': 'title: obese and lean rat islets, diet (K2)',
    'GSE142572': 'title: hPSCs treated with CHIR99021; stem cells, not beta-cells (K1)',
}
rows = []
for r in D.itertuples():
    acc = r.accession; n = pd.to_numeric(r.n_samples, errors='coerce'); org = str(r.organism); assay = str(r.assay); ae = str(r.entrytype) == 'ArrayExpress'
    stem = bool(STEM.search(str(r.title) + ' ' + str(r.summary)))
    base = dict(accession=acc, sources=r.sources, pmids=r.pmids, organism=org, assay=assay, n_samples=r.n_samples, title=r.title)
    if acc in KEEP:
        k = KEEP[acc]; rows.append(dict(base, decision='keep', reason=k['why'], expected_sign=k['sign'], mirror=k['mirror'], arm=k['arm'], weak_match_note=k['weak'], contrast=k['contrast'], pairing=k['pairing'], data_route=k['route'], k4=k['k4']))
        continue
    if acc in DROP and DROP[acc]: dec, why = 'drop', DROP[acc]
    elif not re.search(r'Homo sapiens|Mus musculus|Rattus norvegicus', org): dec, why = 'drop', 'organism is not human, mouse or rat (K1)'
    elif not (ae or 'Expression profiling' in assay): dec, why = 'drop', 'not an expression profile (assay: %s) (K1)' % assay
    elif not ae and pd.notna(n) and n < 6: dec, why = 'drop', 'fewer than 6 samples, so fewer than 3 per arm (K3)'
    elif stem: dec, why = 'drop', 'stem-cell, progenitor or reprogramming design (K1; flag from title or summary, checked against the title)'
    else: dec, why = 'drop', 'no Notch or Wnt perturbation of islets or beta-cells in the title or summary (K2)'
    rows.append(dict(base, decision=dec, reason=why, expected_sign='', mirror='', arm='', weak_match_note='', contrast='', pairing='', data_route='', k4=''))
F = pd.DataFrame(rows); F.to_csv(os.path.join(OUT, 'frozen_candidates.csv'), index=False)
print('candidates', len(F), '| keep', int((F.decision == 'keep').sum()), '| drop', int((F.decision == 'drop').sum()))
print(F[F.decision == 'keep'][['accession', 'organism', 'n_samples', 'expected_sign', 'mirror', 'arm', 'k4']].to_string(index=False))
print(F[F.decision == 'drop'].reason.str.replace(r'\(.*?\)$', '', regex=True).str.slice(0, 60).value_counts().head(14).to_string())
