"""
Audit of the raw rat C-peptide ELISA workbook (2020-06-08_ELISA new.xlsx, sheet 2020-06-05_UG_Rat_C-peptide_2) behind Fig. 3, Fig. 1c (DNA) and
Supplementary Table S7. Read-only. Writes tidy CSV files and a text report.
What it establishes (checked against the workbook's own formulas and against the data files used for the figures):
  * layout of the file: 3 groups x 5 glucose conditions; one DNA value per group x condition (15) and four rows (wells) per condition; the authors state four culture wells per condition (16 Oct 2026),
    which the file itself cannot show (the rows are labelled by condition only). An earlier reading as 'one culture well read four times' was an inference and was wrong.
  * the "ng/ug DNA" columns divide by the GROUP-MEAN DNA of the five wells of that group (cells L20, L25, L30), not by the well's own DNA
  * the four DNA values per group in Fig. 1c are the DNA of the 2.8, 5.5, 22 and 44 mM wells (the basal 'Normal KR' well is not used)
Usage: python -I elisa_raw_audit.py <elisa.xlsx> <fig1c.xlsx> <cpeptide_all_conditions.csv> <out_dir>
"""
import sys, os, csv, statistics as st
import openpyxl
sys.stdout.reconfigure(encoding="utf-8")
XL, F1C, CSVIN, OUT = sys.argv[1:5]
os.makedirs(OUT, exist_ok=True)

wb = openpyxl.load_workbook(XL, data_only=True); ws = wb.worksheets[0]
wf = openpyxl.load_workbook(XL, data_only=False).worksheets[0]
V = lambda a: ws[a].value
GROUPS = {"CTRL": dict(rows=range(20, 40), dna_rows=range(20, 25), mean_cell="L20", denom="$L$20"),
          "DAPT": dict(rows=range(43, 63), dna_rows=range(25, 30), mean_cell="L25", denom="$L$25"),
          "DKK-1": dict(rows=range(66, 86), dna_rows=range(30, 35), mean_cell="L30", denom=None)}
COND = {"Normal KR": 0.0, "Glu2.8": 2.8, "Glu5.5": 5.5, "Glu22": 22.0, "Glu44": 44.0}
tidy, dnas = [], []
for g, d in GROUPS.items():
    mean_dna = V(d["mean_cell"]) if d["mean_cell"] else None
    if g == "DKK-1": mean_dna = V("L30")
    cond_list = list(COND)
    for k, r in enumerate(d["dna_rows"]):
        dnas.append(dict(group=g, condition_mM=COND[cond_list[k]], DNA_ug=V(f"K{r}"), label=V(f"C{r}")))
    for i, r in enumerate(d["rows"]):
        lab = V(f"O{r}")
        assert lab in COND, (g, r, lab)
        tidy.append(dict(group=g, condition_mM=COND[lab], condition=lab, elisa_replicate=(i % 4) + 1, sheet_row=r,
                         OD450=V(f"P{r}"), OD590=V(f"Q{r}"), dOD=V(f"R{r}"), pM_slope0007=V(f"S{r}"), ng_per_300uL=V(f"W{r}"),
                         group_mean_DNA_ug=mean_dna, ng_per_ug_groupmeanDNA=V(f"X{r}"), formula_X=str(wf[f"X{r}"].value)))
# ---- checks -------------------------------------------------------------------------------------------------------------------------
rep = []
def say(s): rep.append(s); print(s)
say("## Raw ELISA audit: " + os.path.basename(XL))
say(f"readings: {len(tidy)} (3 groups x 5 conditions x 4 wells); DNA values: {len(dnas)} (one per group x condition)")
bad = [t for t in tidy if abs(t["ng_per_ug_groupmeanDNA"] - t["ng_per_300uL"] / t["group_mean_DNA_ug"]) > 1e-6 * t["ng_per_ug_groupmeanDNA"]]
say(f"X = W / group-mean DNA holds for {len(tidy) - len(bad)} of {len(tidy)} readings (violations: {len(bad)})")
for g in GROUPS:
    ds = [x["DNA_ug"] for x in dnas if x["group"] == g]
    say(f"{g}: DNA per well (basal, 2.8, 5.5, 22, 44 mM) = {[round(x, 4) for x in ds]}; mean of 5 = {st.mean(ds):.4f}; sheet group mean = {[t for t in tidy if t['group']==g][0]['group_mean_DNA_ug']:.4f}; mean of the four non-basal wells = {st.mean(ds[1:]):.4f}")
# Fig. 1c DNA values
f1 = openpyxl.load_workbook(F1C, data_only=True).worksheets[0]
rows = [r for r in f1.iter_rows(values_only=True)]
hdr = list(rows[0]); it, idn = hdr.index("Treatment"), hdr.index("DNA")
f1c = {}
for r in rows[1:]:
    if r[it] is not None: f1c.setdefault(str(r[it]), []).append(float(r[idn]))
for g in GROUPS:
    a = sorted(round(x, 4) for x in f1c[g]); b = sorted(round(x["DNA_ug"], 4) for x in dnas if x["group"] == g and x["condition_mM"] > 0)
    say(f"Fig. 1c DNA values for {g}: {a} | the four non-basal wells of the ELISA sheet: {b} | identical: {a == b}")
# match with the data file used for the earlier C-peptide analysis
inp = list(csv.DictReader(open(CSVIN, encoding="utf-8")))
say(f"columns of {os.path.basename(CSVIN)}: {list(inp[0].keys())}; rows: {len(inp)}")
# summaries
summ = []
for g in GROUPS:
    for c in sorted(set(t["condition_mM"] for t in tidy)):
        rr = [t for t in tidy if t["group"] == g and t["condition_mM"] == c]
        dna = [x["DNA_ug"] for x in dnas if x["group"] == g and x["condition_mM"] == c][0]
        summ.append(dict(group=g, condition_mM=c, n_elisa_readings=len(rr), dOD_mean=st.mean(t["dOD"] for t in rr), dOD_min=min(t["dOD"] for t in rr), dOD_max=max(t["dOD"] for t in rr),
                         ng300_mean=st.mean(t["ng_per_300uL"] for t in rr), ng300_min=min(t["ng_per_300uL"] for t in rr), ng300_max=max(t["ng_per_300uL"] for t in rr),
                         per_groupmeanDNA_mean=st.mean(t["ng_per_ug_groupmeanDNA"] for t in rr), per_wellDNA_mean=st.mean(t["ng_per_300uL"] for t in rr) / dna,
                         well_DNA_ug=dna, group_mean_DNA_ug=rr[0]["group_mean_DNA_ug"], readings_below_lowest_standard=sum(1 for t in rr if t["dOD"] < 0.044)))
ctrl = {s["condition_mM"]: s for s in summ if s["group"] == "CTRL"}
for s in summ:
    s["ratio_to_CTRL_absolute"] = s["ng300_mean"] / ctrl[s["condition_mM"]]["ng300_mean"]
    s["ratio_to_CTRL_per_groupmeanDNA"] = s["per_groupmeanDNA_mean"] / ctrl[s["condition_mM"]]["per_groupmeanDNA_mean"]
    s["ratio_to_CTRL_dOD"] = s["dOD_mean"] / ctrl[s["condition_mM"]]["dOD_mean"]
say("\nstandards (dOD450-590): blank B0 = " + str(V("E4")) + "; 100 pM = " + str(V("E5")) + "; 200 = " + str(V("E6")) + "; 400 = " + str(V("E7")) + "; 800 = " + str(V("E8")) + "; 1600 = " + str(V("E9")))
say("dOD per pM at the standards: " + ", ".join(f"{p}: {V('E%d' % r) / p:.5f}" for p, r in ((100, 5), (200, 6), (400, 7), (800, 8), (1600, 9))) + "  (workbook uses one slope, 0.0007 per pM, through zero)")
say("conversion chain in the workbook: S = dOD/0.0007 (pM); T = S*5; U = T*300/4; V = U*3.1615e-6; W = V*300 (ng per 300 uL) -> W = 0.3556 * dOD/0.0007 ... (a factor 75 = 300/4 is applied to a concentration)")
say("\n| group | mM | dOD mean (range) | readings < 100-pM standard (dOD 0.044) | ng/300 uL (workbook) mean (range) | per group-mean DNA | well DNA (ug) | group-mean DNA (ug) | ratio to CTRL: absolute / per group-mean DNA |")
say("|---|---|---|---|---|---|---|---|---|")
for s in summ:
    say(f"| {s['group']} | {s['condition_mM']:g} | {s['dOD_mean']:.3f} ({s['dOD_min']:.3f}-{s['dOD_max']:.3f}) | {s['readings_below_lowest_standard']}/4 | {s['ng300_mean']:.1f} ({s['ng300_min']:.1f}-{s['ng300_max']:.1f}) | {s['per_groupmeanDNA_mean']:.1f} | {s['well_DNA_ug']:.4f} | {s['group_mean_DNA_ug']:.4f} | {s['ratio_to_CTRL_absolute']:.2f} / {s['ratio_to_CTRL_per_groupmeanDNA']:.2f} |")
odd = [t for t in tidy if t["OD590"] > 0.05]
say("\nreadings with an unusually high OD590 reference (> 0.05; typical 0.036-0.040): " + "; ".join(f"{t['group']} {t['condition']} rep {t['elisa_replicate']}: OD450 {t['OD450']}, OD590 {t['OD590']}, X = {t['ng_per_ug_groupmeanDNA']:.1f}" for t in odd))
with open(os.path.join(OUT, "Cpeptide_raw_elisa_tidy.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(tidy[0].keys())); w.writeheader(); w.writerows(tidy)
with open(os.path.join(OUT, "Cpeptide_raw_elisa_summary.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(summ[0].keys())); w.writeheader(); w.writerows(summ)
with open(os.path.join(OUT, "Cpeptide_raw_elisa_dna_wells.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(dnas[0].keys())); w.writeheader(); w.writerows(dnas)
open(os.path.join(OUT, "Cpeptide_raw_elisa_audit.md"), "w", encoding="utf-8").write("\n".join(rep) + "\n")
