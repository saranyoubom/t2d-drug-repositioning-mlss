"""
Extracts the melt-peak curves (-d(RFU)/dT versus temperature) of every well from the Bio-Rad CFX Maestro 'Melt Curve Derivative Results'
exports, with gene, treatment and melt temperature (Tm) from the matching 'Melt Curve Peak Results' export, into tidy CSV files.
The raw exports have lower-cased OOXML part names; they are repaired into a scratch copy first (the originals are never touched).
Usage: python -I extract_melt_curves.py <Export_qPCR_dir> <fix_xlsx_case.py dir> <scratch_dir> <out_dir>
Outputs: melt_curves_long.csv (run, well, gene, treatment, temperature, dRFU_dT), melt_peaks.csv (run, well, gene, treatment, Tm, peak_height)
"""
import sys, os, re, csv, collections
import openpyxl
sys.stdout.reconfigure(encoding="utf-8")
EXPORT, FIXDIR, SCR, OUT = sys.argv[1:5]
sys.path.insert(0, FIXDIR)
from fix_xlsx_case import fix_file
os.makedirs(SCR, exist_ok=True); os.makedirs(OUT, exist_ok=True)

RUNS = {"Jun03_Gr1": "2020-Jun-03_Gr1_13genes-Ctrl-DAPT-DKK1", "Jun03_Gr2": "2020-Jun-03_Gr2_13genes-Ctrl-DAPT-DKK1",
        "Jun24": "2020-Jun-24_Round2_13genes-Ctrl-DAPT-DKK1",
        "Aug26": "2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Kcnj11-Cacna1d-Cacna1c-Ptbp1", "Nov09": "2020-11-09_ITPR, RYR1-3"}
TREAT = {"ctrl": "CTRL", "control": "CTRL", "1.ctrl": "CTRL", "dapt": "DAPT", "2. dapt": "DAPT", "dkk-1": "DKK-1", "dkk": "DKK-1", "dkk1": "DKK-1", "3. dkk": "DKK-1"}
GENE = {"Ins-1": "Ins1", "GLUT-2": "Glut2", "ITPR1": "Itpr1", "RYR1": "Ryr1", "RYR2": "Ryr2", "RYR3": "Ryr3", "GAPDH": "Gapdh"}

def norm_well(w):
    m = re.fullmatch(r"([A-Pa-p])0*(\d+)", str(w).strip())
    return f"{m.group(1).upper()}{int(m.group(2))}" if m else str(w).strip()

def clean_gene(t):
    t = re.sub(r"^\d+\.\s*", "", str(t).strip())
    t = re.sub(r"_(Ctrl|DAPT|DKK1|DKK-1)$", "", t)
    return GENE.get(t, t)

def treat_of(sample, target):
    s = str(sample or "").strip().lower()
    if s in TREAT: return TREAT[s]
    m = re.search(r"_(Ctrl|DAPT|DKK1|DKK-1)$", str(target))
    return TREAT[m.group(1).lower()] if m else ""

curves, peaks = [], []
for run, stem in RUNS.items():
    paths = {}
    for kind in ("Derivative", "Peak"):
        src = os.path.join(EXPORT, f"{stem} -  Melt Curve {kind} Results.xlsx")
        dst = os.path.join(SCR, f"{run}__{kind}.xlsx"); fix_file(src, dst); paths[kind] = dst
    pk = openpyxl.load_workbook(paths["Peak"], data_only=True); ws = pk[pk.sheetnames[0]]
    info = {}                                           # well -> (gene, treatment, Tm, height)
    for r in ws.iter_rows(min_row=2, values_only=True):
        if not r[1]: continue
        w = norm_well(r[1]); gene = clean_gene(r[3]); tr = treat_of(r[5], r[3])
        tm = r[6]; ph = r[7]
        try: tm = float(tm)
        except (TypeError, ValueError): tm = None
        try: ph = float(ph)
        except (TypeError, ValueError): ph = None
        info[(gene, w)] = (tr, tm, ph); info.setdefault(("*", w), (gene, tr, tm, ph, str(r[1]).strip()))
    dv = openpyxl.load_workbook(paths["Derivative"], data_only=True)
    for sn in dv.sheetnames:
        if sn == "Run Information": continue
        ws = dv[sn]
        rows = list(ws.iter_rows(values_only=True))
        header = rows[0]; data = rows[1:]
        temps = [r[1] for r in data]
        for j in range(2, len(header)):
            if header[j] is None: continue
            w = norm_well(header[j])
            if sn == "SYBR":                              # wells of the 96-well plate: gene from the peak table
                g = info.get(("*", w))
                if not g: continue
                gene, tr, tm, ph, label = g
            else:                                          # one sheet per gene
                gene = GENE.get(sn, sn)
                v = info.get((gene, w))
                if not v: continue
                tr, tm, ph = v; label = f"{sn}:{header[j]}"      # same well label as the per-well efficiency table
            for t, r in zip(temps, data):
                if t is None or r[j] is None: continue
                curves.append((run, label, w, gene, tr, float(t), float(r[j])))
            peaks.append((run, label, w, gene, tr, tm, ph))
with open(os.path.join(OUT, "melt_curves_long.csv"), "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f); wr.writerow(["run", "well_label", "well", "gene", "treatment", "temperature", "dRFU_dT"]); wr.writerows(curves)
with open(os.path.join(OUT, "melt_peaks.csv"), "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f); wr.writerow(["run", "well_label", "well", "gene", "treatment", "Tm", "peak_height"]); wr.writerows(peaks)
cnt = collections.Counter((p[3]) for p in peaks)
print("wells:", len(peaks), "| curve points:", len(curves), "| genes:", len(cnt))
for g, n in sorted(cnt.items()): print(f"   {g:8s} {n:3d} wells  runs={sorted(set(p[0] for p in peaks if p[3]==g))}  treatments={dict(collections.Counter(p[4] for p in peaks if p[3]==g))}")
