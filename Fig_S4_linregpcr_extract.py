"""
Fig_S4_linregpcr_extract.py
===========================
Post hoc amplification efficiency estimation (LinRegPCR method) for all 21
RT-qPCR gene targets used in this study. Reads per-cycle fluorescence data
from Bio-Rad CFX Manager xlsx exports and outputs three CSV files used by
Fig_S4_linregpcr_plot.R to produce Supplementary Figure 4.

Algorithm
---------
For each well × gene × treatment combination:
  1. Log10-transform baseline-subtracted RFU (shift by min+1 to avoid log(0)).
  2. Scan all 5-cycle sliding windows; retain those where R² >= R2_MIN and
     E = 10^slope >= E_MIN (E_MIN excludes flat baseline-noise windows).
  3. Identify the window whose midpoint is closest to the log-RFU inflection
     point (steepest rise = exponential phase).
  4. Among all wells in the same condition, select the representative well
     whose exponential-phase window has the highest R².

Output CSVs
-----------
  linregpcr_curves.csv   — full log10-RFU traces for representative wells
  linregpcr_windows.csv  — 5-cycle regression fit values for the selected window
  linregpcr_repannot.csv — per-condition E, efficiency (%), R² annotations

Input files required
--------------------
Set BASE to the directory containing the six Bio-Rad xlsx export files:
  2020-Jun-03  : 13 genes (CTRL / DAPT / DKK-1)
  2020-Aug-26  : Kcnj11 + 5 genes (Control / DAPT / DKK)
  2020-Nov-09  : Gapdh + Itpr1 + Ryr2 / Ryr3 (1.Ctrl / 2. DAPT / 3. DKK)

Set OUT to the directory where the three CSV files should be written.

Dependencies: numpy, scipy
"""

import io, sys, zipfile, csv, pathlib, re
import xml.etree.ElementTree as ET
import numpy as np
from scipy.stats import linregress

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# ── User-configurable paths ────────────────────────────────────────────────────
BASE = pathlib.Path("Export_qPCR")          # directory with Bio-Rad xlsx files
OUT  = pathlib.Path(".")                    # directory for output CSVs
# ──────────────────────────────────────────────────────────────────────────────

NS    = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NSREL = "http://schemas.openxmlformats.org/package/2006/relationships"

# ── XML helpers ───────────────────────────────────────────────────────────────
def _get_strs(z):
    si = ET.parse(z.open("xl/sharedstrings.xml")).getroot().findall(f".//{{{NS}}}si")
    return [s.find(f"{{{NS}}}t").text if s.find(f"{{{NS}}}t") is not None else "" for s in si]

def _parse_ws(z, xml_path, strs):
    rows = []
    for row in ET.parse(z.open(xml_path)).getroot().findall(f".//{{{NS}}}row"):
        cells = {}
        for c in row.findall(f"{{{NS}}}c"):
            col = "".join(ch for ch in c.get("r", "") if ch.isalpha())
            v = c.find(f"{{{NS}}}v")
            if v is not None and v.text is not None:
                cells[col] = strs[int(v.text)] if c.get("t") == "s" else v.text
        rows.append(cells)
    return rows

def _sheet_map(z):
    wb   = ET.parse(z.open("xl/workbook.xml")).getroot()
    rels = ET.parse(z.open("xl/_rels/workbook.xml.rels")).getroot()
    rid2tgt = {r.get("Id"): r.get("Target") for r in rels.findall(f"{{{NSREL}}}Relationship")}
    out = {}
    for s in wb.findall(f".//{{{NS}}}sheet"):
        name = s.get(f"{{{NS}}}name") or s.get("name", "?")
        rid  = s.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id") or ""
        tgt  = rid2tgt.get(rid, "")
        if tgt.startswith("worksheets/"): tgt = "xl/" + tgt
        elif tgt.startswith("/xl/"): tgt = tgt[1:]
        out[name] = tgt
    return out

def read_well_map(cq_path, strip_numbers=True):
    with zipfile.ZipFile(str(cq_path)) as z:
        strs = _get_strs(z)
        sm   = _sheet_map(z)
        xml  = sm.get("SYBR") or next(v for k, v in sm.items() if k != "Run Information" and v)
        rows = _parse_ws(z, xml, strs)
    hdr = rows[0]
    rev = {v: k for k, v in hdr.items()}
    wc, tc, sc = rev.get("Well"), rev.get("Target"), rev.get("Sample")
    mapping = {}
    for row in rows[1:]:
        w, t, s = row.get(wc, ""), row.get(tc, ""), row.get(sc, "") or ""
        if not w or not t: continue
        m = re.match(r"([A-H])0?(\d+)", str(w))
        norm_w = (m.group(1) + m.group(2)) if m else str(w)
        gene = re.sub(r"^\d+\.\s*", "", str(t)) if strip_numbers else str(t)
        mapping[norm_w] = (gene, str(s))
    return mapping

def read_amp(amp_path):
    with zipfile.ZipFile(str(amp_path)) as z:
        strs = _get_strs(z)
        sm   = _sheet_map(z)
        cyc_out, rfu_all = None, {}
        for sheet_name, xml_path in sm.items():
            if sheet_name == "Run Information" or not xml_path: continue
            if xml_path not in z.namelist(): continue
            rows = _parse_ws(z, xml_path, strs)
            if not rows: continue
            hdr = rows[0]
            cc = next((c for c, v in hdr.items() if str(v) == "Cycle"), None)
            if cc is None: continue
            wc = {c: str(v) for c, v in hdr.items() if str(v) != "Cycle"}
            cyc_list, rfu_this = [], {w: [] for w in wc.values()}
            for row in rows[1:]:
                if cc not in row: continue
                cyc_list.append(float(row[cc]))
                for col, wpos in wc.items():
                    rfu_this[wpos].append(float(row.get(col, 0) or 0))
            if cyc_out is None and cyc_list:
                cyc_out = np.array(cyc_list)
            for wpos, vals in rfu_this.items():
                rfu_all[wpos] = np.array(vals)
    return cyc_out, rfu_all

# ── LinRegPCR core ────────────────────────────────────────────────────────────
E_MIN   = 1.10   # minimum E to exclude flat baseline-noise windows (E ~ 1.007)
R2_MIN  = 0.99   # minimum R² for a valid exponential-phase window
WIN     = 5      # window width (cycles)

def find_all_windows(rfu_raw, cycles):
    shift = max(0.0, -rfu_raw.min()) + 1.0
    y = np.log10(rfu_raw + shift)
    valid = []
    for i in range(len(cycles) - WIN + 1):
        xw, yw = cycles[i:i + WIN], y[i:i + WIN]
        if np.any(~np.isfinite(yw)): continue
        slope, ic, r, *_ = linregress(xw, yw)
        E = 10 ** slope
        if slope > 0 and r ** 2 >= R2_MIN and E >= E_MIN:
            valid.append((i, slope, ic, r ** 2, E))
    return valid, shift

def inflection_index(rfu_raw, cycles):
    shift = max(0.0, -rfu_raw.min()) + 1.0
    y = np.log10(rfu_raw + shift)
    return int(np.argmax(np.diff(y)))

def best_exp_window(rfu_raw, cycles):
    """Window closest to inflection point; ties broken by highest R²."""
    valid, shift = find_all_windows(rfu_raw, cycles)
    if not valid:
        return None, shift
    infl = inflection_index(rfu_raw, cycles)
    scored = sorted(valid, key=lambda w: (abs(w[0] + (WIN - 1) / 2.0 - infl), -w[3]))
    return scored[0], shift

def select_representative(cond_wells, rfu, cycles):
    """Well whose exponential-phase window has the highest R² (best regression
    in the correct amplification region)."""
    best = (None, None, 1.0, -1.0)   # well, window, shift, r2
    for wp, _ in cond_wells:
        if wp not in rfu: continue
        w, shift = best_exp_window(rfu[wp], cycles)
        if w is None: continue
        if w[3] > best[3]:
            best = (wp, w, shift, w[3])
    return best[0], best[1], best[2]

# ── Load Bio-Rad xlsx exports ─────────────────────────────────────────────────
wmap_gr1 = read_well_map(BASE / "2020-Jun-03_Gr1_13genes-Ctrl-DAPT-DKK1 -  Quantification Cq Results.xlsx")
cyc_gr1, rfu_gr1 = read_amp(BASE / "2020-Jun-03_Gr1_13genes-Ctrl-DAPT-DKK1 -  Quantification Amplification Results.xlsx")

wmap_a26 = read_well_map(
    BASE / "2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Kcnj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Cq Results.xlsx",
    strip_numbers=False)
cyc_a26, rfu_a26 = read_amp(
    BASE / "2020-08-26_INS-1_DAPT-DKK_Tcf7L2-Kcnj11-Cacna1d-Cacna1c-Ptbp1 -  Quantification Amplification Results.xlsx")

wmap_n09 = read_well_map(BASE / "2020-11-09_ITPR, RYR1-3 -  Quantification Cq Results.xlsx")
cyc_n09, rfu_n09 = read_amp(BASE / "2020-11-09_ITPR, RYR1-3 -  Quantification Amplification Results.xlsx")

def g2w_from(wmap):
    d = {}
    for wp, (gene, cond) in wmap.items():
        d.setdefault(gene, []).append((wp, cond))
    return d

g2w_gr1 = g2w_from(wmap_gr1)
g2w_a26 = g2w_from(wmap_a26)
g2w_n09 = g2w_from(wmap_n09)

# ── Gene definitions ──────────────────────────────────────────────────────────
# (map_key, display_label, run_id, [(raw_condition, display_condition), ...])
# Order matches Supplementary Figure 3 (4-column, 6-row grid).
GENES = [
    ("Glp1r",   "Glp1r",   "gr1", [("Ctrl", "CTRL"), ("DAPT", "DAPT"), ("DKK-1", "DKK-1")]),
    ("GLUT-2",  "Glut2",   "gr1", [("Ctrl", "CTRL"), ("DAPT", "DAPT"), ("DKK-1", "DKK-1")]),
    ("Ins-1",   "Ins1",    "gr1", [("Ctrl", "CTRL"), ("DAPT", "DAPT"), ("DKK-1", "DKK-1")]),
    ("Hey1",    "Hey1",    "gr1", [("Ctrl", "CTRL"), ("DAPT", "DAPT"), ("DKK-1", "DKK-1")]),
    ("Hes1",    "Hes1",    "gr1", [("Ctrl", "CTRL"), ("DAPT", "DAPT"), ("DKK-1", "DKK-1")]),
    ("Wnt2",    "Wnt2",    "gr1", [("Ctrl", "CTRL"), ("DAPT", "DAPT"), ("DKK-1", "DKK-1")]),
    ("Wnt2b",   "Wnt2b",   "gr1", [("Ctrl", "CTRL"), ("DAPT", "DAPT"), ("DKK-1", "DKK-1")]),
    ("Wnt5a",   "Wnt5a",   "gr1", [("Ctrl", "CTRL"), ("DAPT", "DAPT"), ("DKK-1", "DKK-1")]),
    ("Wnt5b",   "Wnt5b",   "gr1", [("Ctrl", "CTRL"), ("DAPT", "DAPT"), ("DKK-1", "DKK-1")]),
    ("Wnt9a",   "Wnt9a",   "gr1", [("Ctrl", "CTRL"), ("DAPT", "DAPT"), ("DKK-1", "DKK-1")]),
    ("Lef1",    "Lef1",    "gr1", [("Ctrl", "CTRL"), ("DAPT", "DAPT"), ("DKK-1", "DKK-1")]),
    ("Tcf7",    "Tcf7",    "gr1", [("Ctrl", "CTRL"), ("DAPT", "DAPT"), ("DKK-1", "DKK-1")]),
    ("Tcf7l2",  "Tcf7l2",  "a26", [("Control", "CTRL"), ("DAPT", "DAPT"), ("DKK", "DKK-1")]),
    ("ITPR1",   "Itpr1",   "n09", [("1.Ctrl", "CTRL"), ("2. DAPT", "DAPT"), ("3. DKK", "DKK-1")]),
    ("RYR2",    "Ryr2",    "n09", [("1.Ctrl", "CTRL"), ("2. DAPT", "DAPT"), ("3. DKK", "DKK-1")]),
    ("RYR3",    "Ryr3",    "n09", [("1.Ctrl", "CTRL"), ("2. DAPT", "DAPT"), ("3. DKK", "DKK-1")]),
    ("Cacna1c", "Cacna1c", "a26", [("Control", "CTRL"), ("DAPT", "DAPT"), ("DKK", "DKK-1")]),
    ("Cacna1d", "Cacna1d", "a26", [("Control", "CTRL"), ("DAPT", "DAPT"), ("DKK", "DKK-1")]),
    ("Kcnj11",  "Kcnj11",  "a26", [("Control", "CTRL"), ("DAPT", "DAPT"), ("DKK", "DKK-1")]),
    ("Ptbp1",   "Ptbp1",   "a26", [("Control", "CTRL"), ("DAPT", "DAPT"), ("DKK", "DKK-1")]),
    ("Gapdh",   "Gapdh",   "gr1", [("Ctrl", "CTRL"), ("DAPT", "DAPT"), ("DKK-1", "DKK-1")]),
]

# ── Extract and write CSVs ────────────────────────────────────────────────────
run_map = {
    "gr1": (cyc_gr1, rfu_gr1, g2w_gr1),
    "a26": (cyc_a26, rfu_a26, g2w_a26),
    "n09": (cyc_n09, rfu_n09, g2w_n09),
}

curves_rows, windows_rows, repannot_rows = [], [], []

for gene_key, gene_label, src, cond_pairs in GENES:
    cyc, rfu, g2w = run_map[src]
    wells = g2w.get(gene_key, [])

    for raw_cond, disp_cond in cond_pairs:
        cond_wells = [(wp, c) for wp, c in wells if c == raw_cond and wp in rfu]
        if not cond_wells:
            continue
        best_well, best_win, best_shift = select_representative(cond_wells, rfu, cyc)
        if best_well is None:
            best_well  = cond_wells[0][0]
            best_win   = None
            best_shift = max(0.0, -rfu[best_well].min()) + 1.0

        y_log = np.log10(rfu[best_well] + best_shift)
        for ci, (cv, lv) in enumerate(zip(cyc, y_log)):
            curves_rows.append({
                "gene": gene_label, "condition": disp_cond,
                "cycle": cv, "log10_rfu": round(float(lv), 5)
            })
            if best_win and best_win[0] <= ci <= best_win[0] + WIN - 1:
                si, slope, ic, _, _ = best_win
                fit_y = slope * cv + ic
                windows_rows.append({
                    "gene": gene_label, "condition": disp_cond,
                    "cycle": cv, "fit_log10": round(float(fit_y), 5)
                })

        if best_win:
            _, slope, ic, r2, E_rep = best_win
            repannot_rows.append({
                "gene":      gene_label,
                "condition": disp_cond,
                "E_rep":     round(E_rep, 3),
                "pct_rep":   round((E_rep - 1) * 100, 1),
                "r2":        round(r2, 4),
            })

def write_csv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

write_csv(OUT / "linregpcr_curves.csv",   curves_rows,   ["gene", "condition", "cycle", "log10_rfu"])
write_csv(OUT / "linregpcr_windows.csv",  windows_rows,  ["gene", "condition", "cycle", "fit_log10"])
write_csv(OUT / "linregpcr_repannot.csv", repannot_rows, ["gene", "condition", "E_rep", "pct_rep", "r2"])

print(f"Curves:    {len(curves_rows)} rows")
print(f"Windows:   {len(windows_rows)} rows")
print(f"Rep annot: {len(repannot_rows)} rows")
genes_out = sorted(set(r["gene"] for r in curves_rows))
print(f"Genes ({len(genes_out)}): {genes_out}")
