"""
R2C3: Attempts to reconstruct which raw well corresponds to which technical
replicate label in the already-published fold-change values, by matching
implied dCt values (from raw Cq records) against the ground-truth dCt values
backed out of the published fold-change table.

This is evidence for a disclosed manuscript limitation, not a figure/table
generator: it demonstrates that the well-to-replicate correspondence needed
to retrospectively apply the efficiency correction (R2C3_run_sliwin_all.R)
to the already-published fold-changes cannot be reconstructed for most
genes, because it was not recorded and preserved anywhere that survives.

Inputs required (not included in this deposit):
  ./raw_expression_data_delta_ct.xlsx   Published per-replicate fold-change
                                        table ("Raw_expression_data" sheet:
                                        Treatment column + one column per
                                        gene of 2^-ddCt values), reproduces
                                        Supplementary Data 1 exactly.
  ./all_cq.csv                         Per-well Run/Well/Gene/Treatment/Cq,
                                        extracted from the raw CFX
                                        "Quantification Cq Results" exports
                                        (see script header note below; not
                                        included -- raw instrument exports).

Output: resolved_replicate_mapping.csv -- the well/replicate pairs that
could be matched within tolerance, and a console summary of how many
gene x treatment groups were fully resolved.
"""
import csv, math, sys, io
from collections import defaultdict
import openpyxl
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

GROUND_TRUTH_XLSX = "./raw_expression_data_delta_ct.xlsx"
ALL_CQ_CSV = "./all_cq.csv"
OUT_CSV = "./resolved_replicate_mapping.csv"

wb = openpyxl.load_workbook(GROUND_TRUTH_XLSX, read_only=True, data_only=True)
ws = wb['Raw_expression_data']
rows = list(ws.iter_rows(min_row=1, values_only=True))
header = rows[0]
genes = header[1:]

# ground_truth[gene][treatment] = list of dCt values (one per replicate, in file order)
ground_truth = defaultdict(lambda: defaultdict(list))
for row in rows[1:]:
    treat = row[0].replace('-', '').upper()
    treat = 'DKK1' if 'DKK' in treat else treat
    for gi, gene in enumerate(genes):
        val = row[gi+1]
        if val is None:
            continue
        val = float(val)
        dct = -math.log2(val)
        ground_truth[gene][treat].append(dct)

# load Cq records
cq = defaultdict(list)  # (gene, treatment) -> list of (run, well, Cq)
gapdh_cq = defaultdict(list)  # (run, treatment) -> list of (well, Cq)
with open(ALL_CQ_CSV, encoding='utf-8') as f:
    for row in csv.DictReader(f):
        run, well, gene, treat, cqval = row['Run'], row['Well'], row['Gene'], row['Treatment'], float(row['Cq'])
        if gene == 'Gapdh':
            gapdh_cq[(run, treat)].append((well, cqval))
        cq[(gene, treat)].append((run, well, cqval))

TOL = 0.02  # log2 units; Cq differences preserved to ~6 sig figs in source, so exact matches expected

results = {}  # (gene,treat) -> list of (run, target_well, target_cq, gapdh_well, gapdh_cq, matched_gt_dct)
unresolved = []

for gene in genes:
    for treat in ['CTRL', 'DAPT', 'DKK1']:
        gts = ground_truth[gene].get(treat, [])
        if not gts:
            continue
        candidates = cq.get((gene, treat), [])
        # Build all valid (candidate_idx, gt_idx) edges within tolerance, then
        # greedily assign best (closest) matches first so no well or ground-
        # truth replicate slot is used twice.
        all_pairs = []  # (run, well, tcq, gwell, gcq, dct)
        for run, well, tcq in candidates:
            for gwell, gcq in gapdh_cq.get((run, treat), []):
                all_pairs.append((run, well, tcq, gwell, gcq, tcq - gcq))
        edges = []
        for pi, (run, well, tcq, gwell, gcq, dct) in enumerate(all_pairs):
            for gi, gtval in enumerate(gts):
                diff = abs(dct - gtval)
                if diff < TOL:
                    edges.append((diff, pi, gi))
        edges.sort()
        used_gt_idx = set()
        used_wells = set()
        matches = []
        for diff, pi, gi in edges:
            if gi in used_gt_idx:
                continue
            run, well, tcq, gwell, gcq, dct = all_pairs[pi]
            if (run, well) in used_wells or (run, gwell) in used_wells:
                continue
            matches.append((run, well, tcq, gwell, gcq, gi, dct, gts[gi]))
            used_gt_idx.add(gi)
            used_wells.add((run, well))
            used_wells.add((run, gwell))
        results[(gene, treat)] = matches
        if len(matches) < len(gts):
            unresolved.append((gene, treat, len(matches), len(gts)))

print(f"Total (gene,treatment) groups: {sum(1 for g in genes for t in ['CTRL','DAPT','DKK1'] if ground_truth[g].get(t))}")
print(f"Fully resolved (all replicates matched): {sum(1 for k,v in results.items() if len(v)==len(ground_truth[k[0]][k[1]]))}")
print(f"Unresolved/partial: {len(unresolved)}")
for u in unresolved:
    print("  PARTIAL:", u)

with open(OUT_CSV, "w", newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(["Gene","Treatment","Run","TargetWell","TargetCq","GapdhWell","GapdhCq","ReplicateIdx","MatchedDct","GroundTruthDct"])
    for (gene,treat), matches in results.items():
        for run, well, tcq, gwell, gcq, gi, dct, gtval in matches:
            w.writerow([gene,treat,run,well,tcq,gwell,gcq,gi,dct,gtval])
print(f"\nSaved resolved mapping to {OUT_CSV}")
