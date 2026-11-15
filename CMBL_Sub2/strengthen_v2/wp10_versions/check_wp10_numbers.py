"""check_wp10_numbers.py - independent rebuild of the WP10 numbers from the raw inputs with the csv module only (no pandas, no scipy, no import of the analysis scripts):
the STRING classification (v12.0 export of the reported run, v12.5 network from the API), the DrugCentral counts and potencies (old file, new file), the MLSS
combination counts, and the rank statements taken from results/mlss_comparison.csv; then compares them with the result files and with the numbers stated in README.md.
Prints OK or FAIL per check, '%d checks, %d failures', and writes check_wp10_last_run.txt last. Usage: python -I check_wp10_numbers.py"""
import os, csv, json, gzip, io, hashlib, re, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.dirname(HERE)
PROJECT = os.path.dirname(os.path.dirname(V2))
R = os.path.join(os.path.dirname(V2), 'R_17genes')
RES = os.path.join(HERE, 'results')
fails, count = [], [0]


def ok(c, m):
    count[0] += 1
    print(('OK   ' if c else 'FAIL ') + m, flush=True)
    if not c:
        fails.append(m)


rd = open(os.path.join(HERE, 'README.md'), encoding='utf8').read()


def has(t, m):
    ok(t in rd, 'README states %s: "%s"' % (m, t))


def pn(a, b):
    return '-'.join(sorted([a, b]))


# ---------------- STRING
old = list(csv.DictReader(open(os.path.join(R, 'string_interactions.csv'), encoding='utf-8-sig', newline='')))
new_lines = open(os.path.join(RES, 'string_v12_5_network.tsv'), encoding='utf8').read().strip().splitlines()
hdr = new_lines[0].split('\t')
new = [dict(zip(hdr, l.split('\t'))) for l in new_lines[1:]]
cons = {pn(*r['Consensus_Pairs'].split('_')[:2]) for r in csv.DictReader(open(os.path.join(R, 'Fig_4d_consensus_pairs.csv'), encoding='utf-8-sig', newline=''))}
with open(os.path.join(R, 'Fig_4e_Hes1_significant_correlations_only.csv'), encoding='utf-8-sig', newline='') as f:
    rr = csv.reader(f)
    head = next(rr)[1:]
    pairs = {}
    for row in rr:
        for g2, v in zip(head, row[1:]):
            if v not in ('', 'NA') and row[0] != g2 and abs(float(v)) >= 0.7:
                pairs[pn(row[0], g2)] = float(v)
ok(len(pairs) == 120 and len(cons) == 21, '120 pairs with |r| >= 0.7 and 21 consensus pairs')
ALIAS = {'Slc2a2': 'Glut2'}


def sset(edges, alias):
    out = set()
    for a, b in edges:
        if alias:
            a, b = ALIAS.get(a, a), ALIAS.get(b, b)
        out.add(pn(a, b))
    return out


def cls(ss):
    c = {'Known Validated': 0, 'Novel Discovery': 0, 'STRING Only': 0, 'Not Validated': 0}
    for p in pairs:
        i, v = p in ss, p in cons
        c['Known Validated' if i and v else 'Novel Discovery' if v else 'STRING Only' if i else 'Not Validated'] += 1
    return c


o_e = [(r['node1'], r['node2']) for r in old]
n_e = [(r['preferredName_A'], r['preferredName_B']) for r in new]
res = {}
for name, al in [('lit', False), ('alias', True)]:
    for v, e in [('v120', o_e), ('v125', n_e)]:
        c = cls(sset(e, al))
        res[(name, v)] = c
SM = json.load(open(os.path.join(RES, 'string_summary.json'), encoding='utf8'))
for name, key in [('lit', 'as_reported'), ('alias', 'alias')]:
    for v, k2 in [('v120', 'v12_0'), ('v125', 'v12_5')]:
        c = res[(name, v)]
        s = SM[key][k2]
        ok(all(c[k] == s[k] for k in c) and c['Known Validated'] + c['STRING Only'] == s['in_STRING'], 'STRING %s %s counts equal the result file: %s' % (key, k2, c))
ok(len(sset(o_e, False)) == 59 and len(sset(n_e, False)) == 59, '59 distinct pairs in the v12.0 export and 59 in the v12.5 network')
g17 = {'Wnt5a', 'Wnt5b', 'Wnt9a', 'Lef1', 'Tcf7', 'Tcf7l2', 'Hes1', 'Itpr1', 'Ryr2', 'Ryr3', 'Cacna1c', 'Cacna1d', 'Kcnj11', 'Ptbp1', 'Glut2', 'Ins1', 'Glp1r'}
e17 = lambda E: {pn(ALIAS.get(a, a), ALIAS.get(b, b)) for a, b in E if ALIAS.get(a, a) in g17 and ALIAS.get(b, b) in g17}
a0, a1 = e17(o_e), e17(n_e)
ok(len(a0) == 38 and len(a1) == 39 and sorted(a1 - a0) == ['Kcnj11-Ryr2', 'Wnt5a-Wnt9a', 'Wnt5b-Wnt9a'] and sorted(a0 - a1) == ['Cacna1c-Itpr1', 'Cacna1d-Itpr1'], 'edges among the 17 genes: 38 in 12.0, 39 in 12.5, 3 gained, 2 lost')
glut = [p for p in pairs if 'Glut2' in p.split('-')]
ok(len(glut) == 15 and sorted(p for p in glut if p in sset(o_e, True)) == ['Glp1r-Glut2', 'Glut2-Ins1', 'Glut2-Kcnj11', 'Glut2-Tcf7l2'] and not [p for p in glut if p in sset(o_e, False)], 'Glut2: 15 pairs, 4 of them in STRING only under the alias, none under literal names')
for t in ['Known Validated 5', '34 in STRING (28.3%)']:
    pass
has('29 (24.2%) | 70 (58.3%) | 34 (28.3%)', 'reported STRING classification')
has('30 (25.0%) | 69 (57.5%) | 35 (29.2%)', 'v12.5 classification')
has('33 (27.5%) | 66 (55.0%) | 38 (31.7%)', 'v12.0 alias classification')
has('34 (28.3%) | 65 (54.2%) | 39 (32.5%)', 'v12.5 alias classification')
has('38 in 12.0 and 39 in 12.5', 'edges among 17 genes')
has('Glp1r-Glut2, Glut2-Ins1, Glut2-Kcnj11, Glut2-Tcf7l2; 15 of the 120 pairs involve Glut2', 'Glut2 pairs')

# ---------------- DrugCentral
def read_dc(path, gz):
    raw = open(path, 'rb').read()
    txt = (gzip.decompress(raw) if gz else raw).decode('utf-8-sig')
    rows = list(csv.DictReader(io.StringIO(txt), delimiter='\t' if gz else ','))
    return rows, hashlib.sha256(raw).hexdigest()


odc, _ = read_dc(os.path.join(R, 'drug.target.interaction.csv'), False)
ndc, nsha = read_dc(os.path.join(PROJECT, 'Database', 'drug.target.interaction.tsv.gz'), True)
ok(nsha.startswith('19086849'), 'new DrugCentral file sha256 starts with 19086849')
ok(len(odc) == 19378 and len(ndc) == 22364, 'rows: 19378 before and 22364 after')
ok(len({r['DRUG_NAME'] for r in odc}) == 2587 and len({r['DRUG_NAME'] for r in ndc}) == 2897, 'distinct drugs: 2587 before and 2897 after')
hubs = [r['Gene'].upper() for r in csv.DictReader(open(os.path.join(R, 'Fig_6b_Gene_Centrality_Summary.csv'), encoding='utf-8-sig')) if float(r['Degree']) >= 5]
ok(sorted(hubs) == ['CACNA1D', 'GLP1R', 'PTBP1', 'TCF7L2'], 'hubs (degree >= 5): ' + ', '.join(sorted(hubs)))


def targets(rows):
    t = {}
    for r in rows:
        g = r['GENE'].upper()
        if g in hubs and r['ORGANISM'] == 'Homo sapiens':
            t.setdefault(r['DRUG_NAME'], [])
            if g not in t[r['DRUG_NAME']]:
                t[r['DRUG_NAME']].append(g)
    return t


def pot(rows):
    p = {}
    for r in rows:
        g = r['GENE'].upper()
        if g in hubs and r['ACT_TYPE'] in ('IC50', 'Ki', 'EC50') and r['ORGANISM'] == 'Homo sapiens' and r['RELATION'] == '=' and r['ACT_VALUE'] not in ('', 'NA'):
            k = (r['DRUG_NAME'], g)
            p[k] = max(p.get(k, -1e9), float(r['ACT_VALUE']))
    return p


t0, t1 = targets(odc), targets(ndc)
DS = json.load(open(os.path.join(RES, 'drugcentral_summary.json'), encoding='utf8'))
ok(len(t0) == 18 and len(t1) == 19 and sorted(set(t1) - set(t0)) == ['tirzepatide'] and not (set(t0) - set(t1)), 'MLSS drugs 18 before, 19 after, only tirzepatide added, none removed')
cnt = lambda t, g: sum(1 for d, gs in t.items() if g in gs)
ok((cnt(t0, 'GLP1R'), cnt(t1, 'GLP1R'), cnt(t0, 'CACNA1D'), cnt(t1, 'CACNA1D'), cnt(t0, 'TCF7L2') + cnt(t1, 'TCF7L2') + cnt(t0, 'PTBP1') + cnt(t1, 'PTBP1')) == (7, 8, 11, 11, 0), 'GLP1R drugs 7 then 8, CACNA1D drugs 11 then 11, none for TCF7L2 and PTBP1')
ok(all(len(v) == 1 for v in list(t0.values()) + list(t1.values())), 'every drug has exactly one hub target')
n0, n1 = len(t0), len(t1)
cross0 = sum(1 for a, b in itertools.combinations(t0, 2) if t0[a] != t0[b])
cross1 = sum(1 for a, b in itertools.combinations(t1, 2) if t1[a] != t1[b])
ok((n0 * (n0 - 1) // 2, n1 * (n1 - 1) // 2, cross0, cross1) == (153, 171, 77, 88), 'combinations 153 then 171; cross-target 77 then 88')
p0, p1 = pot(odc), pot(ndc)
diff = sorted(k for k in set(p0) | set(p1) if k[0] in set(t0) | set(t1) and abs(p0.get(k, -99) - p1.get(k, -99)) > 1e-9)
ok(diff == [('tirzepatide', 'GLP1R')] and abs(p1[('tirzepatide', 'GLP1R')] - 8.37) < 1e-9, 'the only potency that differs between the files is tirzepatide at GLP1R = 8.37')
for k, v in [('mlss_drugs_old', 18), ('mlss_drugs_new', 19), ('combinations_old', 153), ('combinations_new', 171), ('cross_target_old', 77), ('cross_target_new', 88), ('combinations_shared', 153), ('score_changed_shared', 0), ('top15_old_still_in_top15_new', 14), ('combinations_with_new_drugs', 18), ('pic50_changes', 1), ('rows_old', 19378), ('rows_new', 22364)]:
    ok(DS[k] == v, 'drugcentral_summary.json %s = %s' % (k, v))
ok(DS['spearman_scores_shared'] == 1.0 and DS['spearman_ranks_shared'] == 1.0 and not DS['tie_straddles_rank15_old'] and not DS['tie_straddles_rank15_new'], 'Spearman 1.0 for scores and ranks of the shared combinations; no tie across rank 15')
cmp = list(csv.DictReader(open(os.path.join(RES, 'mlss_comparison.csv'), encoding='utf-8')))
byc = {r['combination']: r for r in cmp}
old_top = sorted([r for r in cmp if r['rank_old'] and int(r['rank_old']) <= 15], key=lambda r: int(r['rank_old']))
ok(len(old_top) == 15 and sum(1 for r in old_top if int(r['rank_new']) <= 15) == 14, '14 of the 15 reported top combinations are in the new top 15')
ok(all(int(r['rank_new']) == int(r['rank_old']) for r in old_top if int(r['rank_old']) <= 12) and [int(r['rank_new']) for r in old_top[12:]] == [14, 15, 16], 'old ranks 1 to 12 unchanged; old ranks 13, 14, 15 become 14, 15, 16')
ok(byc['isradipine + tirzepatide']['rank_new'] == '13' and byc['isradipine + tirzepatide']['status'] == 'new only' and abs(float(byc['isradipine + tirzepatide']['MLSS_new']) - 9.9158) < 5e-5, 'isradipine + tirzepatide enters at rank 13 with MLSS 9.9158')
tz = sorted(int(r['rank_new']) for r in cmp if 'tirzepatide' in r['combination'])
ok(len(tz) == 18 and tz[:11] == [13, 26, 27, 30, 33, 69, 70, 71, 72, 73, 74] and min(tz[11:]) == 99, 'tirzepatide combinations: 18; cross-target ranks 13, 26, 27, 30, 33, 69 to 74; same-target from rank 99')
same = [r for r in cmp if r['status'] == 'both' and r['MLSS_old'] != '' and abs(float(r['MLSS_old']) - float(r['MLSS_new'])) > 1e-9]
ok(not same, 'no shared combination changed its score')
for t, w in [('Spearman rho 1.0', 'rank correlation'), ('isradipine + tirzepatide enters at rank 13 (MLSS 9.9158)', 'new top-15 entrant'), ('**14 of the 15 reported combinations stay in the top 15**', 'top-15 overlap'), ('153 become 171; cross-target pairs 77 become 88', 'combination counts'),
             ('18 become 19', 'drug count'), ('(8.37)', 'tirzepatide potency'), ('22364 rows against 19378', 'rows'), ('2897 distinct drugs against 2587', 'distinct drugs'),
             ('11 cross-target at ranks 13, 26, 27, 30, 33, 69, 70, 71, 72, 73 and 74; 7 same-target at rank 99 or lower', 'tirzepatide combination ranks'), ('0 changed', 'changed scores')]:
    has(t, w)
ok(len(open(os.path.join(HERE, 'FINGERPRINTS.txt'), encoding='utf8').read().splitlines()) == 10, 'FINGERPRINTS.txt holds 5 fingerprints')
fp = open(os.path.join(HERE, 'FINGERPRINTS.txt'), encoding='utf8').read()
ok(hashlib.sha256(open(os.path.join(HERE, 'PRESPEC_WP10_database_versions.md'), 'rb').read()).hexdigest() in fp, 'the pre-specification equals its fingerprint (unchanged since it was fingerprinted)')
print('\n%d checks, %d failures' % (count[0], len(fails)))
open(os.path.join(HERE, 'check_wp10_last_run.txt'), 'w', encoding='utf8').write('%d checks, %d failures%s\n' % (count[0], len(fails), (': ' + '; '.join(fails)) if fails else ''))
