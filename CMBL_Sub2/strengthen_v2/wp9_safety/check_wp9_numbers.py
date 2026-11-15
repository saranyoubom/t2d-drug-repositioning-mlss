"""check_wp9_numbers.py - independent rebuild of every WP9 number from the saved openFDA responses (results/raw/*.json) and from the pre-specification text.
It has its own parser of the saved responses, its own query builders (not fda_fetch.py) and its own formulas (not wp9_screen.py); wp9_screen.py is imported only to compare its
group definitions and to compare its output files with the rebuilt values. Prints OK or FAIL per check, 'N checks, M failures', and writes check_wp9_last_run.txt (last).
Usage: python -I check_wp9_numbers.py"""
import os, sys, re, csv, json, math, glob, hashlib, urllib.parse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from scipy import stats
import wp9_screen as W

V2 = os.path.dirname(HERE)
RES = os.path.join(HERE, 'results')
RAW = os.path.join(RES, 'raw')
BASE = 'https://api.fda.gov/drug/event.json'
Dd, Rr = 'patient.drug.openfda.substance_name', 'patient.reaction.reactionmeddrapt.exact'
fails, count = [], [0]


def ok(cond, msg):
    count[0] += 1
    print(('OK   ' if cond else 'FAIL ') + msg, flush=True)
    if not cond:
        fails.append(msg)


def close(a, b, tol=6e-4):
    return abs(float(a) - float(b)) <= tol * max(1.0, abs(float(b)))


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


# 1. groups from the pre-specification text
spec = open(os.path.join(HERE, 'PRESPEC_WP9_safety_screen.md'), encoding='utf8').read()
groups = {}
for m in re.finditer(r'^\| ([A-Z]+) \| [^|]+ \| ([^|]+) \|$', spec, flags=re.M):
    groups[m.group(1)] = [t.strip() for t in m.group(2).split(';')]
ok(len(groups) == 8, 'pre-specification table has 8 event groups (%d)' % len(groups))
ok(groups == W.PT, 'group definitions of wp9_screen.py equal the pre-specification table term by term')
ok(sum(len(v) for v in groups.values()) == 40, 'the groups hold 40 preferred terms (%d)' % sum(len(v) for v in groups.values()))

# 2. saved responses
resp = {}
for p in glob.glob(os.path.join(RAW, '*.json')):
    rec = json.load(open(p, encoding='utf8'))
    resp[rec['url']] = rec['response']['meta']['results']['total']
    ok_url = hashlib.sha1(rec['url'].encode()).hexdigest() + '.json' == os.path.basename(p)
    if not ok_url:
        ok(False, 'file name is not the sha1 of its url: ' + os.path.basename(p))


def U(search):
    return BASE + '?' + ('search=' + urllib.parse.quote_plus(search) + '&' if search else '') + 'limit=1'


def T(search):
    return resp[U(search)]


def drug(t):
    return Dd + ':"' + t + '"'


def grp(c):
    return Rr + ':(' + ' OR '.join('"' + x + '"' for x in groups[c]) + ')'


# 3. MLSS inputs
cls = {r['drug']: r['target_gene'] for r in csv.DictReader(open(os.path.join(V2, 'wp3_mlss', 'direction_table.csv'), encoding='utf8')) if r['used_in_mlss'] == 'True'}
cross = []
for r in csv.DictReader(open(os.path.join(V2, 'wp3_mlss', 'mlss_scenario_scores.csv'), encoding='utf8')):
    a, b = r['Combination'].split(' + ')
    if cls[a] != cls[b]:
        cross.append((int(float(r['Rank'])), r['Combination'], a, b))
cross.sort()
ok(len(cls) == 18 and len(cross) == 77 and [c[0] for c in cross] == list(range(1, 78)), '18 drugs, 77 cross-target pairs with ranks 1 to 77')
tm = lambda n: 'nitrendipine' if n == '(S)-nitrendipine' else n.lower()
terms = sorted({tm(n) for n in cls})
ok(len(terms) == 17, '17 distinct search terms')

# 4. rebuild and compare
S = json.load(open(os.path.join(RES, 'summary.json'), encoding='utf8'))
N = T(None)
ok(N == S['N_reports'], 'N all reports %d equals summary' % N)
E = {c: T(grp(c)) for c in groups}
ok(all(E[c] == S['per_group'][c]['E_g'] for c in groups), 'event-group report totals equal summary: ' + ', '.join('%s %d' % (c, E[c]) for c in groups))
nd = {t: T(drug(t)) for t in terms}
ok(all(nd[t] == S['drug_reports'][t] for t in terms), 'drug report totals equal summary')
ev = {t: nd[t] >= 100 for t in terms}
ok(sum(ev.values()) == S['evaluable_drugs'], '%d of 17 drugs evaluable (100 reports or more)' % sum(ev.values()))
a_dg = {(t, c): T(drug(t) + ' AND ' + grp(c)) for t in terms if ev[t] for c in groups}
ok(all(a_dg[(t, c)] <= min(nd[t], E[c]) for (t, c) in a_dg), 'every drug-group count is at most the drug total and the group total')


def lor_ci(a, n, e):
    cells = np.array([a, n - a, e - a, N - n - e + a], dtype=float)
    if cells.min() <= 0:
        cells = cells + 0.5
    lo_r = math.log(cells[0] * cells[3] / (cells[1] * cells[2]))
    se = math.sqrt(float((1.0 / cells).sum()))
    return math.exp(lo_r), math.exp(lo_r - 1.959963984540054 * se), math.exp(lo_r + 1.959963984540054 * se)


rows = list(csv.DictReader(open(os.path.join(RES, 'signals_single.csv'), encoding='utf8')))
ok(len(rows) == sum(ev.values()) * 8, 'signals_single.csv has %d rows (evaluable drugs x 8 groups)' % len(rows))
bad = 0
sig = {}
for r in rows:
    t, c = r['term'], r['group']
    a = a_dg[(t, c)]
    o, lo, hi = lor_ci(a, nd[t], E[c])
    s = a >= 3 and lo > 1
    sig[(t, c)] = s
    if not (int(r['reports_in_group']) == a and int(r['reports_of_drug']) == nd[t] and close(r['ROR'], o) and close(r['ROR_lo95'], lo) and close(r['ROR_hi95'], hi) and int(r['signal']) == int(s)):
        bad += 1
ok(bad == 0, 'every row of signals_single.csv (counts, ROR, interval, signal) equals the independent rebuild (%d rows differ)' % bad)
ok(sum(sig.values()) == S['single_signals_total'], '%d single-drug signals of %d tests equal summary' % (sum(sig.values()), len(rows)))

nab = {}
aab = {}
for (_, _, a, b) in cross:
    k = tuple(sorted((tm(a), tm(b))))
    if k not in nab:
        nab[k] = T(drug(k[0]) + ' AND ' + drug(k[1]))
        if nab[k] >= 20:
            for c in groups:
                aab[(k, c)] = T(drug(k[0]) + ' AND ' + drug(k[1]) + ' AND ' + grp(c))
ok(len(nab) == 70, '70 distinct search-term pairs retrieved')
ok(all(nab[tuple(sorted((tm(a), tm(b))))] <= min(nd[tm(a)], nd[tm(b)]) for (_, _, a, b) in cross), 'co-report counts never exceed either drug total')
pr = {int(r['rank']): r for r in csv.DictReader(open(os.path.join(RES, 'pairs_all77.csv'), encoding='utf8'))}
lg = list(csv.DictReader(open(os.path.join(RES, 'pairs_long.csv'), encoding='utf8')))
bad = 0
nflag = 0
nlong = 0
nevaluable = 0
for (rank, name, a, b) in cross:
    tA, tB = tm(a), tm(b)
    k = tuple(sorted((tA, tB)))
    n_ab = nab[k]
    r = pr[rank]
    status = 'evaluable' if (n_ab >= 20 and ev[tA] and ev[tB]) else ('component drug below 100 reports' if n_ab >= 20 else 'fewer than 20 co-reports')
    sA = (' '.join(c for c in groups if sig.get((tA, c))) or '-') if ev[tA] else 'n/e'
    sB = (' '.join(c for c in groups if sig.get((tB, c))) or '-') if ev[tB] else 'n/e'
    flags = []
    if status == 'evaluable':
        nevaluable += 1
        for c in groups:
            x = aab[(k, c)]
            ci = stats.binomtest(x, n_ab).proportion_ci(confidence_level=0.95, method='exact')
            pA = (a_dg[(tA, c)] - x) / (nd[tA] - n_ab)
            pB = (a_dg[(tB, c)] - x) / (nd[tB] - n_ab)
            fl = x >= 3 and ci.low > max(pA, pB)
            if fl:
                flags.append(c)
                nflag += 1
            nlong += 1
            row = [q for q in lg if int(q['rank']) == rank and q['group'] == c]
            if not (len(row) == 1 and int(row[0]['reports_both_in_group']) == x and int(row[0]['reports_both']) == n_ab and close(row[0]['p_both_lo95'], ci.low, 1e-3) and close(row[0]['p_both_hi95'], ci.high, 1e-3) and int(row[0]['co_report_excess']) == int(fl)):
                bad += 1
    fl_s = (' '.join(flags) or '-') if status == 'evaluable' else 'n/a'
    if not (int(r['reports_listing_both']) == n_ab and r['status'] == status and r['signals_drug_A'] == sA and r['signals_drug_B'] == sB and r['co_report_excess'] == fl_s and r['pair'] == name):
        bad += 1
ok(bad == 0, 'every row of pairs_all77.csv and pairs_long.csv equals the independent rebuild (%d rows differ)' % bad)
ok(nevaluable == S['pairs_evaluable'] and nlong == S['pair_group_tests'] and nflag == S['co_report_flags'] and len(S['co_report_flagged_pairs']) == S['co_report_flagged_pair_count'],
   '%d evaluable pairs, %d pair-group tests and %d co-report flags (in %d pairs) equal summary' % (nevaluable, nlong, nflag, S['co_report_flagged_pair_count']))
ok(len(S['warnings']) == 0, 'analysis raised no consistency warning (%d)' % len(S['warnings']))
top = list(csv.DictReader(open(os.path.join(RES, 'pairs_top15.csv'), encoding='utf8')))
ok([int(t['rank']) for t in top] == list(range(1, 16)), 'pairs_top15.csv holds ranks 1 to 15')
for c in groups:
    t15 = sum(1 for (rank, _, a, b) in cross if rank <= 15 and (sig.get((tm(a), c)) or sig.get((tm(b), c))))
    a77 = sum(1 for (_, _, a, b) in cross if sig.get((tm(a), c)) or sig.get((tm(b), c)))
    ok(t15 == S['per_group'][c]['pairs_top15_with_signalled_drug'] and a77 == S['per_group'][c]['pairs_all77_with_signalled_drug'], 'group %s: pairs with a signalled drug, top 15 %d, all 77 %d' % (c, t15, a77))

# 5. fingerprints
lines = open(os.path.join(V2, 'PROTOCOL_v2_SHA256.txt'), encoding='utf8').read().splitlines()
hashes = [l.split(' *')[0] for l in lines if l.endswith('*wp9_safety/PRESPEC_WP9_safety_screen.md')]
ok(len(hashes) >= 2 and hashes[-1] == sha(os.path.join(HERE, 'PRESPEC_WP9_safety_screen.md')), 'pre-specification file equals its latest fingerprint (%d fingerprints recorded)' % len(hashes))
files = sorted(os.listdir(RAW))
h = hashlib.sha256(''.join('%s:%s\n' % (f, sha(os.path.join(RAW, f))) for f in files).encode()).hexdigest()
rs = open(os.path.join(RES, 'RAW_SHA256.txt'), encoding='utf8').read().split()[0]
ok(h == rs and any(l.startswith(rs) for l in lines),
   'raw files (%d) equal the fingerprint written before the analysis (RAW_SHA256.txt and PROTOCOL_v2_SHA256.txt)' % len(files))
ok(all(f.endswith('.json') for f in files), 'raw folder holds only json files')

# 6. secondary, post hoc free-text fallback (pre-specification section 9)
FBN = {'albiglutide': ['TANZEUM', 'ALBIGLUTIDE'], 'bepridil': ['VASCOR', 'BEPRIDIL'], 'nilvadipine': ['NILVADIPINE'], 'nitrendipine': ['NITRENDIPINE']}
for t, names in FBN.items():
    ok(('%s: %s' % (t, ', '.join(names))) in spec, 'fallback names of %s equal the pre-specification text (%s)' % (t, ', '.join(names)))
import wp9_fallback as FBM
ok(FBM.NAMES == FBN, 'wp9_fallback.py names equal the checker list')
low = [t for t in terms if nd[t] < 100]
ok(sorted(low) == sorted(FBN), 'drugs below 100 reports through substance_name: %s' % ', '.join(sorted(low)))


def tx(t):
    return 'patient.drug.medicinalproduct:(' + ' OR '.join('"%s"' % n for n in FBN[t]) + ')'


F = json.load(open(os.path.join(RES, 'fallback_summary.json'), encoding='utf8'))
fsig = list(csv.DictReader(open(os.path.join(RES, 'fallback_signals_single.csv'), encoding='utf8')))
fpair = list(csv.DictReader(open(os.path.join(RES, 'fallback_pairs.csv'), encoding='utf8')))
bad, n_ev_pairs, n_flags, fsigs, n_fb_rows = 0, 0, 0, {}, 0
for t in sorted(FBN):
    n_t = T(tx(t))
    info = F['drugs'][t]
    if not (info['reports'] == n_t and info['evaluable'] == (n_t >= 100)):
        bad += 1
    if n_t < 100:
        continue
    a_t = {c: T(tx(t) + ' AND ' + grp(c)) for c in groups}
    for c in groups:
        o, lo, hi = lor_ci(a_t[c], n_t, E[c])
        s = a_t[c] >= 3 and lo > 1
        fsigs[(t, c)] = s
        row = [r for r in fsig if r['term'] == t and r['group'] == c]
        if not (len(row) == 1 and int(row[0]['reports_in_group']) == a_t[c] and close(row[0]['ROR'], o) and close(row[0]['ROR_lo95'], lo) and close(row[0]['ROR_hi95'], hi) and int(row[0]['signal']) == int(s)):
            bad += 1
    for (rank, name, a, b) in cross:
        ta, tb = tm(a), tm(b)
        if t not in (ta, tb):
            continue
        p = tb if ta == t else ta
        n_ab = T(tx(t) + ' AND ' + drug(p))
        pev, cev = n_ab >= 20, nd[p] >= 100
        status = 'evaluable' if (pev and cev) else ('partner drug below 100 reports' if pev else 'fewer than 20 co-reports')
        flags = []
        if pev and cev:
            n_ev_pairs += 1
            for c in groups:
                x = T(tx(t) + ' AND ' + drug(p) + ' AND ' + grp(c))
                ci = stats.binomtest(x, n_ab).proportion_ci(confidence_level=0.95, method='exact')
                pA = (a_t[c] - x) / (n_t - n_ab)
                pB = (a_dg[(p, c)] - x) / (nd[p] - n_ab)
                if x >= 3 and ci.low > max(pA, pB):
                    flags.append(c)
                    n_flags += 1
        sA = ' '.join(c for c in groups if fsigs.get((t, c))) or '-'
        sP = (' '.join(c for c in groups if sig.get((p, c))) or '-') if nd[p] >= 100 else 'n/e'
        row = [r for r in fpair if int(r['rank']) == rank and r['fallback_term'] == t]
        n_fb_rows += 1
        if not (len(row) == 1 and int(row[0]['reports_listing_both']) == n_ab and row[0]['status'] == status and row[0]['signals_fallback_drug'] == sA
                and row[0]['signals_partner_drug'] == sP and row[0]['co_report_excess'] == ((' '.join(flags) or '-') if status == 'evaluable' else 'n/a')):
            bad += 1
ok(bad == 0, 'every secondary (fallback) row equals the independent rebuild (%d differ; %d pair rows)' % (bad, n_fb_rows))
ok(F['pairs_evaluable'] == n_ev_pairs and F['co_report_flags'] == n_flags and F['pairs'] == len(fpair) == n_fb_rows,
   'secondary: %d pair rows, %d evaluable pairs and %d co-report flags equal fallback_summary' % (n_fb_rows, n_ev_pairs, n_flags))

# 7. numbers stated in README.md, rebuilt from the saved responses
RD = open(os.path.join(HERE, 'README.md'), encoding='utf8').read()


def has(txt, msg):
    ok(txt in RD, 'README states ' + msg + ': "' + txt + '"')


has('{:,}'.format(N), 'the number of reports in the database')
has(S['last_updated'], 'the snapshot date')
has('%d responses' % len(files), 'the number of saved responses')
has('**%d of %d terms are evaluable**' % (sum(ev.values()), len(terms)), 'evaluable terms')
has('**%d of %d drug-group tests are signals**' % (sum(sig.values()), len(rows)), 'signal count')
for t in terms:
    if ev[t]:
        has('{:,}'.format(nd[t]), 'report total of ' + t)
bad_terms = [t for t in terms if not ev[t]]
ok(sorted(bad_terms) == ['albiglutide', 'bepridil', 'nilvadipine', 'nitrendipine'], 'not evaluable terms are albiglutide, bepridil, nilvadipine, nitrendipine')
inv_pairs = [(rank, name) for (rank, name, a, b) in cross if tm(a) in bad_terms or tm(b) in bad_terms]
has('**%d of the 77 pairs contain such a drug, and %d of the top 15**' % (len(inv_pairs), sum(1 for r, _ in inv_pairs if r <= 15)), 'pairs with an invisible drug')
nn = [nab[tuple(sorted((tm(a), tm(b))))] for (_, _, a, b) in cross]
has('**%d have no report listing both drugs, %d have 1 to 19, and %d have 20 or more**' % (sum(1 for x in nn if x == 0), sum(1 for x in nn if 1 <= x <= 19), sum(1 for x in nn if x >= 20)), 'co-report distribution')
evr = [(rank, nab[tuple(sorted((tm(a), tm(b))))]) for (rank, _, a, b) in cross if nab[tuple(sorted((tm(a), tm(b))))] >= 20]
has('ranks ' + ', '.join(str(r) for r, _ in evr[:-1]) + ' and ' + str(evr[-1][0]), 'ranks of the evaluable pairs')
has('%s to %s co-reports, median %d' % (str(min(x for _, x in evr)), '{:,}'.format(max(x for _, x in evr)), int(np.median([x for _, x in evr]))), 'range and median of the co-reports')
top15n = [nab[tuple(sorted((tm(a), tm(b))))] for (rank, _, a, b) in cross if rank <= 15]
has('**Of the top 15 pairs, %d have fewer than 20 co-reports (%d have none' % (sum(1 for x in top15n if x < 20), sum(1 for x in top15n if x == 0)), 'top 15 co-reports')
has('**%d flags in %d of the %d evaluable pairs, among %d pair-group tests; none in the top 15**' % (S['co_report_flags'], S['co_report_flagged_pair_count'], S['pairs_evaluable'], S['pair_group_tests']), 'flag counts')
for rk, nm, g in S['co_report_flagged_pairs']:
    has('rank %d %s (%s)' % (rk, nm, g.replace(' ', ', ')), 'flagged pair')
grp_counts = sorted([(sum(1 for t in terms if sig.get((t, c))), c) for c in groups], key=lambda z: (-z[0], ['PANC', 'HYPO', 'HYPER', 'GI', 'BILI', 'HYPOT', 'CARD', 'OEDEMA'].index(z[1])))
has('By group, drugs with a signal: ' + ', '.join('%s %d' % (c, k) for k, c in grp_counts), 'drugs with a signal per group')
has(', '.join('%s %d / %d' % (c, S['per_group'][c]['pairs_top15_with_signalled_drug'], S['per_group'][c]['pairs_all77_with_signalled_drug']) for c in groups), 'pairs with a signalled drug (top 15 / all 77)')
by_drug = {t: sum(1 for c in groups if sig.get((t, c))) for t in terms if ev[t]}
has('8 of 8 groups for ' + ', '.join(sorted(t for t, k in by_drug.items() if k == 8)[:-1]) + ' and ' + sorted(t for t, k in by_drug.items() if k == 8)[-1], 'drugs with 8 of 8')
ok(sorted(t for t, k in by_drug.items() if k == 5) == ['dulaglutide', 'exenatide', 'isradipine', 'liraglutide', 'nisoldipine', 'semaglutide'] and sorted(t for t, k in by_drug.items() if k == 3) == ['clevidipine', 'lixisenatide'] and sorted(t for t, k in by_drug.items() if k == 2) == ['nimodipine'], 'drugs with 5, 3 and 2 signalled groups as stated')


def rr(t, c):
    return lor_ci(a_dg[(t, c)], nd[t], E[c])


def fmt(t, c, d):
    o, lo, hi = rr(t, c)
    return ('%.' + str(d) + 'f (%.' + str(d) + 'f to %.' + str(d) + 'f)') % (o, lo, hi)


has('ROR %s (%s' % (fmt('liraglutide', 'PANC', 1).split(' (')[0], fmt('liraglutide', 'PANC', 1).split(' (')[1].replace(' to ', ' to ')), 'liraglutide PANC ROR')
has(fmt('liraglutide', 'PANC', 1), 'liraglutide PANC interval')
has('glucagon HYPO %d of %s, ROR %s' % (a_dg[('glucagon', 'HYPO')], '{:,}'.format(nd['glucagon']), fmt('glucagon', 'HYPO', 1)), 'glucagon HYPO')
has('amlodipine OEDEMA ROR ' + fmt('amlodipine', 'OEDEMA', 2), 'amlodipine OEDEMA')
has('HYPOT ' + fmt('amlodipine', 'HYPOT', 2), 'amlodipine HYPOT')
has('HYPOT ' + fmt('isradipine', 'HYPOT', 1), 'isradipine HYPOT')
has('OEDEMA ' + fmt('isradipine', 'OEDEMA', 1), 'isradipine OEDEMA')
has('HYPOT ' + fmt('liraglutide', 'HYPOT', 2), 'liraglutide HYPOT')
has('OEDEMA ' + fmt('liraglutide', 'OEDEMA', 2), 'liraglutide OEDEMA')
has('%s of %s reports' % ('{:,}'.format(a_dg[('liraglutide', 'PANC')]), '{:,}'.format(nd['liraglutide'])), 'liraglutide PANC counts')
# top-15 table rows
for r in top:
    ca = r['signals_drug_A']
    cb = r['signals_drug_B']
    ca = 'all 8 groups' if len(ca.split()) == 8 else ca
    cb = 'all 8 groups' if len(cb.split()) == 8 else cb
    fl = 'n/a' if r['co_report_excess'] == 'n/a' else ('none' if r['co_report_excess'] == '-' else r['co_report_excess'])
    has('| %s | %s | %s | %s | %s | %s | %s |' % (r['rank'], r['pair'], r['MLSS'], r['reports_listing_both'], ca, cb, fl), 'top-15 row %s' % r['rank'])
# secondary
fb_n = {t: T(tx(t)) for t in FBN}
for t in FBN:
    has('{:,}'.format(fb_n[t]), 'free-text report total of ' + t)
fbrow = {(r['term'], r['group']): r for r in fsig}
for (t, c, d) in [('albiglutide', 'HYPO', 2), ('albiglutide', 'HYPER', 2), ('albiglutide', 'PANC', 2), ('albiglutide', 'GI', 2), ('bepridil', 'HYPOT', 2), ('nilvadipine', 'OEDEMA', 2), ('nilvadipine', 'CARD', 2), ('nitrendipine', 'HYPOT', 2), ('nitrendipine', 'BILI', 2), ('nitrendipine', 'OEDEMA', 2)]:
    has('%s (%s)' % (c, ('%.' + str(d) + 'f') % float(fbrow[(t, c)]['ROR'])), 'free-text ROR %s %s' % (t, c))
has('CARD (%.1f)' % float(fbrow[('bepridil', 'CARD')]['ROR']), 'bepridil CARD')
has('CARD (%.1f)' % float(fbrow[('nitrendipine', 'CARD')]['ROR']), 'nitrendipine CARD')
ok(not fsigs.get(('nitrendipine', 'HYPO')) and all(fsigs.get(('nitrendipine', c)) for c in groups if c != 'HYPO'), 'nitrendipine has signals in 7 of 8 groups, all but HYPO')
has('Of the %d pairs that contain one of these four drugs, **one has 20 or more co-reports (albiglutide + amlodipine, rank 21, 86 co-reports; excess flag GI)**' % len(fpair), 'secondary pair counts')
ok([(r['rank'], r['pair'], r['reports_listing_both'], r['co_report_excess']) for r in fpair if int(r['reports_listing_both']) >= 20] == [('21', 'albiglutide + amlodipine', '86', 'GI')], 'the only secondary pair with 20 or more co-reports is rank 21 albiglutide + amlodipine, 86, flag GI')
sec15 = {(r['rank'], r['fallback_term']): int(r['reports_listing_both']) for r in fpair if int(r['rank']) <= 15}
ok(sec15[('4', 'nitrendipine')] == 5 and sec15[('11', 'nitrendipine')] == 4 and sec15[('15', 'albiglutide')] == 3 and max(v for k, v in sec15.items() if k not in [('4', 'nitrendipine'), ('8', 'nitrendipine'), ('11', 'nitrendipine'), ('12', 'nitrendipine'), ('15', 'albiglutide')]) == 0, 'secondary co-reports in the top 15: nitrendipine + liraglutide 5, nitrendipine + semaglutide 4, albiglutide + nifedipine 3, all others 0')

# 8. WORDING_safety_screen.md: every number of sections 1 to 7 is a number of the README, and the framing words stay inside negations
WD = open(os.path.join(HERE, 'WORDING_safety_screen.md'), encoding='utf8').read()
body = WD.split('## 8. What not to write')[0]
tok = lambda x: set(re.findall(r'(?<![\w.])\d+(?:\.\d+)?', re.sub(r'(?<=\d),(?=\d{3})', '', x)))
rd_tok = tok(RD)
NONRESULT = {'256', '350', '5.5'}      # SHA-256, the 350-word abstract limit, 'part 5.5' of Additional file 5
missing = sorted(t for t in tok(body) if (len(t) >= 2 or '.' in t) and t not in rd_tok and t not in NONRESULT)
ok(not missing, 'every number of the wording (sections 1 to 7) is a number of the README (not found: %s)' % ', '.join(missing))
bad_sent = []
for sent in re.split(r'(?<=[.;:])\s+|\n', body):
    if re.search(r'\b(predict\w*|safe|unsafe|risk|validated|confirms)\b', sent, flags=re.I) and not re.search(r'\b(not|no|never|neither|nor|without)\b', sent, flags=re.I):
        bad_sent.append(sent[:90])
ok(not bad_sent, 'framing words (predict, safe, unsafe, risk, validated, confirms) only inside negations in sections 1 to 7 (offending: %s)' % ' | '.join(bad_sent))
wc = len(re.findall(r'\S+', body.split('## 5. Abstract fragments')[0]))
ok(len(re.findall(r'\S+', re.search(r'## 2\. Results.*?\n(.*?)\n\n## 3', body, flags=re.S).group(1))) <= 330, 'the Results text is at most 330 words')

print('\n%d checks, %d failures' % (count[0], len(fails)))
with open(os.path.join(HERE, 'check_wp9_last_run.txt'), 'w', encoding='utf8') as f:
    f.write('%d checks, %d failures%s\n' % (count[0], len(fails), (': ' + '; '.join(fails)) if fails else ''))
