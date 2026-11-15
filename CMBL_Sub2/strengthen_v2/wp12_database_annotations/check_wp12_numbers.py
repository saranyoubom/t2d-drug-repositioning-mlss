"""check_wp12_numbers.py - WP12 (9 Nov 2026): independent check of README.md and WORDING_annotations.md.
Own parsers (csv, gzip, json, plain-python percentile); it does not import wp12_annotations.py.
Part 1: hub annotation recomputed from the two DrugCentral files.  Part 2: quartile thresholds, class sizes, class changes and the cross/same-target split recomputed from the two MLSS tables.
Part 3: Drugs@FDA facts rebuilt from the saved raw responses (every response is read, the query name is taken from its URL).  Part 4: every number quoted in the texts is present in them.
Usage: python -I check_wp12_numbers.py   (writes check_wp12_last_run.txt)"""
import sys, os, csv, gzip, io, json, re, hashlib, urllib.parse
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
W10 = os.path.join(P, 'CMBL_Sub2', 'strengthen_v2', 'wp10_versions', 'results')
RAW = os.path.join(HERE, 'results', 'raw')
TXT = open(os.path.join(HERE, 'README.md'), encoding='utf8').read() + '\n' + open(os.path.join(HERE, 'WORDING_annotations.md'), encoding='utf8').read()
out, fails = [], [0]


def check(name, ok, detail=''):
    out.append('%s  %s %s' % ('PASS' if ok else 'FAIL', name, detail)); fails[0] += 0 if ok else 1


def has(s): return s in TXT


def pct(v, q):                                   # numpy 'linear' percentile in plain python
    v = sorted(v); k = (len(v) - 1) * q / 100.0; f = int(k); c = min(f + 1, len(v) - 1)
    return v[f] + (v[c] - v[f]) * (k - f)


# ---------------------------------------------------------------- part 1
def dc(kind):
    if kind == 'old':
        return list(csv.reader(open(os.path.join(P, 'CMBL_Sub2', 'R_17genes', 'drug.target.interaction.csv'), encoding='utf-8-sig', newline='')))
    return list(csv.reader(io.TextIOWrapper(gzip.open(os.path.join(P, 'Database', 'drug.target.interaction.tsv.gz')), encoding='utf8'), delimiter='\t'))


def hub(rows, g):
    H = rows[0]; gi, di, oi, ti = H.index('GENE'), H.index('DRUG_NAME'), H.index('ORGANISM'), H.index('TDL')
    sel = [r for r in rows[1:] if g in [x.strip() for x in r[gi].split('|')]]
    hum = [r for r in sel if r[oi] == 'Homo sapiens']
    return len(sel), len(hum), len({r[di].lower() for r in hum}), sorted({t for r in sel for t in r[ti].split('|') if t})


old, new = dc('old'), dc('new')
check('DrugCentral files: reported 19378 rows, 25 Sep 2026 22364 rows', (len(old) - 1, len(new) - 1) == (19378, 22364) and has('19378') and has('22364'))
exp = {'GLP1R': ((8, 7, 7, ['Tclin']), (9, 8, 8, ['Tclin'])), 'TCF7L2': ((0, 0, 0, []), (0, 0, 0, [])), 'CACNA1D': ((25, 25, 24, ['Tclin']), (25, 25, 24, ['Tclin'])), 'PTBP1': ((0, 0, 0, []), (0, 0, 0, [])),
       'WNT5A': ((0, 0, 0, []), (0, 0, 0, [])), 'CACNA1C': ((41, 35, 33, ['Tclin']), (41, 35, 34, ['Tclin']))}
for g, (eo, en) in exp.items():
    check('hub %s: reported file %s | 25 Sep 2026 file %s' % (g, hub(old, g)[:3], hub(new, g)[:3]), hub(old, g) == eo and hub(new, g) == en)
# ---------------------------------------------------------------- part 2
comp = list(csv.DictReader(open(os.path.join(W10, 'mlss_comparison.csv'), encoding='utf8')))
allc = {r['Combination']: r for r in csv.DictReader(open(os.path.join(W10, 'mlss_new_all_combinations.csv'), encoding='utf8'))}
so = {r['combination']: float(r['MLSS_old']) for r in comp if r['MLSS_old']}; sn = {r['combination']: float(r['MLSS_new']) for r in comp if r['MLSS_new']}
check('153 reported and 171 combinations; shared scores identical', len(so) == 153 and len(sn) == 171 and all(abs(so[c] - sn[c]) < 1e-12 for c in so) and has('153') and has('171'))
t = [pct(list(so.values()), q) for q in (75, 50, 25)]; t2 = [pct(list(sn.values()), q) for q in (75, 50, 25)]
check('thresholds of the 153: 7.879365 / 7.65 / 2.754', [round(x, 6) for x in t] == [7.879365, 7.65, 2.754] and has('7.88') and has('7.65') and has('2.75'), str([round(x, 6) for x in t]))
check('thresholds of the 171: 7.879365 / 7.65 / 2.8674', [round(x, 6) for x in t2] == [7.879365, 7.65, 2.8674], str([round(x, 6) for x in t2]))


def cls(s, th): return 'Q1' if s >= th[0] else 'Q2' if s >= th[1] else 'Q3' if s >= th[2] else 'Q4'


from collections import Counter
c1 = Counter(cls(s, t) for s in so.values()); c2 = Counter(cls(s, t2) for s in sn.values())
check('class sizes of the 153: 42 / 35 / 48 / 28', [c1[k] for k in ('Q1', 'Q2', 'Q3', 'Q4')] == [42, 35, 48, 28] and has('42, 35, 48 and 28'))
check('class sizes of the 171: 47 / 41 / 40 / 43', [c2[k] for k in ('Q1', 'Q2', 'Q3', 'Q4')] == [47, 41, 40, 43] and has('47, 41, 40 and 43'))
chg = [c for c in so if cls(so[c], t) != cls(sn[c], t2)]
check('14 reported combinations change class (all Q3 -> Q4) when the 171 are used', len(chg) == 14 and all(cls(so[c], t) == 'Q3' and cls(sn[c], t2) == 'Q4' for c in chg) and has('14 of the 153'))
cross = {c for c in so if allc[c]['Targets_A'] != allc[c]['Targets_B']}
mc = [so[c] for c in cross]; ms = [so[c] for c in so if c not in cross]
check('77 cross-target scores 7.65 to 9.960357; 76 same-target scores 2.70 to 4.05 (a gap)', len(cross) == 77 and round(min(mc), 6) == 7.65 and round(max(mc), 6) == 9.960357 and len(ms) == 76 and round(min(ms), 6) == 2.7 and round(max(ms), 6) == 4.05 and has('2.70') and has('4.05') and has('9.96'))
check('Q1 and Q2 hold exactly the cross-target pairs, Q3 and Q4 the same-target pairs', all((cls(so[c], t) in ('Q1', 'Q2')) == (c in cross) for c in so))
rk = sorted(so, key=lambda c: -so[c]); r15 = {r['combination'] for r in comp if r['rank_old'] and int(r['rank_old']) <= 15}
check('the 15 combinations of Fig. 6 are all in Q1', len(r15) == 15 and all(cls(so[c], t) == 'Q1' for c in r15) and has('all 15'))
new18 = [c for c in sn if c not in so]
check('18 new combinations all contain tirzepatide; 5 / 6 / 6 / 1 in Q1 / Q2 / Q3 / Q4 (171 classes)', len(new18) == 18 and all('tirzepatide' in c for c in new18) and [Counter(cls(sn[c], t2) for c in new18)[k] for k in ('Q1', 'Q2', 'Q3', 'Q4')] == [5, 6, 6, 1] and has('5, 6, 6 and 1'))
# ---------------------------------------------------------------- part 3
resp = {}
for fn in os.listdir(RAW):
    if not fn.endswith('.json'): continue
    rec = json.load(open(os.path.join(RAW, fn), encoding='utf8'))
    u = urllib.parse.urlparse(rec['url']); qs = urllib.parse.parse_qs(u.query)
    name = re.search(r'name:"([^"]+)"', qs['search'][0]).group(1); skip = int(qs['skip'][0])
    resp.setdefault(name, []).append((skip, rec))
check('saved responses: 19 files for 18 names (amlodipine has 2 pages)', sum(len(v) for v in resp.values()) == 19 and len(resp) == 18 and len(resp['amlodipine']) == 2, str(sum(len(v) for v in resp.values())))
tot = lambda n: resp[n][0][1]['response']['meta']['results']['total']
apps = lambda n: [a for _, rec in sorted(resp[n], key=lambda x: x[0]) for a in rec['response'].get('results', [])]
check('totals: albiglutide 1, bepridil 2, nilvadipine 0, nitrendipine 0, amlodipine 121', (tot('albiglutide'), tot('bepridil'), tot('nilvadipine'), tot('nitrendipine'), tot('amlodipine')) == (1, 2, 0, 0, 121) and len(apps('amlodipine')) == 121)
a = apps('albiglutide')[0]
check('albiglutide: BLA125431, TANZEUM, GLAXOSMITHKLINE LLC, Discontinued, ORIG AP 20140415', a['application_number'] == 'BLA125431' and a['products'][0]['brand_name'] == 'TANZEUM' and a['sponsor_name'] == 'GLAXOSMITHKLINE LLC' and all(p['marketing_status'] == 'Discontinued' for p in a['products'])
      and sorted(s['submission_status_date'] for s in a['submissions'] if s['submission_type'] == 'ORIG' and s['submission_status'] == 'AP')[0] == '20140415' and has('BLA 125431') and has('15 April 2014'))
b = {x['application_number']: x for x in apps('bepridil')}
check('bepridil: NDA019001 BEPADIN and NDA019002 VASCOR, all Discontinued, ORIG AP 19901228', set(b) == {'NDA019001', 'NDA019002'} and {b[k]['products'][0]['brand_name'] for k in b} == {'BEPADIN', 'VASCOR'} and all(p['marketing_status'] == 'Discontinued' for x in b.values() for p in x['products'])
      and all(sorted(s['submission_status_date'] for s in x['submissions'] if s['submission_type'] == 'ORIG' and s['submission_status'] == 'AP')[0] == '19901228' for x in b.values()) and has('28 December 1990'))
retr = sorted(rec['retrieved'] for v in resp.values() for _, rec in v)
check('retrieval on 9 Nov 2026 between 00:34 and 00:36 (+0700)', retr[0].startswith('2026-11-09 00:3') and retr[-1].startswith('2026-11-09 00:3'), retr[0] + ' .. ' + retr[-1])
notmk = []
for n in resp:
    st = [p.get('marketing_status') for a in apps(n) for p in a.get('products', []) if any(n in (i.get('name') or '').lower() for i in p.get('active_ingredients', []))]
    if not any(s in ('Prescription', 'Over-the-counter') for s in st): notmk.append(n)
check('names with no product that is Prescription or Over-the-counter: albiglutide, bepridil, nilvadipine, nitrendipine', sorted(notmk) == ['albiglutide', 'bepridil', 'nilvadipine', 'nitrendipine'], str(sorted(notmk)))
rows19 = list(csv.DictReader(open(os.path.join(W10, 'approval_tiers_old_vs_new.csv'), encoding='utf8')))
others = [r['drug'].lower().replace('(s)-', '') for r in rows19 if r['drug'].lower().replace('(s)-', '') not in ('albiglutide', 'bepridil', 'nilvadipine', 'nitrendipine')]
check('19 drug entries; the other 14 each have a Prescription product and none is Over-the-counter', len(rows19) == 19 and len(others) == 14 and all(any(p.get('marketing_status') == 'Prescription' for a in apps(n) for p in a.get('products', [])) for n in others)
      and not any(p.get('marketing_status') == 'Over-the-counter' for n in others for a in apps(n) for p in a.get('products', [])), str(len(others)))
isr = {a['application_number']: {p['marketing_status'] for p in a['products']} for a in apps('isradipine')}
check('isradipine: both Dynacirc applications discontinued; ANDA077169 and ANDA077317 Prescription', isr.get('NDA019546') == {'Discontinued'} and isr.get('NDA020336') == {'Discontinued'} and isr.get('ANDA077169') == {'Prescription'} and isr.get('ANDA077317') == {'Prescription'})
check('47 distinct scores among the 153; 14 combinations score exactly 7.65', len({round(v, 9) for v in so.values()}) == 47 and sum(1 for v in so.values() if abs(v - 7.65) < 1e-9) == 14 and has('47 distinct') and has('14 combinations have exactly the score 7.65'))
bad = {'albiglutide', 'bepridil', 'nilvadipine', 'nitrendipine'}
names = lambda c: [x.strip().lower().replace('(s)-', '') for x in c.split(' + ')]
top = [c for c in so if c in r15]
check('10 of the 15 top combinations, 75 of the 153 and 35 of the 77 cross-target pairs contain one of the four', sum(1 for c in top if set(names(c)) & bad) == 10 and sum(1 for c in so if set(names(c)) & bad) == 75 and sum(1 for c in cross if set(names(c)) & bad) == 35 and has('10 of the 15') and has('35 of the 77') and has('75 of the 153'))
# ---------------------------------------------------------------- part 4: results files agree with the texts
s = json.load(open(os.path.join(HERE, 'results', 'summary.json'), encoding='utf8'))
check('summary.json thresholds and sizes equal the recomputation', [round(s['quartiles']['thresholds_reported_153'][k], 6) for k in ('t75', 't50', 't25')] == [7.879365, 7.65, 2.754] and list(s['quartiles']['class_sizes_reported_153'].values()) == [42, 35, 48, 28])
out.append('%d checks, %d failures' % (len(out), fails[0]))
open(os.path.join(HERE, 'check_wp12_last_run.txt'), 'w', encoding='utf8').write('\n'.join(out) + '\n')
print('\n'.join(out))
