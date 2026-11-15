"""wp9_screen.py - WP9: pharmacovigilance safety screen (an annotation, NOT a prediction) of the MLSS drugs and cross-target drug pairs.
Pre-specification: PRESPEC_WP9_safety_screen.md (fingerprinted before any report count of an MLSS drug was retrieved); every rule below is from it.
Usage (always python -I):
  python -I wp9_screen.py retrieve   needs the network; cached and resumable; at most 850 new requests per run (openFDA allows 1000 per day without a key)
  python -I wp9_screen.py analyze    reads only results/raw/*.json (offline); writes results/drug_counts.csv, signals_single.csv, pairs_all77.csv,
                                     pairs_long.csv, pairs_top15.csv and summary.json
Inputs (read only): ../wp3_mlss/direction_table.csv, ../wp3_mlss/mlss_scenario_scores.csv."""
import os, sys, csv, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)              # python -I does not add the script folder itself
from scipy import stats
import fda_fetch as ff

V2 = os.path.dirname(HERE)
RES = os.path.join(HERE, 'results')
D = 'patient.drug.openfda.substance_name'
R = 'patient.reaction.reactionmeddrapt.exact'
MIN_DRUG, MIN_PAIR, MIN_A, MAX_NEW = 100, 20, 3, 850
Z = 1.959963984540054
GROUPS = [
    ('HYPO', 'hypoglycemia', ['HYPOGLYCAEMIA', 'HYPOGLYCAEMIC COMA', 'HYPOGLYCAEMIC UNCONSCIOUSNESS', 'BLOOD GLUCOSE DECREASED']),
    ('HYPER', 'hyperglycemia, loss of glycemic control', ['HYPERGLYCAEMIA', 'BLOOD GLUCOSE INCREASED', 'DIABETES MELLITUS INADEQUATE CONTROL', 'DIABETIC KETOACIDOSIS', 'GLYCOSYLATED HAEMOGLOBIN INCREASED']),
    ('PANC', 'pancreatic events', ['PANCREATITIS', 'PANCREATITIS ACUTE', 'PANCREATITIS CHRONIC', 'PANCREATITIS NECROTISING', 'LIPASE INCREASED', 'AMYLASE INCREASED']),
    ('GI', 'gastrointestinal intolerance', ['NAUSEA', 'VOMITING', 'DIARRHOEA', 'ABDOMINAL PAIN', 'ABDOMINAL PAIN UPPER', 'DECREASED APPETITE', 'CONSTIPATION']),
    ('BILI', 'biliary events', ['CHOLELITHIASIS', 'CHOLECYSTITIS', 'CHOLECYSTITIS ACUTE', 'BILIARY COLIC']),
    ('HYPOT', 'hypotension, syncope', ['HYPOTENSION', 'ORTHOSTATIC HYPOTENSION', 'SYNCOPE', 'PRESYNCOPE', 'BLOOD PRESSURE DECREASED']),
    ('OEDEMA', 'edema', ['OEDEMA PERIPHERAL', 'OEDEMA', 'PERIPHERAL SWELLING']),
    ('CARD', 'conduction and rhythm', ['BRADYCARDIA', 'ATRIOVENTRICULAR BLOCK', 'ATRIOVENTRICULAR BLOCK COMPLETE', 'ELECTROCARDIOGRAM QT PROLONGED', 'TORSADE DE POINTES', 'VENTRICULAR TACHYCARDIA']),
]
CODES = [g[0] for g in GROUPS]
PT = {g[0]: g[2] for g in GROUPS}


def term(name):
    return 'nitrendipine' if name == '(S)-nitrendipine' else name.lower()


def dq(t):
    return '%s:"%s"' % (D, t)


def gq(c):
    return '%s:(%s)' % (R, ' OR '.join('"%s"' % t for t in PT[c]))


def pq(a, b):
    return dq(a) + ' AND ' + dq(b)


def load_mlss():
    cls = {}
    with open(os.path.join(V2, 'wp3_mlss', 'direction_table.csv'), encoding='utf8', newline='') as f:
        for r in csv.DictReader(f):
            if r['used_in_mlss'] == 'True':
                cls[r['drug']] = r['target_gene']
    rows = []
    with open(os.path.join(V2, 'wp3_mlss', 'mlss_scenario_scores.csv'), encoding='utf8', newline='') as f:
        for r in csv.DictReader(f):
            a, b = r['Combination'].split(' + ')
            rows.append({'rank': int(float(r['Rank'])), 'pair': r['Combination'], 'mlss': float(r['reported_MLSS']), 'A': a, 'B': b, 'cA': cls[a], 'cB': cls[b]})
    cross = sorted([r for r in rows if r['cA'] != r['cB']], key=lambda r: r['rank'])
    assert len(cls) == 18 and len(rows) == 153 and len(cross) == 77, (len(cls), len(rows), len(cross))
    assert [r['rank'] for r in cross] == list(range(1, 78))
    return cls, rows, cross


def pair_keys(cross):
    seen, keys = set(), []
    for r in cross:
        k = tuple(sorted((term(r['A']), term(r['B']))))
        if k not in seen:
            seen.add(k)
            keys.append(k)
    return keys


def budget():
    if ff._state['requests'] >= MAX_NEW:
        print('request budget of this run reached (%d new requests); run again later (the daily limit without a key is 1000)' % ff._state['requests'], flush=True)
        sys.exit(2)


def retrieve():
    cls, rows, cross = load_mlss()
    terms = sorted({term(n) for n in cls})
    keys = pair_keys(cross)
    assert len(terms) == 17 and len(keys) == 70, (len(terms), len(keys))
    print('database: all reports', ff.total(None), ', last updated', ff.last_updated(), flush=True)
    for c in CODES:                       # spelling check of every preferred term, independent of the drugs
        for t in PT[c]:
            n = ff.total('%s:"%s"' % (R, t))
            print('  PT %-7s %-38s %9d' % (c, t, n), flush=True)
            if n == 0:
                print('PT without any report: fix the spelling BEFORE any drug count is retrieved and log it in section 9 of the pre-specification')
                sys.exit(3)
            budget()
    for c in CODES:
        print('group', c, ff.total(gq(c)), flush=True)
        budget()
    for t in terms:
        n = ff.total(dq(t))
        print('drug %-14s %9d reports' % (t, n), flush=True)
        if n >= MIN_DRUG:
            for c in CODES:
                ff.total(dq(t) + ' AND ' + gq(c))
                budget()
        budget()
    for i, (a, b) in enumerate(keys, 1):
        n = ff.total(pq(a, b))
        print('pair %2d/70 %-14s + %-14s %8d reports%s' % (i, a, b, n, '' if n >= MIN_PAIR else '  (not evaluable)'), flush=True)
        if n >= MIN_PAIR:
            for c in CODES:
                ff.total(pq(a, b) + ' AND ' + gq(c))
                budget()
        budget()
    print('retrieval complete: %d new requests in this run, %d answered from the saved files' % (ff._state['requests'], ff._state['cached']))


def ror(a, n_d, E, N):
    cells = [a, n_d - a, E - a, N - n_d - E + a]
    if min(cells) <= 0:
        cells = [x + 0.5 for x in cells]
    a_, b_, c_, d_ = cells
    lor = math.log(a_ * d_ / (b_ * c_))
    se = math.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
    return math.exp(lor), math.exp(lor - Z * se), math.exp(lor + Z * se)


def cp(a, n):
    lo = 0.0 if a == 0 else float(stats.beta.ppf(0.025, a, n - a + 1))
    hi = 1.0 if a == n else float(stats.beta.ppf(0.975, a + 1, n - a))
    return lo, hi


def write_csv(name, header, data):
    with open(os.path.join(RES, name), 'w', encoding='utf8', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(data)


def analyze():
    ff.set_offline(True)
    cls, rows, cross = load_mlss()
    N = ff.total(None)
    LU = ff.last_updated()
    E = {c: ff.total(gq(c)) for c in CODES}
    entries = {}
    for n in cls:
        entries.setdefault(term(n), []).append(n)
    terms = sorted(entries)
    tclass = {t: cls[entries[t][0]] for t in terms}
    nd = {t: ff.total(dq(t)) for t in terms}
    ev = {t: nd[t] >= MIN_DRUG for t in terms}
    a_dg = {(t, c): ff.total(dq(t) + ' AND ' + gq(c)) for t in terms if ev[t] for c in CODES}
    warn = []
    # single-drug table
    sig, single_rows = {}, []
    for t in terms:
        if not ev[t]:
            continue
        for c in CODES:
            a = a_dg[(t, c)]
            if a > min(nd[t], E[c]):
                warn.append('a_dg > min(n_d, E_g): %s %s' % (t, c))
            o, lo, hi = ror(a, nd[t], E[c], N)
            s = a >= MIN_A and lo > 1
            sig[(t, c)] = s
            single_rows.append([t, tclass[t], c, a, nd[t], '%.5f' % (a / nd[t]), E[c], '%.3f' % o, '%.3f' % lo, '%.3f' % hi, int(s)])
    write_csv('signals_single.csv', ['term', 'target_gene', 'group', 'reports_in_group', 'reports_of_drug', 'proportion', 'reports_in_group_all_drugs', 'ROR', 'ROR_lo95', 'ROR_hi95', 'signal'], single_rows)
    write_csv('drug_counts.csv', ['term', 'mlss_entries', 'target_gene', 'reports', 'evaluable'], [[t, ';'.join(entries[t]), tclass[t], nd[t], int(ev[t])] for t in terms])
    # pairs
    nab, aab = {}, {}
    for k in pair_keys(cross):
        nab[k] = ff.total(pq(*k))
        if nab[k] >= MIN_PAIR:
            for c in CODES:
                aab[(k, c)] = ff.total(pq(*k) + ' AND ' + gq(c))
    pair_rows, long_rows = [], []
    for r in cross:
        tA, tB = term(r['A']), term(r['B'])
        k = tuple(sorted((tA, tB)))
        n_ab = nab[k]
        pev = n_ab >= MIN_PAIR
        cev = ev[tA] and ev[tB]
        sA = [c for c in CODES if sig.get((tA, c))]
        sB = [c for c in CODES if sig.get((tB, c))]
        flags = []
        if n_ab > min(nd[tA], nd[tB]):
            warn.append('n_AB > min(n_A, n_B): %s' % (k,))
        if pev and cev:
            for c in CODES:
                a = aab[(k, c)]
                if a > min(a_dg[(tA, c)], a_dg[(tB, c)], n_ab):
                    warn.append('a_AB > min(a_A, a_B, n_AB): %s %s' % (k, c))
                p = a / n_ab
                lo, hi = cp(a, n_ab)
                pA = (a_dg[(tA, c)] - a) / (nd[tA] - n_ab) if nd[tA] > n_ab else float('nan')
                pB = (a_dg[(tB, c)] - a) / (nd[tB] - n_ab) if nd[tB] > n_ab else float('nan')
                fl = a >= MIN_A and lo > max(pA, pB)
                if fl:
                    flags.append(c)
                long_rows.append([r['rank'], r['pair'], c, a, n_ab, '%.5f' % p, '%.5f' % lo, '%.5f' % hi, '%.5f' % pA, '%.5f' % pB, int(fl)])
        status = 'evaluable' if (pev and cev) else ('component drug below %d reports' % MIN_DRUG if pev else 'fewer than %d co-reports' % MIN_PAIR)
        pair_rows.append([r['rank'], r['pair'], '%.4f' % r['mlss'], tA, tB, n_ab, status, (' '.join(sA) or '-') if ev[tA] else 'n/e', (' '.join(sB) or '-') if ev[tB] else 'n/e', (' '.join(flags) or '-') if (pev and cev) else 'n/a'])
    hdr = ['rank', 'pair', 'MLSS', 'term_A', 'term_B', 'reports_listing_both', 'status', 'signals_drug_A', 'signals_drug_B', 'co_report_excess']  # '-' = none, 'n/e' = drug not evaluable, 'n/a' = pair not evaluable
    write_csv('pairs_all77.csv', hdr, pair_rows)
    write_csv('pairs_top15.csv', hdr, [p for p in pair_rows if p[0] <= 15])
    write_csv('pairs_long.csv', ['rank', 'pair', 'group', 'reports_both_in_group', 'reports_both', 'p_both', 'p_both_lo95', 'p_both_hi95', 'p_A_without_B', 'p_B_without_A', 'co_report_excess'], long_rows)
    # summary
    n_sig = sum(1 for v in sig.values() if v)
    per_group = {}
    for c in CODES:
        drugs_c = sorted(t for t in terms if sig.get((t, c)))
        top = sum(1 for r in cross if r['rank'] <= 15 and (sig.get((term(r['A']), c)) or sig.get((term(r['B']), c))))
        allp = sum(1 for r in cross if (sig.get((term(r['A']), c)) or sig.get((term(r['B']), c))))
        per_group[c] = {'E_g': E[c], 'signalled_drugs': drugs_c, 'pairs_top15_with_signalled_drug': top, 'pairs_all77_with_signalled_drug': allp}
    flagged = [(p[0], p[1], p[9]) for p in pair_rows if p[9] not in ('-', 'n/a')]
    summary = {'last_updated': LU, 'N_reports': N, 'terms': len(terms), 'evaluable_drugs': sum(ev.values()), 'not_evaluable_drugs': [t for t in terms if not ev[t]],
               'drug_reports': nd, 'single_signals_total': n_sig, 'single_tests_total': len(single_rows), 'per_group': per_group,
               'pairs_total': len(cross), 'pairs_with_20_or_more_co_reports': sum(1 for p in pair_rows if not p[6].startswith('fewer')),
               'pairs_evaluable': sum(1 for p in pair_rows if p[6] == 'evaluable'), 'pairs_evaluable_top15': sum(1 for p in pair_rows if p[0] <= 15 and p[6] == 'evaluable'),
               'co_report_flags': sum(len(f[2].split()) for f in flagged), 'co_report_flagged_pair_count': len(flagged), 'co_report_flagged_pairs': flagged, 'pair_group_tests': len(long_rows), 'warnings': warn}
    with open(os.path.join(RES, 'summary.json'), 'w', encoding='utf8') as f:
        json.dump(summary, f, indent=1)
    print('N reports %d (last updated %s); drugs evaluable %d of %d; single-drug signals %d of %d tests; pairs evaluable %d of %d; co-report flags %d (in %d pairs) of %d pair-group tests; warnings %d'
          % (N, LU, summary['evaluable_drugs'], len(terms), n_sig, len(single_rows), summary['pairs_evaluable'], len(cross), summary['co_report_flags'], len(flagged), len(long_rows), len(warn)))
    for w in warn:
        print('WARNING', w)


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    if cmd == 'retrieve':
        retrieve()
    elif cmd == 'analyze':
        analyze()
    else:
        print(__doc__)
