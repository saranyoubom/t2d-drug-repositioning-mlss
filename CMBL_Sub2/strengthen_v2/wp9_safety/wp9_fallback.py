"""wp9_fallback.py - SECONDARY, POST HOC analysis of the WP9 safety screen (pre-specification, section 9, entry of 7 Nov 15:22).
Drugs with fewer than 100 reports through patient.drug.openfda.substance_name (discontinued products carry no harmonized openfda fields) are searched again through the
free-text field patient.drug.medicinalproduct with names fixed in advance (NAMES below). Same statistics, same thresholds and same flags as the primary analysis; the partner drug of a pair keeps its
primary counts. Everything this script writes is labelled secondary and post hoc and is never merged into the primary counts.
Usage (python -I): python -I wp9_fallback.py retrieve   (network, cached, at most 120 new requests per run)
                   python -I wp9_fallback.py analyze    (offline; writes results/fallback_signals_single.csv, fallback_pairs.csv, fallback_summary.json)"""
import os, sys, csv, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fda_fetch as ff
import wp9_screen as W

TEXT = 'patient.drug.medicinalproduct'
NAMES = {'albiglutide': ['TANZEUM', 'ALBIGLUTIDE'], 'bepridil': ['VASCOR', 'BEPRIDIL'], 'nilvadipine': ['NILVADIPINE'], 'nitrendipine': ['NITRENDIPINE']}
MAX_NEW = 120
RES = W.RES


def tq(t):
    return '%s:(%s)' % (TEXT, ' OR '.join('"%s"' % n for n in NAMES[t]))


def partners(cross, t):
    out = []
    for r in cross:
        a, b = W.term(r['A']), W.term(r['B'])
        if t in (a, b):
            out.append((r, b if a == t else a))
    return out


def low_terms(terms):
    return [t for t in terms if ff.total(W.dq(t)) < W.MIN_DRUG]


def budget():
    if ff._state['requests'] >= MAX_NEW:
        print('request budget of this run reached (%d new requests); run again later' % ff._state['requests'], flush=True)
        sys.exit(2)


def retrieve():
    cls, rows, cross = W.load_mlss()
    terms = sorted({W.term(n) for n in cls})
    low = low_terms(terms)
    print('drugs below %d reports through substance_name:' % W.MIN_DRUG, low, flush=True)
    for t in low:
        if t not in NAMES:
            print('no fallback name was fixed in advance for', t, '(not searched)')
            continue
        n = ff.total(tq(t))
        print('fallback %-12s %8d reports (free text %s)' % (t, n, NAMES[t]), flush=True)
        budget()
        if n < W.MIN_DRUG:
            continue
        for c in W.CODES:
            ff.total(tq(t) + ' AND ' + W.gq(c))
            budget()
        for r, p in sorted({(x['rank'], p): (x, p) for x, p in partners(cross, t)}.values(), key=lambda z: z[0]['rank']):
            nab = ff.total(tq(t) + ' AND ' + W.dq(p))
            print('  pair %-12s + %-14s %7d co-reports' % (t, p, nab), flush=True)
            budget()
            if nab >= W.MIN_PAIR:
                for c in W.CODES:
                    ff.total(tq(t) + ' AND ' + W.dq(p) + ' AND ' + W.gq(c))
                    budget()
    print('fallback retrieval complete: %d new requests, %d from the saved files' % (ff._state['requests'], ff._state['cached']))


def analyze():
    ff.set_offline(True)
    cls, rows, cross = W.load_mlss()
    terms = sorted({W.term(n) for n in cls})
    N = ff.total(None)
    E = {c: ff.total(W.gq(c)) for c in W.CODES}
    nd = {t: ff.total(W.dq(t)) for t in terms}
    ev = {t: nd[t] >= W.MIN_DRUG for t in terms}
    a_dg = {(t, c): ff.total(W.dq(t) + ' AND ' + W.gq(c)) for t in terms if ev[t] for c in W.CODES}
    sig_primary = {}
    for (t, c), a in a_dg.items():
        o, lo, hi = W.ror(a, nd[t], E[c], N)
        sig_primary[(t, c)] = a >= W.MIN_A and lo > 1
    single, pairs_out, summary = [], [], {'note': 'SECONDARY, POST HOC, free-text names; not merged into the primary counts', 'drugs': {}}
    for t in low_terms(terms):
        if t not in NAMES:
            summary['drugs'][t] = {'searched': False}
            continue
        n_t = ff.total(tq(t))
        info = {'searched': True, 'names': NAMES[t], 'reports': n_t, 'evaluable': n_t >= W.MIN_DRUG}
        summary['drugs'][t] = info
        if n_t < W.MIN_DRUG:
            continue
        a_t = {c: ff.total(tq(t) + ' AND ' + W.gq(c)) for c in W.CODES}
        sig_t = {}
        for c in W.CODES:
            o, lo, hi = W.ror(a_t[c], n_t, E[c], N)
            s = a_t[c] >= W.MIN_A and lo > 1
            sig_t[c] = s
            single.append([t, cls[[n for n in cls if W.term(n) == t][0]], c, a_t[c], n_t, '%.5f' % (a_t[c] / n_t), E[c], '%.3f' % o, '%.3f' % lo, '%.3f' % hi, int(s)])
        info['signals'] = [c for c in W.CODES if sig_t[c]]
        for r, p in partners(cross, t):
            nab = ff.total(tq(t) + ' AND ' + W.dq(p))
            pev = nab >= W.MIN_PAIR
            cev = ev[p]
            sP = [c for c in W.CODES if sig_primary.get((p, c))]
            flags = []
            if pev and cev:
                for c in W.CODES:
                    x = ff.total(tq(t) + ' AND ' + W.dq(p) + ' AND ' + W.gq(c))
                    lo, hi = W.cp(x, nab)
                    pA = (a_t[c] - x) / (n_t - nab) if n_t > nab else float('nan')
                    pB = (a_dg[(p, c)] - x) / (nd[p] - nab) if nd[p] > nab else float('nan')
                    if x >= W.MIN_A and lo > max(pA, pB):
                        flags.append(c)
            status = 'evaluable' if (pev and cev) else ('partner drug below %d reports' % W.MIN_DRUG if pev else 'fewer than %d co-reports' % W.MIN_PAIR)
            pairs_out.append([r['rank'], r['pair'], '%.4f' % r['mlss'], t, p, nab, status, ' '.join(c for c in W.CODES if sig_t[c]) or '-', (' '.join(sP) or '-') if cev else 'n/e', (' '.join(flags) or '-') if (pev and cev) else 'n/a'])
    pairs_out.sort(key=lambda z: z[0])
    W.write_csv('fallback_signals_single.csv', ['term', 'target_gene', 'group', 'reports_in_group', 'reports_of_drug', 'proportion', 'reports_in_group_all_drugs', 'ROR', 'ROR_lo95', 'ROR_hi95', 'signal'], single)
    W.write_csv('fallback_pairs.csv', ['rank', 'pair', 'MLSS', 'fallback_term', 'partner_term', 'reports_listing_both', 'status', 'signals_fallback_drug', 'signals_partner_drug', 'co_report_excess'], pairs_out)
    summary['pairs'] = len(pairs_out)
    summary['pairs_evaluable'] = sum(1 for p in pairs_out if p[6] == 'evaluable')
    summary['pairs_evaluable_top15'] = sum(1 for p in pairs_out if p[0] <= 15 and p[6] == 'evaluable')
    summary['pairs_top15'] = sum(1 for p in pairs_out if p[0] <= 15)
    summary['co_report_flags'] = sum(0 if p[9] in ('-', 'n/a') else len(p[9].split()) for p in pairs_out)
    with open(os.path.join(RES, 'fallback_summary.json'), 'w', encoding='utf8') as f:
        json.dump(summary, f, indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != 'drugs'}), {t: (v.get('reports'), v.get('signals')) for t, v in summary['drugs'].items()})


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    if cmd == 'retrieve':
        retrieve()
    elif cmd == 'analyze':
        analyze()
    else:
        print(__doc__)
