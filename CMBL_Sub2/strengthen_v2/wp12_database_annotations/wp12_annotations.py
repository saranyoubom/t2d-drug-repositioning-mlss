"""wp12_annotations.py - WP12 (9 Nov 2026): annotations from the updated databases, as pre-specified in PRESPEC_WP12_annotations.md (fingerprinted before the retrieval).
  part 1  hub-gene target annotation from the DrugCentral file of 25 Sep 2026 (record counts, drugs, TDL)
  part 2  score-quartile classes of the MLSS combinations (153 reported, 171 of the 25 Sep 2026 release)
  part 3  regulatory / market status of the 19 MLSS drugs: DrugCentral approval-list membership and openFDA Drugs@FDA (retrieved on the day of the run)
Usage (python -I):  python -I wp12_annotations.py retrieve   (network: openFDA drugsfda, cached in results/raw/, safe to repeat)
                    python -I wp12_annotations.py analyze    (offline: reads saved files only; writes results/*.csv and results/summary.json)
Read only with respect to every other folder of the project."""
import os, sys, csv, gzip, io, json, time, hashlib, urllib.request, urllib.parse, urllib.error
import numpy as np
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
W10 = os.path.join(P, 'CMBL_Sub2', 'strengthen_v2', 'wp10_versions', 'results')
RAW = os.path.join(HERE, 'results', 'raw')
RES = os.path.join(HERE, 'results')
os.makedirs(RAW, exist_ok=True)
BASE = 'https://api.fda.gov/drug/drugsfda.json'
HUBS = ['GLP1R', 'TCF7L2', 'CACNA1D', 'PTBP1']
CONTEXT = ['WNT5A', 'CACNA1C']
MIN_GAP = 0.5
_last = [0.0]


def url_for(name, skip=0, limit=100):
    return BASE + '?search=' + urllib.parse.quote_plus('products.active_ingredients.name:"%s"' % name) + '&limit=%d&skip=%d' % (limit, skip)


def fetch(url, tries=8):
    path = os.path.join(RAW, hashlib.sha1(url.encode()).hexdigest() + '.json')
    if os.path.exists(path):
        return json.load(open(path, encoding='utf8'))['response']
    delay = 3.0
    for k in range(tries):
        gap = time.time() - _last[0]
        if gap < MIN_GAP: time.sleep(MIN_GAP - gap)
        _last[0] = time.time()
        err = None
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'research-screen/1.0'})
            with urllib.request.urlopen(req, timeout=120) as r: resp = json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                resp = {'meta': {'results': {'total': 0, 'skip': 0, 'limit': 100}}, 'results': [], 'note': 'HTTP 404: no match'}
            else:
                err = 'HTTP %s' % e.code
                if e.code not in (429, 500, 502, 503, 504): raise
        except Exception as e:
            err = '%s: %s' % (type(e).__name__, str(e)[:100])
        if err is None:
            tmp = path + '.tmp'
            json.dump({'url': url, 'retrieved': time.strftime('%Y-%m-%d %H:%M:%S %z'), 'response': resp}, open(tmp, 'w', encoding='utf8'))
            os.replace(tmp, path)
            return resp
        print('  retry %d/%d after %s' % (k + 1, tries, err), flush=True); time.sleep(delay)
    raise RuntimeError('failed: ' + url)


def drug_list():
    rows = list(csv.DictReader(open(os.path.join(W10, 'approval_tiers_old_vs_new.csv'), encoding='utf8')))
    return rows


def query_name(drug):
    d = drug.lower()
    return d[4:] if d.startswith('(s)-') else d


def listing(name):
    """all applications for an active-ingredient name (pages of 100, at most 6 pages)"""
    out, skip, total = [], 0, None
    while True:
        r = fetch(url_for(name, skip))
        total = r.get('meta', {}).get('results', {}).get('total', 0)
        out += r.get('results', [])
        skip += 100
        if skip >= total or skip >= 600: break
    return total, out


def retrieve():
    names = sorted({query_name(r['drug']) for r in drug_list()})
    print('queries for %d names' % len(names))
    for n in names:
        t, res = listing(n)
        print('  %-14s total applications %d (retrieved %d)' % (n, t, len(res)), flush=True)


# --------------------------------------------------------------------------- analysis
def read_dc(kind):
    """kind 'reported': Analysis_outputs/R_17genes/drug.target.interaction.csv (the file of the reported run); kind 'new': Database/drug.target.interaction.tsv.gz (dump of 25 Sep 2026)"""
    if kind == 'reported':
        return list(csv.reader(open(os.path.join(P, 'CMBL_Sub2', 'R_17genes', 'drug.target.interaction.csv'), encoding='utf-8-sig', newline='')))
    return list(csv.reader(io.TextIOWrapper(gzip.open(os.path.join(P, 'Database', 'drug.target.interaction.tsv.gz')), encoding='utf8'), delimiter='	'))


def annot(rows, g):
    h = rows[0]; gi, di, ti, oi, mi = h.index('GENE'), h.index('DRUG_NAME'), h.index('TDL'), h.index('ORGANISM'), h.index('MOA')
    sel = [r for r in rows[1:] if g in r[gi].replace(' ', '').split('|')]
    hum = [r for r in sel if r[oi] == 'Homo sapiens']
    tdl = sorted({x for r in sel for x in r[ti].split('|') if x})
    return {'records_all_species': len(sel), 'records_human': len(hum), 'distinct_drugs_human': len({r[di].lower() for r in hum}),
            'moa_records_human': sum(1 for r in hum if r[mi] == '1'), 'TDL_values': '|'.join(tdl) if tdl else ('no record in the file' if not sel else 'none stated')}


def part1():
    summ = json.load(open(os.path.join(W10, 'drugcentral_summary.json'), encoding='utf8'))
    old, new = read_dc('reported'), read_dc('new')
    out = []
    for g in HUBS + CONTEXT:
        a, b = annot(old, g), annot(new, g)
        r = {'gene': g, 'role': 'hub (degree >= 5)' if g in HUBS else 'context (degree 4)'}
        for k, v in a.items(): r['reported_file_' + k] = v
        for k, v in b.items(): r['file_25Sep2026_' + k] = v
        r['MLSS_drugs_reported_file'] = summ['drugs_per_hub_old'].get(g, '')
        r['MLSS_drugs_25Sep2026_file'] = summ['drugs_per_hub_new'].get(g, '')
        out.append(r)
    with open(os.path.join(RES, 'hub_annotation.csv'), 'w', encoding='utf8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    return out


def classes(scores):
    t75, t50, t25 = (float(np.percentile(scores, q)) for q in (75, 50, 25))
    def cl(s): return 'Q1' if s >= t75 else 'Q2' if s >= t50 else 'Q3' if s >= t25 else 'Q4'
    return {'t75': t75, 't50': t50, 't25': t25}, cl


def part2():
    comp = list(csv.DictReader(open(os.path.join(W10, 'mlss_comparison.csv'), encoding='utf8')))
    allc = {r['Combination']: r for r in csv.DictReader(open(os.path.join(W10, 'mlss_new_all_combinations.csv'), encoding='utf8'))}
    old = {r['combination']: float(r['MLSS_old']) for r in comp if r['MLSS_old'] != ''}
    new = {r['combination']: float(r['MLSS_new']) for r in comp if r['MLSS_new'] != ''}
    ro = {r['combination']: int(r['rank_old']) for r in comp if r['rank_old'] != ''}
    rn = {r['combination']: int(r['rank_new']) for r in comp if r['rank_new'] != ''}
    th_o, cl_o = classes(list(old.values())); th_n, cl_n = classes(list(new.values()))
    rows, size_o, size_n, changes = [], {}, {}, []
    for c in sorted(new, key=lambda k: rn[k]):
        a = allc[c]; cross = a['Targets_A'] != a['Targets_B']
        co = cl_o(old[c]) if c in old else ''; cn = cl_n(new[c])
        if co: size_o[co] = size_o.get(co, 0) + 1
        size_n[cn] = size_n.get(cn, 0) + 1
        if co and co != cn: changes.append((c, co, cn))
        rows.append({'combination': c, 'in_reported_153': int(c in old), 'MLSS': round(new[c], 6), 'rank_reported': ro.get(c, ''), 'rank_25Sep2026': rn[c], 'class_reported_153': co, 'class_171': cn,
                     'cross_target_pair': int(cross), 'targets': '%s + %s' % (a['Targets_A'], a['Targets_B']), 'reported_top15': int(c in ro and ro[c] <= 15)})
    with open(os.path.join(RES, 'mlss_quartile_classes.csv'), 'w', encoding='utf8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    top15 = [(r['combination'], r['rank_reported'], r['class_reported_153']) for r in rows if r['reported_top15']]
    return {'n_reported': len(old), 'n_171': len(new), 'thresholds_reported_153': th_o, 'thresholds_171': th_n, 'class_sizes_reported_153': dict(sorted(size_o.items())), 'class_sizes_171': dict(sorted(size_n.items())),
            'reported_combinations_changing_class_when_171_used': changes, 'top15_reported_classes': top15, 'distinct_scores_reported_153': len(set(round(v, 9) for v in old.values())),
            'max_abs_score_difference_shared': max(abs(old[c] - new[c]) for c in old)}


def part3():
    dl = drug_list(); apps_rows, drug_rows = [], []
    cache = {}
    for r in dl:
        d = r['drug']; q = query_name(d)
        if q not in cache:
            cache[q] = listing(q)
        total, res = cache[q]
        nb = [a for a in res if a['application_number'].startswith(('NDA', 'BLA'))]
        def first_ap(a):
            ds = sorted(s['submission_status_date'] for s in a.get('submissions', []) if s.get('submission_type') == 'ORIG' and s.get('submission_status') == 'AP')
            return ds[0] if ds else ''
        stat = {}
        for a in res:
            for p in a.get('products', []):
                if any(q in (i.get('name', '') or '').lower() for i in p.get('active_ingredients', [])):
                    stat[p.get('marketing_status', 'not stated')] = stat.get(p.get('marketing_status', 'not stated'), 0) + 1
        fa = sorted(x for x in (first_ap(a) for a in nb) if x)
        for a in nb:
            sts = sorted({p.get('marketing_status', 'not stated') for p in a.get('products', [])})
            apps_rows.append({'drug': d, 'application': a['application_number'], 'sponsor': a.get('sponsor_name', ''), 'brands': '|'.join(sorted({p.get('brand_name', '') for p in a.get('products', [])})),
                              'marketing_status_of_products': '|'.join(sts), 'first_ORIG_AP_date': first_ap(a)})
        marketed = sum(v for k, v in stat.items() if k in ('Prescription', 'Over-the-counter'))
        tf = r['FDA_EMA_PMDA_new'].strip('()').replace(' ', '').split(',')
        drug_rows.append({'drug': d, 'queried_name': q, 'applications_total_in_Drugs_at_FDA': total, 'NDA_or_BLA_applications': len(nb), 'products_by_marketing_status': json.dumps(stat, sort_keys=True),
                          'products_not_discontinued': marketed, 'US_status': ('no FDA application' if total == 0 else ('listed as marketed (not discontinued)' if marketed else 'application(s) present, all products discontinued or tentative')),
                          'first_US_NDA_BLA_approval': fa[0] if fa else '', 'NDA_BLA_brands': '|'.join(sorted({p.get('brand_name', '') for a in nb for p in a.get('products', [])})),
                          'DrugCentral_FDA': tf[0], 'DrugCentral_EMA': tf[1], 'DrugCentral_PMDA': tf[2]})
    with open(os.path.join(RES, 'market_status_by_drug.csv'), 'w', encoding='utf8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(drug_rows[0].keys())); w.writeheader(); w.writerows(drug_rows)
    with open(os.path.join(RES, 'market_status_nda_bla_applications.csv'), 'w', encoding='utf8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(apps_rows[0].keys())); w.writeheader(); w.writerows(apps_rows)
    dates = [json.load(open(os.path.join(RAW, fn), encoding='utf8'))['retrieved'] for fn in os.listdir(RAW) if fn.endswith('.json')]
    return drug_rows, apps_rows, {'responses_saved': len(dates), 'first_retrieved': min(dates), 'last_retrieved': max(dates)}


def analyze():
    h = part1(); q = part2(); d, a, meta = part3()
    out = {'hubs': h, 'quartiles': q, 'market': d, 'market_retrieval': meta}
    json.dump(out, open(os.path.join(RES, 'summary.json'), 'w', encoding='utf8'), indent=1, ensure_ascii=False)
    print(json.dumps({'hubs': h, 'quartiles': {k: v for k, v in q.items()}, 'retrieval': meta}, indent=1, ensure_ascii=False)[:5000])
    print('\nmarket status:')
    for r in d:
        print('  %-18s apps %3s NDA/BLA %2s  %-60s first NDA/BLA approval %s  brands %s' % (r['drug'], r['applications_total_in_Drugs_at_FDA'], r['NDA_or_BLA_applications'], r['US_status'], r['first_US_NDA_BLA_approval'], r['NDA_BLA_brands'][:60]))


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    if cmd == 'retrieve': retrieve()
    elif cmd == 'analyze': analyze()
    else: print(__doc__)
