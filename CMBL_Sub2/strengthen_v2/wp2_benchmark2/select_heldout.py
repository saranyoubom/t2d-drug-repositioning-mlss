"""select_heldout.py - protocol v2 section 3.5, held-out microarrays (WP2 deliverable D). Blind selection: only sample METADATA (series-matrix headers) are read, no expression value.
Copy of the selection logic of ../../benchmark_route3/select_datasets.py (same query, numeric sort and shuffle with random.Random(2026), same parsing, perturbation keywords, title fallback,
control-like first group, 3 groups of 4 samples) with ONE change: the screening window is positions 151-450 of the same shuffled lists (GPL570 and GPL1261) instead of 1-150.
Added rules (declared before the lists were retrieved): (1) the re-run query must reproduce the part 1 lists - positions 1-150 of each shuffled list must equal the accessions in
benchmark_route3/eligibility.csv, otherwise the script stops; (2) a series that shares at least one sample (GSM accession, any sample of the series) with one of the 37 analysed development
series is excluded and listed with that reason; (3) NCBI requests are throttled to at most 3 per second.
Writes eligibility_heldout.csv (every screened series with status or reason), frozen_heldout_arrays.csv (ordered eligible list, groups, 12 sample accessions) and appends the SHA-256 of the frozen list
to ../PROTOCOL_v2_SHA256.txt. Headers are cached in <data root>/v2_data/arrays_heldout/headers (untrusted data, never executed).
Usage: python -I select_heldout.py [--check-only | --download-only]   (check-only: re-run the query and compare with part 1; download-only: cache the metadata headers of the window without screening or freezing;
no flag: screen from the cached headers, write the two csv files and the fingerprint)"""
import sys, os, re, gzip, csv, json, time, random, hashlib, urllib.parse, urllib.request, urllib.error
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__)); AO = os.path.dirname(os.path.dirname(HERE)); P1 = os.path.join(AO, 'benchmark_route3'); V2 = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import ncbi_http as H
DEFAULT_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_cache')
def data_root(): return os.environ.get('V2_DATA_ROOT', DEFAULT_ROOT)
EU ='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'; FTP = 'https://ftp.ncbi.nlm.nih.gov/geo/series/'
PLATFORMS = ['GPL570', 'GPL1261']; WINDOW = (150, 450)                               # python slice of the shuffled list: positions 151-450
QUERY = '{gpl}[ACCN] AND gse[ETYP] AND 12:36[n_samples] AND (treated OR treatment OR knockdown OR inhibitor OR stimulated OR overexpression OR siRNA OR shRNA OR agonist OR cytokine)'
KEYWORDS = ['treatment', 'agent', 'compound', 'drug', 'stimulus', 'stimulation', 'condition', 'genotype', 'knockdown', 'knockout', 'sirna', 'shrna', 'transfection', 'infection', 'cytokine', 'ligand', 'inhibitor', 'vehicle', 'dose', 'protocol', 'variation']
CONTROL = re.compile(r'\b(control|ctrl|vehicle|dmso|untreated|mock|nt|scrambled|wild|wt|naive|none|sham|pbs|saline|nontarget)\b', re.I)
TITLE_TAIL = re.compile(r'(?i)[\s_\-\.,#]*(rep(licate)?|sample|replic)?[\s_\-\.#]*\d{1,2}\s*$')


def get(url, binary=False, stream=False):
    return H.http_get(url, binary=binary, timeout=90, stream=stream)


def candidates(gpl):
    term = QUERY.format(gpl=gpl)
    j = json.loads(get(EU + 'esearch.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'term': term, 'retmax': 5000, 'retmode': 'json'})))
    ids = j['esearchresult']['idlist']; accs = []
    for i in range(0, len(ids), 150):
        s = json.loads(get(EU + 'esummary.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'id': ','.join(ids[i:i + 150]), 'retmode': 'json'})))['result']
        accs += [s[x]['accession'] for x in ids[i:i + 150] if x in s and str(s[x].get('accession', '')).startswith('GSE')]
    accs = sorted(set(accs), key=lambda a: int(a[3:])); random.Random(2026).shuffle(accs)
    return int(j['esearchresult']['count']), accs


def header(gse, gpl, hdir):
    """metadata lines of the series matrix (Sample_*, Series_title, Series_geo_accession); same result as select_datasets.header, one request instead of two when the plain file name exists"""
    path = os.path.join(hdir, f'{gse}_{gpl}.txt')
    if os.path.exists(path): return open(path, encoding='utf8').read().split('\n')
    n = int(gse[3:]); folder = f'GSE{n // 1000}nnn' if n >= 1000 else 'GSEnnn'; base = f'{FTP}{folder}/{gse}/matrix/'
    try: r = get(base + f'{gse}_series_matrix.txt.gz', stream=True)
    except urllib.error.HTTPError:
        idx = get(base); names = re.findall(r'href="([^"]*series_matrix\.txt\.gz)"', idx); hit = [x for x in names if gpl in x]
        if not hit: return None
        r = get(base + hit[0], stream=True)
    lines = []
    with gzip.GzipFile(fileobj=r) as g:
        for raw in g:
            line = raw.decode('utf8', 'replace').rstrip('\r\n')
            if line.startswith('!series_matrix_table_begin'): break
            if line.startswith(('!Sample_', '!Series_title', '!Series_geo_accession')): lines.append(line)
    r.close(); open(path, 'w', encoding='utf8').write('\n'.join(lines)); return lines


def parse(lines):
    rows = {}
    for l in lines:
        parts = l.split('\t'); rows.setdefault(parts[0], []).append([x.strip('"') for x in parts[1:]])
    gsm = rows.get('!Sample_geo_accession', [[]])[0]; title = rows.get('!Sample_title', [[]])[0]; plat = set(rows.get('!Sample_platform_id', [[]])[0])
    keys, per = [], [dict() for _ in gsm]
    for line in rows.get('!Sample_characteristics_ch1', []):
        for i, v in enumerate(line[:len(gsm)]):
            k, _, val = v.partition(': ')
            if not _: k, val = 'unnamed', v
            k = k.strip().lower()
            if k not in keys: keys.append(k)
            per[i][k] = val.strip()
    return gsm, title, plat, keys, per


def groups_from(labels):
    """label -> sample indices (file order), labels in order of first appearance"""
    d = {}
    for i, lab in enumerate(labels):
        if lab: d.setdefault(lab, []).append(i)
    return d


def pick_groups(gsm, title, keys, per):
    """perturbation key (first qualifying in file order) or title fallback; 3 groups of >= 4 samples; control-like group first. Returns (info dict or None)"""
    chosen, how = None, None
    for k in keys:
        if any(w in k for w in KEYWORDS):
            d = groups_from([p.get(k) for p in per]); q = {l: ix for l, ix in d.items() if len(ix) >= 4}
            if len(q) >= 3: chosen, how = q, 'characteristic: ' + k; break
    if chosen is None:
        d = groups_from([TITLE_TAIL.sub('', t) for t in title]); q = {l: ix for l, ix in d.items() if len(ix) >= 4}
        if len(q) >= 3: chosen, how = q, 'title'
    if chosen is None: return None
    labels = list(chosen); ctrl = [l for l in labels if CONTROL.search(l.replace('_', ' '))]
    first = ctrl[0] if ctrl else sorted(labels)[0]; sel = [first] + [l for l in labels if l != first][:2]
    return dict(how=how, n_groups=len(chosen), groups=sel, gsms=[[gsm[i] for i in chosen[l][:4]] for l in sel])


def screen(lines, gpl):
    gsm, title, plat, keys, per = parse(lines); n = len(gsm)
    if not (12 <= n <= 36): return 'ineligible: %d samples' % n, None
    if plat != {gpl}: return 'ineligible: platforms %s' % sorted(plat), None
    info = pick_groups(gsm, title, keys, per)
    if info is None: return 'ineligible: no grouping with >= 3 groups of >= 4 samples', None
    return 'eligible', info


def all_gsm(lines):
    for l in lines:
        if l.startswith('!Sample_geo_accession'): return [x.strip('"') for x in l.split('\t')[1:]]
    return []


def dev_sample_sets(root):
    """GSM accessions (all samples of the series header) of the 37 analysed development series -> {gsm: gse}"""
    ref = [json.loads(l) for l in open(os.path.join(P1, 'results_per_dataset.jsonl'), encoding='utf8')]; out = {}
    for r in ref:
        if r['status'] != 'analysed': continue
        p = os.path.join(root, 'bench_data', 'headers', '%s_%s.txt' % (r['gse'], r['platform']))
        for g in all_gsm(open(p, encoding='utf8').read().split('\n')): out.setdefault(g, r['gse'])
    return out


def main():
    root = data_root(); hdir =os.path.join(root, 'v2_data', 'arrays_heldout', 'headers'); os.makedirs(hdir, exist_ok=True)
    check_only = '--check-only' in sys.argv; dev = dev_sample_sets(root); print('development samples (analysed series, all samples):', len(dev), flush=True)
    p1 = {}
    for r in csv.DictReader(open(os.path.join(P1, 'eligibility.csv'), encoding='utf8')): p1.setdefault(r['platform'], []).append(r['gse'])
    elig, frozen, lists = [], [], {}
    for gpl in PLATFORMS:
        count, accs = candidates(gpl); lists[gpl] = accs; same = accs[:150] == p1[gpl]
        print(f'{gpl}: {count} series returned ({len(accs)} accessions); positions 1-150 identical to part 1: {same}', flush=True)
        if not same: raise SystemExit('the re-run query does not reproduce the part 1 list for %s: stop (document and decide before continuing)' % gpl)
    if check_only: return
    os.makedirs(os.path.join(HERE, 'compact'), exist_ok=True)
    json.dump({g: lists[g][WINDOW[0]:WINDOW[1]] for g in PLATFORMS}, open(os.path.join(HERE, 'compact', 'heldout_arrays_screened_order.json'), 'w'))
    if '--download-only' in sys.argv:                                      # cache the metadata headers only (no screening, no freezing)
        for gpl in PLATFORMS:
            for pos in range(WINDOW[0], min(WINDOW[1], len(lists[gpl]))):
                try: header(lists[gpl][pos], gpl, hdir)
                except Exception as e: print('header failed', lists[gpl][pos], type(e).__name__, flush=True)
                if (pos + 1 - WINDOW[0]) % 25 == 0: print(f'  {gpl}: headers cached up to position {pos + 1}', flush=True)
        print('headers cached', flush=True); return
    for gpl in PLATFORMS:
        accs = lists[gpl]
        for pos in range(WINDOW[0], min(WINDOW[1], len(accs))):
            gse = accs[pos]; extra = ''
            try:
                lines = header(gse, gpl, hdir)
                if lines is None: status, info = 'ineligible: no series matrix for the platform', None
                else:
                    status, info = screen(lines, gpl)
                    if info is not None:
                        shared = sorted({dev[g] for g in all_gsm(lines) if g in dev}); nsh = sum(1 for g in all_gsm(lines) if g in dev)
                        if shared: status, info = 'excluded: shares %d samples with development series %s' % (nsh, ','.join(shared)), None
            except Exception as e: status, info = 'ineligible: error %s' % type(e).__name__, None
            elig.append(dict(platform=gpl, order=pos + 1, gse=gse, status=status, how=(info or {}).get('how', ''), n_groups=(info or {}).get('n_groups', ''), groups='|'.join((info or {}).get('groups', []))))
            if info: frozen.append(dict(platform=gpl, order=pos + 1, gse=gse, how=info['how'], n_groups=info['n_groups'], groups='|'.join(info['groups']), gsms=';'.join(sum(info['gsms'], []))))
            if (pos + 1 - WINDOW[0]) % 25 == 0: print(f'  {gpl}: positions {WINDOW[0] + 1}-{pos + 1} screened, {sum(1 for f in frozen if f["platform"] == gpl)} eligible', flush=True)
    for name, rows in (('eligibility_heldout.csv', elig), ('frozen_heldout_arrays.csv', frozen)):
        with open(os.path.join(HERE, name), 'w', newline='', encoding='utf8') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    h = hashlib.sha256(open(os.path.join(HERE, 'frozen_heldout_arrays.csv'), 'rb').read()).hexdigest(); stamp = time.strftime('%Y-%m-%d %H:%M:%S %z')
    with open(os.path.join(V2, 'PROTOCOL_v2_SHA256.txt'), 'a', encoding='utf8') as f: f.write(f'{h} *wp2_benchmark2/frozen_heldout_arrays.csv\nWP2 held-out microarray list frozen BEFORE any expression value was read: {stamp}\n')
    print('eligible per platform:', {g: sum(1 for f in frozen if f['platform'] == g) for g in PLATFORMS}, '| sha256', h, '|', stamp, flush=True)
    print('status counts:', {s: sum(1 for e in elig if e['status'] == s) for s in sorted({e['status'] for e in elig})}, flush=True)


if __name__ == '__main__':
    main()
