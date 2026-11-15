"""Blind selection of the benchmark datasets (PROTOCOL sections 3 and 9). Only sample METADATA (series-matrix header) are read here; no expression value is touched.
Writes eligibility.csv (every screened series with its status or reason) and frozen_datasets.csv (the ordered eligible list with the chosen groups and the 12 sample accessions), then their SHA-256.
Downloaded headers are cached in <data dir>/headers (untrusted data, never executed). Usage: python -I select_datasets.py <data dir>"""
import sys, os, re, gzip, csv, json, time, random, hashlib, urllib.parse, urllib.request, urllib.error
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = sys.argv[1]; os.makedirs(os.path.join(DATA, 'headers'), exist_ok=True)
EU = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'; FTP = 'https://ftp.ncbi.nlm.nih.gov/geo/series/'
PLATFORMS = ['GPL570', 'GPL1261']; N_SCREEN = 150
QUERY = '{gpl}[ACCN] AND gse[ETYP] AND 12:36[n_samples] AND (treated OR treatment OR knockdown OR inhibitor OR stimulated OR overexpression OR siRNA OR shRNA OR agonist OR cytokine)'
KEYWORDS = ['treatment', 'agent', 'compound', 'drug', 'stimulus', 'stimulation', 'condition', 'genotype', 'knockdown', 'knockout', 'sirna', 'shrna', 'transfection', 'infection', 'cytokine', 'ligand', 'inhibitor', 'vehicle', 'dose', 'protocol', 'variation']
CONTROL = re.compile(r'\b(control|ctrl|vehicle|dmso|untreated|mock|nt|scrambled|wild|wt|naive|none|sham|pbs|saline|nontarget)\b', re.I)
TITLE_TAIL = re.compile(r'(?i)[\s_\-\.,#]*(rep(licate)?|sample|replic)?[\s_\-\.#]*\d{1,2}\s*$')


def get(url, binary=False, tries=4, stream=False):
    for k in range(tries):
        try:
            r = urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'geo-benchmark/1.0'}), timeout=90)
            return r if stream else (r.read() if binary else r.read().decode('utf8', 'replace'))
        except urllib.error.HTTPError as e:
            if e.code == 404: raise
            time.sleep(2 + 2 * k)
        except Exception: time.sleep(2 + 2 * k)
    raise RuntimeError('failed ' + url)


def candidates(gpl):
    term = QUERY.format(gpl=gpl)
    j = json.loads(get(EU + 'esearch.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'term': term, 'retmax': 5000, 'retmode': 'json'}))); time.sleep(0.4)
    ids = j['esearchresult']['idlist']; accs = []
    for i in range(0, len(ids), 150):
        s = json.loads(get(EU + 'esummary.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'id': ','.join(ids[i:i + 150]), 'retmode': 'json'})))['result']; time.sleep(0.4)
        accs += [s[x]['accession'] for x in ids[i:i + 150] if x in s and str(s[x].get('accession', '')).startswith('GSE')]
    accs = sorted(set(accs), key=lambda a: int(a[3:])); random.Random(2026).shuffle(accs)
    return int(j['esearchresult']['count']), accs


def matrix_url(gse, gpl):
    n = int(gse[3:]); folder = f'GSE{n // 1000}nnn' if n >= 1000 else 'GSEnnn'; base = f'{FTP}{folder}/{gse}/matrix/'
    try:
        get(base + f'{gse}_series_matrix.txt.gz', stream=True).close(); return base + f'{gse}_series_matrix.txt.gz'
    except urllib.error.HTTPError:
        idx = get(base); names = re.findall(r'href="([^"]*series_matrix\.txt\.gz)"', idx); hit = [x for x in names if gpl in x]
        return base + hit[0] if hit else None


def header(gse, gpl):
    path = os.path.join(DATA, 'headers', f'{gse}_{gpl}.txt')
    if os.path.exists(path): return open(path, encoding='utf8').read().split('\n')
    url = matrix_url(gse, gpl)
    if url is None: return None
    r = get(url, stream=True); lines = []
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


def screen(lines, gpl):
    gsm, title, plat, keys, per = parse(lines); n = len(gsm)
    if not (12 <= n <= 36): return 'ineligible: %d samples' % n, None
    if plat != {gpl}: return 'ineligible: platforms %s' % sorted(plat), None
    chosen, how = None, None
    for k in keys:
        if any(w in k for w in KEYWORDS):
            d = groups_from([p.get(k) for p in per]); q = {l: ix for l, ix in d.items() if len(ix) >= 4}
            if len(q) >= 3: chosen, how = q, 'characteristic: ' + k; break
    if chosen is None:
        d = groups_from([TITLE_TAIL.sub('', t) for t in title]); q = {l: ix for l, ix in d.items() if len(ix) >= 4}
        if len(q) >= 3: chosen, how = q, 'title'
    if chosen is None: return 'ineligible: no grouping with >= 3 groups of >= 4 samples', None
    labels = list(chosen); ctrl = [l for l in labels if CONTROL.search(l.replace('_', ' '))]
    first = ctrl[0] if ctrl else sorted(labels)[0]; sel = [first] + [l for l in labels if l != first][:2]
    return 'eligible', dict(how=how, n_groups=len(chosen), groups=sel, gsms=[[gsm[i] for i in chosen[l][:4]] for l in sel])


def main():
    elig, frozen = [], []
    for gpl in PLATFORMS:
        count, accs = candidates(gpl); print(f'{gpl}: {count} series returned; screening the first {N_SCREEN} of the shuffled list', flush=True)
        for pos, gse in enumerate(accs[:N_SCREEN]):
            try:
                lines = header(gse, gpl)
                if lines is None: status, info = 'ineligible: no series matrix for the platform', None
                else: status, info = screen(lines, gpl)
            except Exception as e: status, info = 'ineligible: error %s' % type(e).__name__, None
            elig.append(dict(platform=gpl, order=pos + 1, gse=gse, status=status, how=(info or {}).get('how', ''), n_groups=(info or {}).get('n_groups', ''), groups='|'.join((info or {}).get('groups', []))))
            if info: frozen.append(dict(platform=gpl, order=pos + 1, gse=gse, how=info['how'], n_groups=info['n_groups'], groups='|'.join(info['groups']), gsms=';'.join(sum(info['gsms'], []))))
            if (pos + 1) % 25 == 0: print(f'  {gpl}: {pos + 1} screened, {sum(1 for f in frozen if f["platform"] == gpl)} eligible', flush=True)
            time.sleep(0.3)
    for name, rows in (('eligibility.csv', elig), ('frozen_datasets.csv', frozen)):
        with open(os.path.join(HERE, name), 'w', newline='', encoding='utf8') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    h = hashlib.sha256(open(os.path.join(HERE, 'frozen_datasets.csv'), 'rb').read()).hexdigest()
    open(os.path.join(HERE, 'FROZEN_SHA256.txt'), 'w').write(f'{h} *frozen_datasets.csv\n{time.strftime("frozen %Y-%m-%d %H:%M:%S %z")}\n')
    print('eligible per platform:', {g: sum(1 for f in frozen if f['platform'] == g) for g in PLATFORMS}, '| sha256', h, flush=True)


if __name__ == '__main__':
    main()
