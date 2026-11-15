"""select_rnaseq.py - protocol v2 section 3.5, held-out RNA-seq (WP2 deliverable E): candidate lists, pattern check of the NCBI-generated count download, metadata screening and freezing.
Query (verbatim from the protocol): ("Homo sapiens"[Organism] OR "Mus musculus"[Organism]) AND "Expression profiling by high throughput sequencing"[DataSet Type] AND gse[ETYP] AND 12:36[n_samples]
AND (treated OR treatment OR knockdown OR inhibitor OR stimulated OR overexpression OR siRNA OR shRNA OR agonist OR cytokine). The set has more than 30000 records, above the 10000 that one
esearch call can return, so it is retrieved with the same query split by n_samples (12, 13, ... 36; the union is the same set). A GEO series record has UID 200000000 + GSE number (checked on the records of
the earlier steps), so accessions come from the UIDs; the organism comes from running the query once per organism; series found under both organisms are excluded as mixed. Per organism the accessions
are sorted numerically and shuffled with random.Random(2026); the first 300 are screened on metadata by the rules of select_heldout.py (12-36 samples, ONE platform for the whole series, same
perturbation keys / title fallback, 3 groups of 4, control-like group first).
Steps (python -I select_rnaseq.py <step>):
  lists           build and save the shuffled lists (compact/rnaseq_lists.json)
  pattern-check   status and content type (never the content) of the NCBI-generated count URL for 20 series per organism taken from positions 301 and later of the shuffled lists, i.e. outside the screened
                  window, plus the annotation file; written to rnaseq_pattern_check.csv. If the downloads are blocked or unavailable the arm is reported as not done and nothing further runs.
  screen          (only if the pattern check shows plain downloads work) screen the first 300 per organism on metadata, write eligibility_heldout_rnaseq.csv and frozen_heldout_rnaseq.csv and append the
                  SHA-256 to ../PROTOCOL_v2_SHA256.txt BEFORE any count is read."""
import sys, os, re, gzip, csv, json, time, random, hashlib, urllib.parse, urllib.error
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__)); V2 = os.path.dirname(HERE); sys.path.insert(0, HERE)
import ncbi_http as H
import select_heldout as SH
EU = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
TERMS = '"Expression profiling by high throughput sequencing"[DataSet Type] AND gse[ETYP] AND {lo}:{hi}[n_samples] AND (treated OR treatment OR knockdown OR inhibitor OR stimulated OR overexpression OR siRNA OR shRNA OR agonist OR cytokine)'
ORGS = {'human': 'Homo sapiens', 'mouse': 'Mus musculus'}; BUILD = {'human': 'GRCh38.p13', 'mouse': 'GRCm39'}; N_SCREEN = 300
COUNT_URL = 'https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc={acc}&format=file&file={acc}_raw_counts_{build}_NCBI.tsv.gz'
ANNOT_URL = {'human': 'https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz', 'mouse': 'https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Mouse.GRCm39.annot.tsv.gz'}
F_LISTS = os.path.join(HERE, 'compact', 'rnaseq_lists.json')


def uids(org_name, lo, hi):
    term = '"%s"[Organism] AND ' % org_name + TERMS.format(lo=lo, hi=hi)
    j = json.loads(H.http_get(EU + 'esearch.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'term': term, 'retmax': 10000, 'retmode': 'json'}), binary=False, timeout=90))['esearchresult']
    n = int(j['count']); ids = j['idlist']
    if n > len(ids): raise SystemExit('n_samples %d-%d for %s has %d records, more than one call returns: split further' % (lo, hi, org_name, n))
    return ids


def build_lists():
    out = {}; raw = {}
    for key, name in ORGS.items():
        s = set()
        for n in range(12, 37): s |= {int(u) - 200000000 for u in uids(name, n, n)}
        raw[key] = s; print(key, 'records', len(s), flush=True)
    mixed = raw['human'] & raw['mouse']
    for key in ORGS:
        accs = sorted(raw[key] - mixed); lst = ['GSE%d' % a for a in accs]; random.Random(2026).shuffle(lst); out[key] = lst
    out['n_records'] = {k: len(v) for k, v in raw.items()}; out['n_mixed'] = len(mixed); out['queried'] = time.strftime('%Y-%m-%d %H:%M:%S %z')
    os.makedirs(os.path.dirname(F_LISTS), exist_ok=True); json.dump(out, open(F_LISTS, 'w')); print('lists saved; mixed-organism series excluded:', len(mixed), flush=True); return out


def probe(url, nbytes=2):
    """status, content type and (if gzip) the first two bytes only: the content is never read"""
    try:
        r = H.http_get(url, stream=True, timeout=60, tries=2); ct = r.headers.get('Content-Type', ''); first = r.read(nbytes) if 'html' not in ct else r.read(60); r.close()
        kind = 'gzip file' if first[:2] == b'\x1f\x8b' else ('HTML (challenge or error page)' if b'<' in first[:60] else 'other')
        return dict(status=r.status, content_type=ct, kind=kind)
    except urllib.error.HTTPError as e: return dict(status=e.code, content_type='', kind='HTTP error')
    except Exception as e: return dict(status=-1, content_type='', kind=type(e).__name__)


def pattern_check(lists):
    rows = []
    for key in ORGS:
        for i, acc in enumerate(lists[key][N_SCREEN:N_SCREEN + 20]):
            p = probe(COUNT_URL.format(acc=acc, build=BUILD[key])); rows.append(dict(organism=key, position=N_SCREEN + i + 1, what='counts ' + acc, **p)); print(rows[-1], flush=True)
        p = probe(ANNOT_URL[key]); rows.append(dict(organism=key, position=0, what='annotation', **p)); print(rows[-1], flush=True)
    with open(os.path.join(HERE, 'rnaseq_pattern_check.csv'), 'w', newline='', encoding='utf8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    n_gz = sum(r['kind'] == 'gzip file' for r in rows); n_html = sum('HTML' in r['kind'] for r in rows); n_404 = sum(r['status'] == 404 for r in rows)
    print('pattern check: %d probes; plain gzip files %d, HTML pages %d, HTTP 404 %d' % (len(rows), n_gz, n_html, n_404), flush=True); return rows


def header_rna(gse, hdir):
    """metadata lines; None if the series has no single series-matrix file (several platforms)"""
    path = os.path.join(hdir, gse + '.txt')
    if os.path.exists(path):
        t = open(path, encoding='utf8').read(); return None if t == 'MULTI' else t.split('\n')
    n = int(gse[3:]); folder = f'GSE{n // 1000}nnn' if n >= 1000 else 'GSEnnn'; url = f'{SH.FTP}{folder}/{gse}/matrix/{gse}_series_matrix.txt.gz'
    try: r = H.http_get(url, stream=True, timeout=90)
    except urllib.error.HTTPError as e:
        if e.code == 404: open(path, 'w').write('MULTI'); return None
        raise
    lines = []
    with gzip.GzipFile(fileobj=r) as g:
        for raw in g:
            line = raw.decode('utf8', 'replace').rstrip('\r\n')
            if line.startswith('!series_matrix_table_begin'): break
            if line.startswith(('!Sample_', '!Series_title', '!Series_geo_accession')): lines.append(line)
    r.close(); open(path, 'w', encoding='utf8').write('\n'.join(lines)); return lines


def screen_rna(lines):
    gsm, title, plat, keys, per = SH.parse(lines); n = len(gsm)
    if not (12 <= n <= 36): return 'ineligible: %d samples' % n, None
    if len(plat) != 1: return 'ineligible: platforms %s' % sorted(plat), None
    info = SH.pick_groups(gsm, title, keys, per)
    if info is None: return 'ineligible: no grouping with >= 3 groups of >= 4 samples', None
    return 'eligible', info


def screen(lists, orgs=None):
    orgs = orgs or list(ORGS)
    root = SH.data_root(); hdir = os.path.join(root, 'v2_data', 'rnaseq', 'headers'); os.makedirs(hdir, exist_ok=True); dev = SH.dev_sample_sets(root); elig, frozen = [], []
    arr = os.path.join(HERE, 'frozen_heldout_arrays.csv'); held = {}
    if os.path.exists(arr):
        for r in csv.DictReader(open(arr, encoding='utf8')):
            p = os.path.join(root, 'v2_data', 'arrays_heldout', 'headers', '%s_%s.txt' % (r['gse'], r['platform']))
            if os.path.exists(p):
                for g in SH.all_gsm(open(p, encoding='utf8').read().split('\n')): held.setdefault(g, r['gse'])
    for key in orgs:
        for pos, gse in enumerate(lists[key][:N_SCREEN]):
            try:
                lines = header_rna(gse, hdir)
                if lines is None: status, info = 'ineligible: several platforms (no single series matrix)', None
                else:
                    status, info = screen_rna(lines)
                    if info is not None:
                        g_all = SH.all_gsm(lines); shared = sorted({dev[g] for g in g_all if g in dev} | {held[g] for g in g_all if g in held})
                        if shared: status, info = 'excluded: shares samples with series %s' % ','.join(shared), None
            except Exception as e: status, info = 'ineligible: error %s' % type(e).__name__, None
            elig.append(dict(organism=key, order=pos + 1, gse=gse, status=status, how=(info or {}).get('how', ''), n_groups=(info or {}).get('n_groups', ''), groups='|'.join((info or {}).get('groups', []))))
            if info: frozen.append(dict(organism=key, order=pos + 1, gse=gse, how=info['how'], n_groups=info['n_groups'], groups='|'.join(info['groups']), gsms=';'.join(sum(info['gsms'], []))))
            if (pos + 1) % 25 == 0: print(f'  {key}: {pos + 1} screened, {sum(1 for f in frozen if f["organism"] == key)} eligible', flush=True)
    for name, rows in (('eligibility_heldout_rnaseq.csv', elig), ('frozen_heldout_rnaseq.csv', frozen)):
        with open(os.path.join(HERE, name), 'w', newline='', encoding='utf8') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    h = hashlib.sha256(open(os.path.join(HERE, 'frozen_heldout_rnaseq.csv'), 'rb').read()).hexdigest(); stamp = time.strftime('%Y-%m-%d %H:%M:%S %z')
    with open(os.path.join(V2, 'PROTOCOL_v2_SHA256.txt'), 'a', encoding='utf8') as f: f.write(f'{h} *wp2_benchmark2/frozen_heldout_rnaseq.csv\nWP2 held-out RNA-seq list frozen BEFORE any count was read: {stamp}\n')
    print('eligible per organism:', {k: sum(1 for f in frozen if f['organism'] == k) for k in orgs}, '| sha256', h, '|', stamp, flush=True)


if __name__ == '__main__':
    step = sys.argv[1] if len(sys.argv) > 1 else 'lists'
    if step == 'lists': build_lists()
    elif step == 'pattern-check': pattern_check(json.load(open(F_LISTS)))
    elif step == 'screen': screen(json.load(open(F_LISTS)), sys.argv[2].split(',') if len(sys.argv) > 2 else None)         # e.g. screen human
    else: raise SystemExit('lists | pattern-check | screen [human,mouse]')
