"""wider_search.py - B1 of PRESPEC_WP4b_positive_control_and_wider_search.md (queries fixed before this ran): three sources, each run once, every returned record listed with the source(s) that returned it.
 S1  NCBI GEO DataSets (db=gds):  gse[ETYP] AND (NOTCH OR WNT) AND (CELL)
 S2  EBI BioStudies, ArrayExpress collection (api/v1/arrayexpress/search), the same three term groups as one query (E-GEOD-n records are the GEO series GSEn and are merged with S1)
 S3  PubMed: (CELL) AND (NOTCH OR WNT) AND (transcriptom* OR "RNA-seq" OR "RNA sequencing" OR microarray OR "expression profiling" OR "gene expression"), all PMIDs (371 <= the cap of 1500), linked to GEO by elink (pubmed_gds)
Plain HTTPS, at most 3 requests per second, plain tool User-Agent, nothing about the user is sent; a refusal or challenge page would be reported, not worked around. No expression value is read.
Output: wider_search/queries.json (queries, counts, date), wider_search/candidates_all.csv (accession, sources, pmids, entry type, title, organism, assay, n_samples, platform, summary).   Usage: python -I wider_search.py"""
import os, sys, re, json, time, urllib.parse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import pandas as pd
import wp4_http as H
from attempts_log import record

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'wider_search'); os.makedirs(OUT, exist_ok=True)
CELL = '(islet OR islets OR "beta cell" OR "beta cells" OR "beta-cell" OR "INS-1" OR INS1 OR "INS-1E" OR MIN6 OR EndoC OR betaTC OR "BRIN-BD11" OR RINm5F OR insulinoma OR "pancreatic beta")'
NOTCH = '(Notch OR "gamma-secretase" OR "gamma secretase" OR DAPT OR LY411575 OR RO4929097 OR "compound E" OR "MK-0752" OR RBPJ OR Jagged OR Jag1 OR Dll1 OR Dll4 OR Notch1 OR Notch2)'
WNT = '(Wnt OR DKK1 OR Dickkopf OR "beta-catenin" OR CTNNB1 OR TCF7L2 OR TCF4 OR IWP OR "IWP-2" OR IWR OR XAV939 OR LGK974 OR "ICG-001" OR CHIR99021 OR GSK3 OR Wnt3a OR Axin OR SFRP OR Frizzled)'
PERT = f'({NOTCH} OR {WNT})'
Q1 = f'gse[ETYP] AND {PERT} AND {CELL}'; Q2 = f'{CELL} AND {PERT}'
Q3 = f'{CELL} AND {PERT} AND (transcriptom* OR "RNA-seq" OR "RNA sequencing" OR microarray OR "expression profiling" OR "gene expression")'
EU = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'; BS = 'https://www.ebi.ac.uk/biostudies/api/v1/'
cand = {}; qlog = dict(date=time.strftime('%Y-%m-%d %H:%M:%S %z'), terms=dict(CELL=CELL, NOTCH=NOTCH, WNT=WNT), S1=dict(query=Q1), S2=dict(query=Q2), S3=dict(query=Q3, cap=1500))


def jget(url): return json.loads(H.get_text(url))


def add(acc, src, **kw):
    d = cand.setdefault(acc, dict(accession=acc, sources=set(), pmids=set(), entrytype=None, title=None, organism=None, assay=None, n_samples=None, platform=None, summary=None))
    d['sources'].add(src)
    for k, v in kw.items():
        if k == 'pmid': d['pmids'].add(str(v))
        elif v not in (None, '') and not d.get(k): d[k] = v


def gds_records(uids):
    """esummary of GEO DataSets uids -> add() for series (GSE); curated datasets (GDS) are mapped to their series"""
    got = 0
    for k in range(0, len(uids), 100):
        r = jget(EU + 'esummary.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'id': ','.join(uids[k:k + 100]), 'retmode': 'json'}))['result']
        for uid in r['uids']:
            x = r[uid]; acc = x.get('accession', '')
            if acc.startswith('GDS') and x.get('gse'): acc = 'GSE' + str(x['gse'])
            if not acc.startswith('GSE'): continue
            yield acc, dict(entrytype=x.get('entrytype'), title=x.get('title'), organism=x.get('taxon'), assay=x.get('gdstype'), n_samples=x.get('n_samples'), platform=x.get('gpl'), summary=(x.get('summary') or '').replace('\n', ' ')[:700]); got += 1


# ---- S1 GEO
j = jget(EU + 'esearch.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'term': Q1, 'retmax': 5000, 'retmode': 'json'}))['esearchresult']; ids1 = j['idlist']; qlog['S1'].update(count=int(j['count']), ids=len(ids1))
for acc, kw in gds_records(ids1): add(acc, 'S1', **kw)
record('NCBI GEO DataSets', 'wider search S1', EU + 'esearch.fcgi (db=gds)', None, 'E-utilities', 'queried', 'count %s' % j['count']); print('S1 GEO: count', j['count'], '| series', sum(1 for d in cand.values() if 'S1' in d['sources']), flush=True)
# ---- S3 PubMed -> GEO
j = jget(EU + 'esearch.fcgi?' + urllib.parse.urlencode({'db': 'pubmed', 'term': Q3, 'retmax': 5000, 'retmode': 'json'}))['esearchresult']; pm = sorted(j['idlist'], key=int)[:1500]; qlog['S3'].update(count=int(j['count']), pmids=len(pm))
link = {}
for k in range(0, len(pm), 50):
    batch = pm[k:k + 50]; r = jget(EU + 'elink.fcgi?dbfrom=pubmed&db=gds&linkname=pubmed_gds&retmode=json&' + '&'.join('id=' + p for p in batch))
    for ls in r.get('linksets', []):
        for db in ls.get('linksetdbs', []):
            if db.get('linkname') == 'pubmed_gds':
                for u in db.get('links', []): link.setdefault(str(u), set()).update(str(i) for i in ls['ids'])
uids = sorted(link); qlog['S3'].update(linked_gds_uids=len(uids))
for acc, kw in gds_records(uids):
    pass
rec3 = {}
for k in range(0, len(uids), 100):
    r = jget(EU + 'esummary.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'id': ','.join(uids[k:k + 100]), 'retmode': 'json'}))['result']
    for uid in r['uids']:
        x = r[uid]; acc = x.get('accession', '')
        if acc.startswith('GDS') and x.get('gse'): acc = 'GSE' + str(x['gse'])
        if not acc.startswith('GSE'): continue
        add(acc, 'S3', entrytype=x.get('entrytype'), title=x.get('title'), organism=x.get('taxon'), assay=x.get('gdstype'), n_samples=x.get('n_samples'), platform=x.get('gpl'), summary=(x.get('summary') or '').replace('\n', ' ')[:700])
        for p in link[uid]: add(acc, 'S3', pmid=p)
record('PubMed -> GEO', 'wider search S3', EU + 'esearch.fcgi (db=pubmed) + elink pubmed_gds', None, 'E-utilities', 'queried', 'PMIDs %d, linked GEO uids %d' % (len(pm), len(uids))); print('S3 PubMed: count', j['count'], '| PMIDs', len(pm), '| linked GEO records', len(uids), '| series in the table', sum(1 for d in cand.values() if 'S3' in d['sources']), flush=True)
# ---- S2 BioStudies / ArrayExpress
hits, page, total = [], 1, None
while True:
    try: r = jget(BS + 'arrayexpress/search?' + urllib.parse.urlencode({'query': Q2, 'pageSize': 100, 'page': page}))
    except Exception as e: qlog['S2']['error'] = repr(e)[:300]; print('S2 FAILED (reported, not worked around):', repr(e)[:200]); break
    total = r.get('totalHits'); hits += r.get('hits', [])
    if len(hits) >= (total or 0) or not r.get('hits'): break
    page += 1
qlog['S2'].update(count=total, hits=len(hits)); ae = []
for h in hits:
    a = h['accession']; acc = 'GSE' + a[len('E-GEOD-'):] if a.startswith('E-GEOD-') else a
    add(acc, 'S2', title=h.get('title'));
    if not acc.startswith('GSE'): ae.append(acc)
for a in ae:                                                                    # non-GEO ArrayExpress studies: description from the study record
    try:
        s = jget(BS + 'studies/' + a); attrs = {x.get('name'): x.get('value') for x in (s.get('section', {}).get('attributes', []) or [])}
        add(a, 'S2', entrytype='ArrayExpress', title=attrs.get('Title') or next((x.get('value') for x in s.get('attributes', []) if x.get('name') == 'Title'), None), summary=((attrs.get('Description') or attrs.get('Abstract') or '')).replace('\n', ' ')[:700],
            organism=attrs.get('Organism'), assay=attrs.get('Study type') or attrs.get('Assay'))
    except Exception as e: add(a, 'S2', summary='(study record not retrievable: %s)' % repr(e)[:80])
# GEO series reached only through ArrayExpress (E-GEOD) need their GEO metadata
need = [acc for acc, d in cand.items() if acc.startswith('GSE') and not d['title'] or acc.startswith('GSE') and not d['organism']]
for acc in need:
    j = jget(EU + 'esearch.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'term': acc + '[ACCN] AND gse[ETYP]', 'retmode': 'json'}))['esearchresult']['idlist']
    if j:
        for a2, kw in gds_records(j): add(a2, 'S2', **kw)
record('EBI BioStudies', 'wider search S2 (arrayexpress collection)', BS + 'arrayexpress/search', None, 'BioStudies API', 'queried', 'hits %s' % total); print('S2 BioStudies: hits', total, '| non-GEO studies', len(ae), flush=True)
rows = []
for acc, d in sorted(cand.items()):
    rows.append(dict(accession=acc, sources='+'.join(sorted(d['sources'])), pmids=';'.join(sorted(d['pmids'], key=int)), entrytype=d['entrytype'], title=d['title'], organism=d['organism'], assay=d['assay'], n_samples=d['n_samples'], platform=d['platform'], summary=d['summary']))
D = pd.DataFrame(rows); D.to_csv(os.path.join(OUT, 'candidates_all.csv'), index=False); qlog['unique_records'] = len(D); qlog['by_source_combination'] = D.sources.value_counts().to_dict()
json.dump(qlog, open(os.path.join(OUT, 'queries.json'), 'w'), indent=1); print(len(D), 'unique records;', qlog['by_source_combination'])
