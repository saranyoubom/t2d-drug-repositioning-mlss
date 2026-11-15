"""fetch_dataset_citations.py - the publications attached to the public datasets used by WP4, taken from the records themselves: the PubMed ids in the GEO series records (E-utilities esearch/esummary, db=gds) and the
'PubMed ID' line of the E-MTAB-5061 IDF, then title, authors, journal, year and DOI from PubMed (esummary, db=pubmed). Plain HTTPS, at most 3 requests per second, nothing about the user is sent. Output:
dataset_citations.csv. The references still have to be checked by the manuscript revision (WP6b) before they are cited.
Usage: python -I fetch_dataset_citations.py"""
import os, sys, json, re, urllib.parse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import pandas as pd
import wp4_http as H
import wp4_core as C
from attempts_log import record

HERE = os.path.dirname(os.path.abspath(__file__)); EU = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
SER = ['GSE81608', 'GSE84133', 'GSE86469', 'GSE81547', 'GSE60803']
pm = {}
for acc in SER:
    j = json.loads(H.get_text(EU + 'esearch.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'term': acc + '[ACCN] AND gse[ETYP]', 'retmode': 'json'})))['esearchresult']['idlist']
    r = json.loads(H.get_text(EU + 'esummary.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'id': ','.join(j), 'retmode': 'json'})))['result']
    pm[acc] = [x for uid in r['uids'] for x in r[uid].get('pubmedids', [])]
idf = open(os.path.join(os.path.dirname(C.PARSED), 'E-MTAB-5061', 'E-MTAB-5061.idf.txt'), encoding='utf8', errors='replace').read()
m = re.search(r'^PubMed ID\t([^\n]*)', idf, flags=re.M); pm['E-MTAB-5061'] = [x.strip() for x in m.group(1).split('\t') if x.strip()] if m else []
ids = sorted({p for v in pm.values() for p in v}); rows = []
S = json.loads(H.get_text(EU + 'esummary.fcgi?' + urllib.parse.urlencode({'db': 'pubmed', 'id': ','.join(ids), 'retmode': 'json'})))['result'] if ids else {}
for acc, lst in pm.items():
    for p in lst:
        r = S.get(p, {}); doi = next((a['value'] for a in r.get('articleids', []) if a.get('idtype') == 'doi'), '')
        rows.append(dict(dataset=acc, pmid=p, first_author=(r.get('authors') or [{}])[0].get('name', ''), title=r.get('title', ''), journal=r.get('source', ''), year=(r.get('pubdate') or '')[:4], volume=r.get('volume', ''), pages=r.get('pages', ''), doi=doi))
D = pd.DataFrame(rows); D.to_csv(os.path.join(HERE, 'dataset_citations.csv'), index=False); pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 80); print(D.to_string(index=False))
record('NCBI GEO / PubMed', 'publications attached to the dataset records', EU + '(esearch/esummary db=gds, db=pubmed)', None, 'GEO series records, E-MTAB-5061 IDF', 'queried', 'citations for the manuscript revision; to be verified there')
