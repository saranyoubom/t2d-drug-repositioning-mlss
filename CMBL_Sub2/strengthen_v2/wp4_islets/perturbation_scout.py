"""perturbation_scout.py - protocol 4.6, step 1: run the two fixed GEO DataSets queries (PRESPEC_WP4_choices.md section D) once through the NCBI E-utilities (db=gds) and list EVERY candidate series with
title, organism, assay type, number of samples, platform and the first 400 characters of the summary, plus an automatic stem-cell-design flag (regex on title + summary). Nothing is dropped here; the keep / drop
decision with its reason is written afterwards by hand from title, summary and design text only (no expression value is read at this step) into candidates_decisions.csv.
Plain HTTPS, at most 3 requests per second, plain tool User-Agent. Output: perturbation/candidates.csv, perturbation/queries.json.
Usage: python -I perturbation_scout.py"""
import os, sys, re, json, urllib.parse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import pandas as pd
import wp4_http as H
from attempts_log import record

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'perturbation'); os.makedirs(OUT, exist_ok=True)
CELL = '(islet OR islets OR "beta cell" OR "beta cells" OR "beta-cell" OR INS-1 OR INS1 OR MIN6 OR "pancreatic beta")'
QUERIES = {
    'notch': 'gse[ETYP] AND (Notch OR "gamma-secretase" OR "gamma secretase" OR DAPT) AND ' + CELL,
    'wnt': 'gse[ETYP] AND (DKK1 OR Dickkopf OR "Wnt inhibitor" OR "beta-catenin" OR TCF7L2 OR TCF4 OR "Wnt signaling" OR "Wnt signalling") AND ' + CELL,
}
EU = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
STEM = re.compile(r'differentiat|iPSC|\bESC\b|embryonic stem|pluripoten|stem cell|stem-cell|progenitor|reprogram', re.I)
rows, qlog = {}, {}
for name, term in QUERIES.items():
    url = EU + 'esearch.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'term': term, 'retmax': 500, 'retmode': 'json'})
    j = json.loads(H.get_text(url))['esearchresult']; ids = j['idlist']; qlog[name] = dict(term=term, count=int(j['count']), ids=ids, translation=j.get('querytranslation'))
    record('NCBI GEO DataSets', 'query ' + name, url[:250], None, 'E-utilities esearch db=gds', 'queried', 'count %s' % j['count'])
    print(name, 'count', j['count'], flush=True)
    for k in range(0, len(ids), 100):
        u = EU + 'esummary.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'id': ','.join(ids[k:k + 100]), 'retmode': 'json'}); res = json.loads(H.get_text(u))['result']
        for uid in res['uids']:
            r = res[uid]; acc = r.get('accession', '')
            d = rows.setdefault(acc, dict(accession=acc, entrytype=r.get('entrytype'), title=r.get('title'), organism=r.get('taxon'), assay=r.get('gdstype'), n_samples=r.get('n_samples'), platform=r.get('gpl'),
                                         queries=set(), summary=(r.get('summary') or '').replace('\n', ' ')))
            d['queries'].add(name)
D = pd.DataFrame(rows.values()); D['queries'] = D.queries.apply(lambda s: '+'.join(sorted(s)))
D['stem_cell_flag'] = [bool(STEM.search((a or '') + ' ' + (b or ''))) for a, b in zip(D.title, D.summary)]; D['summary'] = D.summary.str.slice(0, 400)
D = D.sort_values('accession'); D.to_csv(os.path.join(OUT, 'candidates.csv'), index=False); json.dump(qlog, open(os.path.join(OUT, 'queries.json'), 'w'), indent=1)
print(len(D), 'unique records;', D.entrytype.value_counts().to_dict()); print(D[['accession', 'queries', 'organism', 'assay', 'n_samples', 'stem_cell_flag', 'title']].to_string(index=False, max_colwidth=90))
