"""fetch_ensembl_ids.py - Ensembl gene ids of the 17 human genes, read from the NCBI Datasets gene report (api.ncbi.nlm.nih.gov/datasets/v2/gene/id/<Entrez ids>, field ensembl_gene_ids) for the Entrez
ids already used in wp4_core.py. GSE86469 labels its genes by Ensembl id, so the 17 genes (and the eligible-gene pool) need this key. Plain HTTPS, one request, plain tool User-Agent, nothing about the user is sent.
(An efetch XML request of the Gene records was tried first and refused by the 20 MB response limit of wp4_http.get_text; it was not repeated.)
Output: ensembl_ids_17genes.csv (gene, entrez, ncbi_symbol, ensembl). The attempt is logged in attempts.csv.
Usage: python -I fetch_ensembl_ids.py"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import wp4_http as H
from attempts_log import record
import wp4_core as C
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ids = [C.ENTREZ[g] for g in C.GENES]
url = 'https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/' + '%2C'.join(str(i) for i in ids)
txt = H.get_text(url); rep = json.loads(txt)['reports']
by_id = {int(r['gene']['gene_id']): r['gene'] for r in rep}
rows = [dict(gene=g, entrez=C.ENTREZ[g], ncbi_symbol=by_id[C.ENTREZ[g]]['symbol'], tax_id=by_id[C.ENTREZ[g]]['tax_id'], ensembl=';'.join(by_id[C.ENTREZ[g]].get('ensembl_gene_ids', []))) for g in C.GENES]
D = pd.DataFrame(rows); D.to_csv(os.path.join(HERE, 'ensembl_ids_17genes.csv'), index=False)
print(D.to_string(index=False))
assert len(D) == 17 and all(D.ensembl.str.count('ENSG') == 1), 'every gene needs exactly one Ensembl id'
record('NCBI Gene', 'Ensembl ids of the 17 genes', url[:200], len(txt), 'NCBI Datasets gene report (ensembl_gene_ids)', 'downloaded', 'annotation only; key for GSE86469 gene labels')
