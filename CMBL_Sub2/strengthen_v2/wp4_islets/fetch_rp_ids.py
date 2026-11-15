"""fetch_rp_ids.py - Entrez and Ensembl gene ids of the ribosomal protein genes of rp_list.py, for the datasets whose genes are labelled by id (GSE81608: Entrez, GSE86469: Ensembl). One plain HTTPS request to the NCBI
Datasets gene report (api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/<symbols>/taxon/9606); nothing about the user is sent. Only exact symbol matches are kept; symbols that do not come back are listed.
Output: rp_ids_human.csv (symbol, entrez, ensembl); the attempt is logged in attempts.csv.   Usage: python -I fetch_rp_ids.py"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import pandas as pd
import wp4_http as H
from attempts_log import record
from rp_list import RP

HERE = os.path.dirname(os.path.abspath(__file__))
url = 'https://api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/' + '%2C'.join(RP) + '/taxon/9606?page_size=500'      # the default page holds only 20 reports
txt = H.get_text(url); rep = json.loads(txt).get('reports', [])
rows = {}
for r in rep:
    g = r.get('gene', {}); sym = g.get('symbol')
    if sym in RP and sym not in rows and g.get('tax_id') == '9606': rows[sym] = dict(symbol=sym, entrez=int(g['gene_id']), ensembl=';'.join(g.get('ensembl_gene_ids', [])))
D = pd.DataFrame([rows[s] for s in RP if s in rows]); D.to_csv(os.path.join(HERE, 'rp_ids_human.csv'), index=False)
missing = [s for s in RP if s not in rows]; noens = D[D.ensembl.str.count('ENSG') != 1].symbol.tolist()
print('requested', len(RP), '| returned with an exact symbol match', len(D), '| missing:', missing, '| without exactly one Ensembl id:', noens)
record('NCBI Gene', 'ids of the ribosomal protein genes (positive control A2)', url[:200], len(txt), 'NCBI Datasets gene report', 'downloaded', 'annotation only; keys for GSE81608 (Entrez) and GSE86469 (Ensembl)')
