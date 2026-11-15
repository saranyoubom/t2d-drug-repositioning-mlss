"""fetch_gse60803.py - protocol 4.6: download the series matrix of GSE60803 (1.4 MB) and the platform annotation GPL11532.annot.gz (5.6 MB; the 2 GB family file is NOT used) into the scratch data root
(v2_data/wp4/GSE60803/), one source per folder, sizes checked with HEAD, every attempt logged in attempts.csv. Downloaded files are untrusted data and are only read.
Usage: python -I fetch_gse60803.py"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import wp4_http as H
from attempts_log import record
import wp4_core as C

ROOT = os.path.join(os.path.dirname(C.PARSED), 'GSE60803'); G = 'https://ftp.ncbi.nlm.nih.gov/geo/'
for item, rel, fn in (('series matrix', 'series/GSE60nnn/GSE60803/matrix/GSE60803_series_matrix.txt.gz', 'GSE60803_series_matrix.txt.gz'), ('platform annotation GPL11532.annot.gz', 'platforms/GPL11nnn/GPL11532/annot/GPL11532.annot.gz', 'GPL11532.annot.gz')):
    url = G + rel; h = H.head(url)
    try:
        n = H.download(url, os.path.join(ROOT, fn), max_mb=300); record('GSE60803', item, url, n, 'GEO sample records (authors)', 'downloaded', 'size checked by HEAD before download'); print('ok', item, n)
    except Exception as e:
        record('GSE60803', item, url, h.get('size'), '', 'failed or refused', repr(e)[:300]); print('FAILED', item, repr(e)[:300])
