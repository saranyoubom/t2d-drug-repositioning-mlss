"""probe_sources.py - read-only reconnaissance of the public download locations of the four additional human datasets of protocol 4.3 (GSE86469, GSE81547, E-MTAB-5061, GSE85241):
lists the GEO ftp folders (suppl/ and matrix/) and HEADs every file for its size. Nothing is downloaded. Output is printed (and tee'd into logs/probe_sources.txt by the caller).
Usage: python -I probe_sources.py"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import wp4_http as H

GEO = 'https://ftp.ncbi.nlm.nih.gov/geo/series/'
for gse in ('GSE86469', 'GSE81547', 'GSE85241'):
    n = int(gse[3:]); folder = f'GSE{n // 1000}nnn'
    for sub in ('matrix/', 'suppl/'):
        url = f'{GEO}{folder}/{gse}/{sub}'
        try:
            items = H.list_index(url)
        except Exception as e:
            print(gse, sub, 'LISTING FAILED', repr(e)); continue
        for name, href in items:
            if name.lower().startswith('parent'): continue
            h = H.head(url + href)
            print(f'{gse}\t{sub}{href}\t{h["status"]}\t{h["size"]}\t{h["type"]}')
