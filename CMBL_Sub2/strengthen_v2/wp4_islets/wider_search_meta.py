"""wider_search_meta.py - metadata of the candidates that passed the reading of titles and summaries (B2 of PRESPEC_WP4b...): sample titles and characteristics from the GEO series matrix headers (only lines starting
with '!' are read; the values table is never parsed), the supplementary file lists, and for ArrayExpress the SDRF and the study record. Used to write the contrast of each kept series BEFORE the freeze (B3).
Downloads go to <V2_DATA_ROOT>/v2_data/wp4/wider/<accession>/ (untrusted data, only read), each logged in attempts.csv.
Output: wider_search/sample_records_<accession>.tsv (one row per sample: title, characteristics) and a short printout.   Usage: python -I wider_search_meta.py ACC [ACC ...]"""
import os, sys, gzip, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding='utf8')
import pandas as pd
import wp4_http as H
import wp4_core as C
from attempts_log import record

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'wider_search'); os.makedirs(OUT, exist_ok=True)
ROOT = os.path.join(os.path.dirname(C.PARSED), 'wider'); G = 'https://ftp.ncbi.nlm.nih.gov/geo/series/'


def geo_dir(acc): n = int(acc[3:]); return f'{G}GSE{n // 1000}nnn/{acc}/'


def geo_meta(acc):
    d = os.path.join(ROOT, acc); os.makedirs(d, exist_ok=True); base = geo_dir(acc)
    for item, rel, fn in (('series matrix', f'matrix/{acc}_series_matrix.txt.gz', f'{acc}_series_matrix.txt.gz'), ('suppl filelist', 'suppl/filelist.txt', 'filelist.txt')):
        url = base + rel; h = H.head(url)
        try: n = H.download(url, os.path.join(d, fn), max_mb=300); record(acc, item, url, n, 'GEO sample records (authors)', 'downloaded (metadata read only)', 'size checked by HEAD; only lines starting with ! are read')
        except Exception as e: record(acc, item, url, h.get('size'), '', 'failed or refused', repr(e)[:200]); print(acc, item, 'FAILED', repr(e)[:120])
    rows = {}
    p = os.path.join(d, f'{acc}_series_matrix.txt.gz')
    if os.path.exists(p):
        with gzip.open(p, 'rt', encoding='utf8', errors='replace') as f:
            for l in f:
                if l.startswith('!series_matrix_table_begin'): break
                if l.startswith(('!Sample_title', '!Sample_geo_accession', '!Sample_source_name_ch1', '!Sample_characteristics_ch1', '!Sample_description', '!Sample_platform_id', '!Series_overall_design')):
                    k = l.split('\t')[0]; rows.setdefault(k, []).append([x.strip('"') for x in l.rstrip('\n').split('\t')[1:]])
    if '!Sample_title' not in rows: return None, rows
    n = len(rows['!Sample_title'][0]); T = pd.DataFrame({'gsm': rows.get('!Sample_geo_accession', [[''] * n])[0], 'title': rows['!Sample_title'][0], 'source': rows.get('!Sample_source_name_ch1', [[''] * n])[0]})
    for i, ch in enumerate(rows.get('!Sample_characteristics_ch1', [])): T[f'ch{i}'] = ch
    if '!Sample_platform_id' in rows: T['platform'] = rows['!Sample_platform_id'][0]
    fl = os.path.join(d, 'filelist.txt'); files = open(fl, encoding='utf8', errors='replace').read().splitlines()[:40] if os.path.exists(fl) else []
    return T, dict(design=rows.get('!Series_overall_design', [['']])[0][0][:600] if '!Series_overall_design' in rows else '', files=files)


def ae_meta(acc):
    d = os.path.join(ROOT, acc); os.makedirs(d, exist_ok=True); nn = acc.split('-')[-1]; base = f'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/{nn[-3:]}/{acc}/Files/'
    out = {}
    for fn in (f'{acc}.sdrf.txt', f'{acc}.idf.txt'):
        url = base + fn; h = H.head(url)
        try: n = H.download(url, os.path.join(d, fn), max_mb=50); record(acc, fn, url, n, 'ArrayExpress SDRF/IDF (authors)', 'downloaded (metadata)', ''); out[fn] = os.path.join(d, fn)
        except Exception as e: record(acc, fn, url, h.get('size'), '', 'failed or refused', repr(e)[:200]); print(acc, fn, 'FAILED', repr(e)[:120])
    s = json.loads(H.get_text('https://www.ebi.ac.uk/biostudies/api/v1/studies/' + acc)); files = []
    def walk(x):
        if isinstance(x, dict):
            if 'files' in x and isinstance(x['files'], list):
                for f in x['files']: files.append((f[0].get('path') if isinstance(f, list) and f else f.get('path'), f[0].get('size') if isinstance(f, list) and f else f.get('size')))
            for v in x.values(): walk(v)
        elif isinstance(x, list):
            for v in x: walk(v)
    walk(s); return out, files


if __name__ == '__main__':
    for acc in sys.argv[1:]:
        print('=' * 100); print(acc)
        if acc.startswith('GSE'):
            T, extra = geo_meta(acc)
            if T is None: print('no sample records'); continue
            T.to_csv(os.path.join(OUT, f'sample_records_{acc}.tsv'), sep='\t', index=False); pd.set_option('display.width', 280); pd.set_option('display.max_colwidth', 70)
            print('samples', len(T), '| platform', T.get('platform', pd.Series(['?'])).iloc[0], '| design:', extra['design'][:400]); print(T.drop(columns=[c for c in ('gsm', 'platform') if c in T.columns]).to_string(index=False, max_rows=40))
            print('suppl files:', extra['files'][:12])
        else:
            out, files = ae_meta(acc)
            if 'sdrf' in ' '.join(out):
                S = pd.read_csv(out[f'{acc}.sdrf.txt'], sep='\t', dtype=str); keep = [c for c in S.columns if c.startswith(('Source Name', 'Characteristics', 'Factor Value', 'Comment [ENA_RUN]', 'Derived Array Data File', 'Array Data File', 'Comment [Derived'))]
                pd.set_option('display.width', 280); pd.set_option('display.max_colwidth', 50); print(S[keep].to_string(index=False, max_rows=40)); S.to_csv(os.path.join(OUT, f'sample_records_{acc}.tsv'), sep='\t', index=False)
            print('files in the study record:', files[:30])
