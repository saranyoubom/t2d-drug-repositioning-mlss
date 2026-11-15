"""Run the pre-specified benchmark (PROTOCOL sections 2-5) on the frozen dataset list. One JSON line per processed series in results_per_dataset.jsonl; every series, analysable or not, is logged with its status.
Series are processed in the frozen order; a series that is not analysable (fewer than 34 DE genes, unreadable file, fewer than 12 usable samples, platform not annotated) is logged and the next eligible series is used, until 20 series per platform are analysable or the list is exhausted.
Downloaded files are untrusted data (kept in <data dir>, read only as text/numbers). Usage: python -I run_benchmark.py <data dir> [max_per_platform=20] [dry]
'dry' only decides analysability (DE gene count and file integrity, no correlation is computed) and writes dry_run_per_dataset.jsonl; it lets the list be extended blind to any benchmark result."""
import sys, os, re, io, gzip, csv, json, time, urllib.request, urllib.error
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
sys.stdout.reconfigure(encoding='utf8')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bench_pipeline as B
DATA = sys.argv[1]; MAXP = int(sys.argv[2]) if len(sys.argv) > 2 else 20
for d in ('matrices', 'platforms'): os.makedirs(os.path.join(DATA, d), exist_ok=True)
FTP = 'https://ftp.ncbi.nlm.nih.gov/geo/'
NPERM = 2000; MIN_DE = 34; NGENE = 17; DRY = len(sys.argv) > 3 and sys.argv[3] == 'dry'


def fetch(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 0: return path
    for k in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'geo-benchmark/1.0'}), timeout=180) as r, open(path + '.part', 'wb') as f: f.write(r.read())
            os.replace(path + '.part', path); return path
        except urllib.error.HTTPError as e:
            if e.code == 404: raise
            time.sleep(3 + 3 * k)
        except Exception: time.sleep(3 + 3 * k)
    raise RuntimeError('download failed ' + url)


def series_path(gse, gpl):
    n = int(gse[3:]); folder = f'GSE{n // 1000}nnn' if n >= 1000 else 'GSEnnn'; base = f'{FTP}series/{folder}/{gse}/matrix/'
    path = os.path.join(DATA, 'matrices', f'{gse}_{gpl}.txt.gz')
    try: return fetch(base + f'{gse}_series_matrix.txt.gz', path)
    except urllib.error.HTTPError:
        idx = urllib.request.urlopen(base, timeout=60).read().decode('utf8', 'replace'); hit = [x for x in re.findall(r'href="([^"]*series_matrix\.txt\.gz)"', idx) if gpl in x]
        return fetch(base + hit[0], path)


def read_matrix(path):
    with gzip.open(path, 'rt', encoding='utf8', errors='replace') as f: lines = f.read().split('\n')
    b = next(i for i, l in enumerate(lines) if l.startswith('!series_matrix_table_begin')); e = next(i for i, l in enumerate(lines) if l.startswith('!series_matrix_table_end'))
    return pd.read_csv(io.StringIO('\n'.join(lines[b + 1:e])), sep='\t', index_col=0, quotechar='"', low_memory=False)


def annotation(gpl):
    n = int(gpl[3:]); folder = f'GPL{n // 1000}nnn' if n >= 1000 else 'GPLnnn'
    path = fetch(f'{FTP}platforms/{folder}/{gpl}/annot/{gpl}.annot.gz', os.path.join(DATA, 'platforms', f'{gpl}.annot.gz'))
    with gzip.open(path, 'rt', encoding='utf8', errors='replace') as f: lines = f.read().split('\n')
    b = next(i for i, l in enumerate(lines) if l.startswith('!platform_table_begin')); e = next(i for i, l in enumerate(lines) if l.startswith('!platform_table_end'))
    t = pd.read_csv(io.StringIO('\n'.join(lines[b + 1:e])), sep='\t', quotechar='"', low_memory=False, dtype=str)
    col = next(c for c in t.columns if 'symbol' in c.lower()); out = t[['ID', col]].dropna(); out.columns = ['probe', 'symbol']
    out['symbol'] = out.symbol.str.split(' /// ').str[0].str.strip(); return out[(out.symbol != '') & (out.symbol != '---')].set_index('probe').symbol


ANNOT = {}


def gene_matrix(gse, gpl, gsms):
    """12 samples x genes on the log2 scale (probe with the highest mean linear expression per gene) and the scale decision"""
    M = read_matrix(series_path(gse, gpl)); miss = [g for g in gsms if g not in M.columns]
    if miss: raise ValueError('samples missing: %s' % miss[:3])
    M = M[gsms].apply(pd.to_numeric, errors='coerce')
    q99 = float(np.nanpercentile(M.values, 99)); is_log = q99 < 30; lin = 2.0 ** M if is_log else M
    ok = (lin > 0).all(axis=1) & lin.notna().all(axis=1); lin = lin[ok]
    if gpl not in ANNOT: ANNOT[gpl] = annotation(gpl)
    sym = ANNOT[gpl].reindex(lin.index); lin = lin[sym.notna().values]; sym = sym.dropna()
    mean = lin.mean(axis=1); best = mean.groupby(sym.values).idxmax()               # one probe per gene symbol
    G = lin.loc[best.values]; G.index = best.index
    return np.log2(G.T.values), list(G.index), is_log, q99


def anova_fdr(L):
    """one-way ANOVA across the 3 groups for every gene (columns of L = genes), BH FDR"""
    gm = np.vstack([L[B.GROUP == g].mean(0) for g in range(B.K)]); grand = L.mean(0)
    ssb = 4 * ((gm - grand) ** 2).sum(0); ssw = ((L - gm[B.GROUP]) ** 2).sum(0)
    with np.errstate(divide='ignore', invalid='ignore'): F = (ssb / (B.K - 1)) / (ssw / (B.N - B.K))
    p = stats.f.sf(F, B.K - 1, B.N - B.K); p = np.nan_to_num(p, nan=1.0); return p, multipletests(p, method='fdr_bh')[1]


def run_set(X, den, seed):
    return B.analyse(X, den, nperm=NPERM, seed=seed)


def spike(Llog, den, rng):
    """positive control: shared latent factor added to 6 of the 14 non-denominator genes (log scale)"""
    cand = [i for i in range(Llog.shape[1]) if i not in den]; pick = rng.choice(cand, 6, replace=False)
    gm = np.vstack([Llog[B.GROUP == g].mean(0) for g in range(B.K)]); sd = np.sqrt(((Llog - gm[B.GROUP]) ** 2).sum(0) / (B.N - B.K)); z = rng.standard_normal(B.N)
    out = Llog.copy(); out[:, pick] += 1.5 * z[:, None] * sd[pick][None, :]; return out, pick


def process(row):
    gse, gpl = row['gse'], row['platform']; n = int(gse[3:]); gsms = row['gsms'].split(';'); rec = dict(gse=gse, platform=gpl, order=int(row['order']), how=row['how'], groups=row['groups'])
    L2, genes, is_log, q99 = gene_matrix(gse, gpl, gsms); rec.update(n_genes_mapped=len(genes), scale='log2' if is_log else 'linear', q99=q99)
    expressed = L2.mean(0) > np.median(L2.mean(0)); p, fdr = anova_fdr(L2); de = np.where(expressed & (fdr < 0.05))[0]; rec['n_expressed'] = int(expressed.sum()); rec['n_de'] = int(len(de))
    if len(de) < MIN_DE: rec['status'] = 'not analysable: %d DE genes (< %d)' % (len(de), MIN_DE); return rec
    if DRY: rec['status'] = 'analysable (dry run)'; return rec
    order = np.argsort(np.array(genes)[de], kind='stable'); de_sorted = de[order]
    drawA = np.random.default_rng(2026 + n).choice(de_sorted, NGENE, replace=False)
    ex_sorted = np.where(expressed)[0][np.argsort(np.array(genes)[np.where(expressed)[0]], kind='stable')]
    drawB = np.random.default_rng([2026 + n, 1]).choice(ex_sorted, NGENE, replace=False)
    for tag, draw in (('A', drawA), ('B', drawB)):
        X = 2.0 ** L2[:, draw]; den = [0, 1, 2]; res = run_set(X, den, 2026 + n); res['genes'] = [genes[i] for i in draw]; res['n_de_in_set'] = int((fdr[draw] < 0.05).sum()); rec[tag] = res
    LA = np.log(2.0 ** L2[:, drawA]); sp, pick = spike(LA, [0, 1, 2], np.random.default_rng([2026 + n, 2])); res = run_set(np.exp(sp), [0, 1, 2], 2026 + n)
    res['spiked_genes'] = [genes[drawA[i]] for i in pick]; rec['A_spike'] = res; rec['status'] = 'analysed'; return rec


def main():
    rows = list(csv.DictReader(open(os.path.join(HERE, 'frozen_datasets.csv'), encoding='utf8'))); done = {g: 0 for g in ('GPL570', 'GPL1261')}; out = open(os.path.join(HERE, 'dry_run_per_dataset.jsonl' if DRY else 'results_per_dataset.jsonl'), 'w', encoding='utf8')
    used = {}                                                       # selected sample -> series, for the duplicate rule (protocol deviation D1)
    for row in rows:
        gpl = row['platform']
        if done[gpl] >= MAXP: continue
        t0 = time.time(); gs = set(row['gsms'].split(';')); clash = sorted({used[g] for g in gs if g in used})
        if clash: rec = dict(gse=row['gse'], platform=gpl, order=int(row['order']), status='not analysable: duplicate of ' + clash[0])
        else:
            try: rec = process(row)
            except Exception as e: rec = dict(gse=row['gse'], platform=gpl, order=int(row['order']), status='not analysable: %s %s' % (type(e).__name__, str(e)[:120]))
        if rec['status'] in ('analysed', 'analysable (dry run)'):
            done[gpl] += 1; used.update({g: row['gse'] for g in gs})
        out.write(json.dumps(rec, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n'); out.flush()
        print(f"{gpl} {row['gse']} (#{row['order']}): {rec['status']} | DE {rec.get('n_de', '-')} | done {done} | {time.time() - t0:.0f}s", flush=True)
    print('finished', done, flush=True)


if __name__ == '__main__':
    main()
