"""bench_data_io.py - importable loaders of benchmark part 2 (WP2, protocol v2 section 3), with no side effects at import.
Copies of the loaders of ../../benchmark_route3/run_benchmark.py (that script reads sys.argv at import and cannot be imported; the original is not modified): fetch, series_path, read_matrix,
annotation, gene_matrix, anova_fdr, and the code path that turns one frozen series into the two 12 x 17 analysis matrices of the protocol (gene set A = 17 random differentially expressed
genes, gene set B = 17 random expressed genes). Development series, held-out arrays and RNA-seq series use the SAME functions (prepare_from_log -> prepare_array / prepare_rnaseq).
The logic of the functions is the part 1 logic line for line, so that the seeds give the same draws (checked by gate_part1.py).
Downloaded files are untrusted data: they are read as text or numbers only, never executed; run Python with -I. NCBI politeness: throttle() keeps requests at most 3 per second.
Usage: import bench_data_io as io_; io_.configure(matrices=..., platforms=...)."""
import os, re, io, gzip, time, urllib.request, urllib.error
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

FTP = 'https://ftp.ncbi.nlm.nih.gov/geo/'
# scratch data root outside the project (downloads are untrusted data, kept apart per source); override with the environment variable V2_DATA_ROOT
DEFAULT_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_cache')
NPERM = 2000; MIN_DE = 34; NGENE = 17
N, K = 12, 3
GROUP = np.repeat(np.arange(K), 4)                                        # 3 groups x 4 replicates in file order (same as bench_pipeline.GROUP)
MATRICES = None; PLATFORMS = None; ANNOT = {}
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ncbi_http import throttle, http_get                                    # at most 3 NCBI requests per second (ncbi_http.py)


def data_root(): return os.environ.get('V2_DATA_ROOT', DEFAULT_ROOT)


def configure(matrices, platforms):
    """folders (created if missing) for the series matrices and the platform annotations"""
    global MATRICES, PLATFORMS
    MATRICES, PLATFORMS = matrices, platforms
    for d in (matrices, platforms): os.makedirs(d, exist_ok=True)


def fetch(url, path):
    """download to path (atomic .part rename) unless it exists; copy of run_benchmark.fetch with the throttle"""
    if os.path.exists(path) and os.path.getsize(path) > 0: return path
    for k in range(4):
        try:
            throttle()
            with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'geo-benchmark/1.0'}), timeout=180) as r, open(path + '.part', 'wb') as f: f.write(r.read())
            os.replace(path + '.part', path); return path
        except urllib.error.HTTPError as e:
            if e.code == 404: raise
            time.sleep(3 + 3 * k)
        except Exception: time.sleep(3 + 3 * k)
    raise RuntimeError('download failed ' + url)


def series_path(gse, gpl):
    n = int(gse[3:]); folder = f'GSE{n // 1000}nnn' if n >= 1000 else 'GSEnnn'; base = f'{FTP}series/{folder}/{gse}/matrix/'
    path = os.path.join(MATRICES, f'{gse}_{gpl}.txt.gz')
    try: return fetch(base + f'{gse}_series_matrix.txt.gz', path)
    except urllib.error.HTTPError:
        throttle(); idx = urllib.request.urlopen(base, timeout=60).read().decode('utf8', 'replace'); hit = [x for x in re.findall(r'href="([^"]*series_matrix\.txt\.gz)"', idx) if gpl in x]
        return fetch(base + hit[0], path)


def _table_bounds(path, begin, end):
    """0-based line numbers of the begin and end marker lines of a table inside a gzip text file, found by streaming (the text is not kept in memory)"""
    b = e = None
    with gzip.open(path, 'rt', encoding='utf8', errors='replace') as f:
        for i, l in enumerate(f):
            if b is None:
                if l.startswith(begin): b = i
            elif l.startswith(end): e = i; break
    if b is None or e is None: raise ValueError('table markers not found in ' + os.path.basename(path))
    return b, e


def read_matrix(path):
    """the expression table of a series matrix. LEAN VERSION (5 Nov 2026, memory): same table as the original, which held the whole text, its line list and two copies in memory (about 4 times the text);
    this one finds the marker lines by streaming and lets the pandas C parser read the table from the stream. Equality with the original is tested on the cached series (gate re-run)."""
    b, e = _table_bounds(path, '!series_matrix_table_begin', '!series_matrix_table_end')
    with gzip.open(path, 'rt', encoding='utf8', errors='replace') as f:
        return pd.read_csv(f, sep='\t', index_col=0, quotechar='"', low_memory=False, skiprows=b + 1, nrows=e - b - 2)


def annotation(gpl):
    n = int(gpl[3:]); folder = f'GPL{n // 1000}nnn' if n >= 1000 else 'GPLnnn'
    path = fetch(f'{FTP}platforms/{folder}/{gpl}/annot/{gpl}.annot.gz', os.path.join(PLATFORMS, f'{gpl}.annot.gz'))
    b, e = _table_bounds(path, '!platform_table_begin', '!platform_table_end')                       # LEAN VERSION (5 Nov 2026): only the ID and symbol columns are read, from the stream
    with gzip.open(path, 'rt', encoding='utf8', errors='replace') as f: names = list(pd.read_csv(f, sep='\t', quotechar='"', low_memory=False, dtype=str, skiprows=b + 1, nrows=0).columns)
    col = next(c for c in names if 'symbol' in c.lower())
    with gzip.open(path, 'rt', encoding='utf8', errors='replace') as f: t = pd.read_csv(f, sep='\t', quotechar='"', low_memory=False, dtype=str, skiprows=b + 1, nrows=e - b - 2, usecols=['ID', col])
    out = t[['ID', col]].dropna(); out.columns = ['probe', 'symbol']
    out['symbol'] = out.symbol.str.split(' /// ').str[0].str.strip(); return out[(out.symbol != '') & (out.symbol != '---')].set_index('probe').symbol


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
    gm = np.vstack([L[GROUP == g].mean(0) for g in range(K)]); grand = L.mean(0)
    ssb = 4 * ((gm - grand) ** 2).sum(0); ssw = ((L - gm[GROUP]) ** 2).sum(0)
    with np.errstate(divide='ignore', invalid='ignore'): F = (ssb / (K - 1)) / (ssw / (N - K))
    p = stats.f.sf(F, K - 1, N - K); p = np.nan_to_num(p, nan=1.0); return p, multipletests(p, method='fdr_bh')[1]


def prepare_from_log(rec, L2, genes, n, xfun, dry=False):
    """shared code path of part 1 (run_benchmark.process) after the log-scale gene matrix exists: expressed genes (mean log expression above the median), differential expression (ANOVA, BH
    FDR < 0.05), analysability (>= MIN_DE DE genes) and the two random draws of 17 genes (seed 2026 + GSE number; set B seed [2026 + n, 1]). xfun(draw) gives the 12 x 17 linear matrix.
    Returns (rec, sets) with sets = {'A': dict(X, genes, n_de_in_set, idx), 'B': ...} or None when not analysable / dry run."""
    expressed = L2.mean(0) > np.median(L2.mean(0)); p, fdr = anova_fdr(L2); de = np.where(expressed & (fdr < 0.05))[0]; rec['n_expressed'] = int(expressed.sum()); rec['n_de'] = int(len(de))
    if len(de) < MIN_DE: rec['status'] = 'not analysable: %d DE genes (< %d)' % (len(de), MIN_DE); return rec, None
    if dry: rec['status'] = 'analysable (dry run)'; return rec, None
    order = np.argsort(np.array(genes)[de], kind='stable'); de_sorted = de[order]
    drawA = np.random.default_rng(2026 + n).choice(de_sorted, NGENE, replace=False)
    ex_sorted = np.where(expressed)[0][np.argsort(np.array(genes)[np.where(expressed)[0]], kind='stable')]
    drawB = np.random.default_rng([2026 + n, 1]).choice(ex_sorted, NGENE, replace=False)
    sets = {}
    for tag, draw in (('A', drawA), ('B', drawB)):
        sets[tag] = dict(X=np.asarray(xfun(draw), float), genes=[genes[i] for i in draw], n_de_in_set=int((fdr[draw] < 0.05).sum()), idx=np.asarray(draw))
    rec['status'] = 'prepared'; return rec, sets


def prepare_array(row, dry=False):
    """one frozen microarray series (row of a frozen csv) -> (rec, sets); same record fields as part 1"""
    gse, gpl = row['gse'], row['platform']; n = int(gse[3:]); gsms = row['gsms'].split(';')
    rec = dict(gse=gse, platform=gpl, order=int(row['order']), how=row['how'], groups=row['groups'])
    L2, genes, is_log, q99 = gene_matrix(gse, gpl, gsms); rec.update(n_genes_mapped=len(genes), scale='log2' if is_log else 'linear', q99=q99)
    return prepare_from_log(rec, L2, genes, n, lambda d: 2.0 ** L2[:, d], dry)


# ------------------------------------------------------------------------------------------------------------------ RNA-seq (protocol 3.5)
def rnaseq_cpm(counts):
    """counts: genes x 12 DataFrame/array of raw counts for the selected samples. CPM on all genes of the file; universe = genes with count >= 1 in all 12 samples.
    Returns (cpm of the universe genes, 12 x U; log2(cpm + 1), 12 x U; universe row positions)"""
    C = np.asarray(counts, float); lib = C.sum(axis=0); cpm = C / lib[None, :] * 1e6; keep = np.where((C >= 1).all(axis=1))[0]
    return cpm[keep].T, np.log2(cpm[keep].T + 1.0), keep


def prepare_rnaseq(rec, counts, ids, symbols, n, dry=False):
    """RNA-seq series: counts (genes x 12), gene ids (list) and symbols (list, same order; fall back to ids). Linear scale = CPM on all genes; genes need count >= 1 in all 12 samples;
    expressed and DE as in part 1 on log2(CPM + 1); X for the analysis = CPM of the drawn genes"""
    cpm, L2, keep = rnaseq_cpm(counts); genes = [str(symbols[i]) if symbols[i] not in (None, '') and not (isinstance(symbols[i], float) and np.isnan(symbols[i])) else str(ids[i]) for i in keep]
    rec.update(n_genes_mapped=len(genes), scale='counts->CPM', n_genes_total=int(np.asarray(counts).shape[0]), min_lib=float(np.asarray(counts, float).sum(axis=0).min()))
    return prepare_from_log(rec, L2, genes, n, lambda d: cpm[:, d], dry)
