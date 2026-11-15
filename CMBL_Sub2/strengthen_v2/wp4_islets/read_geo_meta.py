"""read_geo_meta.py - parse the sample-level header block of a GEO series matrix (lines starting with !Sample_) into a table; no expression values are read (the table section is skipped).
Functions: sample_table(path) -> DataFrame with one row per sample (columns: accession, title, source, ch1..chN split into key/value when written as 'key: value', supplementary files).
CLI: python -I read_geo_meta.py <series_matrix.txt.gz> [n_rows_to_print]  prints the characteristics and their value counts."""
import sys, gzip
import pandas as pd


def sample_table(path):
    fields = {}
    with gzip.open(path, 'rt', encoding='utf8', errors='replace') as f:
        for l in f:
            if l.startswith('!series_matrix_table_begin'): break
            if l.startswith('!Sample_'):
                k = l.split('\t')[0].strip(); vals = [x.strip().strip('"') for x in l.rstrip('\n').split('\t')[1:]]
                fields.setdefault(k, []).append(vals)
    cols = {}
    for k, lst in fields.items():
        for i, vals in enumerate(lst):
            name = k[len('!Sample_'):] + ('' if len(lst) == 1 else f'_{i + 1}')
            if k == '!Sample_characteristics_ch1':
                keys = {v.split(': ', 1)[0] for v in vals if ': ' in v}
                if len(keys) == 1:
                    name = 'ch:' + keys.pop(); vals = [v.split(': ', 1)[1] if ': ' in v else v for v in vals]
            cols[name] = vals
    return pd.DataFrame(cols)


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    T = sample_table(sys.argv[1]); nshow = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    print('samples', len(T)); print('columns', list(T.columns))
    print(T.head(nshow).T.to_string(max_colwidth=90))
    for c in T.columns:
        if c.startswith(('ch', 'source', 'molecule', 'extract', 'library', 'platform', 'description', 'data_processing')) or T[c].nunique() <= 12:
            vc = T[c].value_counts()
            if len(vc) <= 25: print('\n--', c, '\n', vc.to_string())
