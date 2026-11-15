"""inspect_gz.py - streaming peek at a (possibly large) gzip text table without loading it: number of lines, number of columns of the first lines, first column labels,
and the first characters of the first lines. Used to check file shapes before any parsing plan is made (8 GB machine: nothing is read densely).
Usage: python -I inspect_gz.py <file.gz|file.txt> [n_lines_to_show] [sep]"""
import sys, gzip, io
sys.stdout.reconfigure(encoding='utf8')
path = sys.argv[1]; show = int(sys.argv[2]) if len(sys.argv) > 2 else 3; sep = sys.argv[3] if len(sys.argv) > 3 else None
opener = gzip.open if path.endswith('.gz') else open
n = 0; first = []
with opener(path, 'rt', encoding='utf8', errors='replace', newline=None) as f:
    for line in f:
        n += 1
        if n <= show: first.append(line.rstrip('\n'))
print('file', path, '| lines', n)
for i, l in enumerate(first, 1):
    s = sep if sep else ('\t' if '\t' in l else ',')
    parts = l.split(s)
    print(f'line {i}: {len(parts)} fields (sep={s!r}) | first 6: {parts[:6]} | last 3: {parts[-3:]} | chars {len(l)}')
