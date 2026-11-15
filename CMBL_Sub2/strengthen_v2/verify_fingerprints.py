"""verify_fingerprints.py - recompute the SHA-256 of the files named in the fingerprint lists of this folder (release copy, November 2026).
Usage (python -I, from this folder): python -I verify_fingerprints.py
For every line "<sha256> *<path>" of PROTOCOL_v2_SHA256.txt and of the per-folder lists, the file is looked up relative to this folder, to its parent CMBL_Sub2, to the folder of the list and to the repository root.
Outcomes: MATCH (the file equals the last recorded hash of its path, or one of them), CHANGED (the file exists but equals none of the recorded hashes: an earlier version was recorded and the file was edited later,
as the protocol file was after each entry of its deviation log), NOT SHIPPED (input data that this repository does not contain, for example raw instrument exports and database files).
The raw openFDA responses are fingerprinted as a set; unpack results/raw_responses.zip into results/raw/ first (see results/README_raw_responses.md) to check them."""
import sys, os, re, glob, hashlib
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.abspath(__file__)); CMBL = os.path.dirname(HERE); ROOT = os.path.dirname(CMBL)
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
LINE = re.compile(r'^\s*([0-9a-fA-F]{64})\s+\*?(.+?)\s*$')
lists = [os.path.join(HERE, 'PROTOCOL_v2_SHA256.txt')] + [p for pat in ('*/*FINGERPRINT*.txt', '*/*SHA256*.txt', '*/*/*SHA256*.txt') for p in glob.glob(os.path.join(HERE, pat))]
lists += glob.glob(os.path.join(CMBL, 'benchmark_route3', '*SHA256*.txt'))
rec = {}; notshipped = {}; rawsets = []
for lf in dict.fromkeys(lists):
    if not os.path.exists(lf): continue
    for line in open(lf, encoding='utf-8', errors='replace'):
        m = LINE.match(line)
        if not m: continue
        h, rel = m.group(1).lower(), m.group(2).strip().strip('"')
        if re.search(r'\(\d+ files\)', rel) or rel.endswith('RAW_SHA256.txt'): rawsets.append((os.path.relpath(lf, HERE), h, rel)); continue
        relp = rel.replace('/', os.sep)
        for b in (HERE, CMBL, os.path.dirname(lf), ROOT):
            c = os.path.normpath(os.path.join(b, relp))
            if os.path.isfile(c): rec.setdefault(c, set()).add(h); break
        else: notshipped.setdefault(rel, set()).add(h)
def state(p, hs):
    b = open(p, 'rb').read()
    if hashlib.sha256(b).hexdigest() in hs: return 'MATCH'
    if hashlib.sha256(b.replace(b'\r\n', b'\n').replace(b'\n', b'\r\n')).hexdigest() in hs: return 'CRLF'
    return 'CHANGED'
st = {p: state(p, hs) for p, hs in rec.items()}
match = [p for p in st if st[p] == 'MATCH']; crlf = [p for p in st if st[p] == 'CRLF']; changed = [p for p in st if st[p] == 'CHANGED']
print('files with a recorded fingerprint in this repository: %d' % len(rec)); print('  MATCH   %d' % len(match))
print('  MATCH after restoring Windows line endings (CRLF)   %d  (files of older folders that git stored with LF)' % len(crlf))
print('  CHANGED %d' % len(changed))
for p in sorted(changed): print('     CHANGED', os.path.relpath(p, HERE), '(recorded hashes: %d)' % len(rec[p]))
print('NOT SHIPPED (inputs that this repository does not contain): %d' % len(notshipped))
for rel in sorted(notshipped)[:60]: print('     ', rel)
seen_raw = set()
for lf, h, rel in rawsets:
    wp = 'wp9_safety' if 'wp9' in rel else 'wp12_database_annotations' if 'wp12' in rel else None
    if wp is None or (wp, h) in seen_raw: continue
    seen_raw.add((wp, h)); raw = os.path.join(HERE, wp, 'results', 'raw')
    if os.path.isdir(raw):
        names = sorted(os.listdir(raw)); r = hashlib.sha256(''.join('%s:%s\n' % (n, sha(os.path.join(raw, n))) for n in names).encode()).hexdigest()
        print('raw openFDA responses of %s (%d files): recipe hash %s, recorded %s -> %s' % (wp, len(names), r[:16], h[:16], 'MATCH' if r == h else 'DIFFERENT'))
    else: print('raw openFDA responses of %s: not unpacked (results/raw_responses.zip), not checked' % wp)
