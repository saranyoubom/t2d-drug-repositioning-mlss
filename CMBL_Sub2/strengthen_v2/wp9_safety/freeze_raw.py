"""freeze_raw.py - fingerprint of all saved openFDA responses (results/raw/*.json), written BEFORE the analysis scripts are run.
Recipe: sha256 over the text made of one line 'file name:sha256 of the file' for every file in results/raw, files sorted by name.
Writes results/RAW_SHA256.txt and appends the fingerprint (with the clock time) to ../PROTOCOL_v2_SHA256.txt.
Usage: python -I freeze_raw.py"""
import os, hashlib, datetime, json
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, 'results', 'raw')
files = sorted(f for f in os.listdir(RAW) if f.endswith('.json'))
assert not [f for f in os.listdir(RAW) if not f.endswith('.json')], 'unexpected non-json file in raw/'
lines = ''.join('%s:%s\n' % (f, hashlib.sha256(open(os.path.join(RAW, f), 'rb').read()).hexdigest()) for f in files)
h = hashlib.sha256(lines.encode()).hexdigest()
times = []
lu = set()
for f in files:
    rec = json.load(open(os.path.join(RAW, f), encoding='utf8'))
    times.append(rec['retrieved'])
    lu.add(rec['response'].get('meta', {}).get('last_updated'))
stamp = datetime.datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S %z')
text = '%s *wp9_safety/results/raw (%d files)\nfirst retrieval %s, last retrieval %s, openFDA last_updated %s; written %s BEFORE the analysis\n' % (h, len(files), min(times), max(times), ','.join(sorted(x for x in lu if x)), stamp)
open(os.path.join(HERE, 'results', 'RAW_SHA256.txt'), 'w', encoding='utf8').write(text)
p = os.path.join(os.path.dirname(HERE), 'PROTOCOL_v2_SHA256.txt')
b = open(p, 'rb').read()
if not b.endswith(b'\n'):
    b += b'\n'
b += ('%s *wp9_safety/results/RAW_SHA256.txt\nWP9 raw openFDA responses (%d files) fingerprinted BEFORE the analysis: %s\n' % (h, len(files), stamp)).encode('utf8')
open(p, 'wb').write(b)
print(text)
