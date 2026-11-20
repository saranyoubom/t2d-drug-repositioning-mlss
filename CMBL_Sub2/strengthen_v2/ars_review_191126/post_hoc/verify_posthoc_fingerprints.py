"""Recomputes the SHA-256 of the two method notes of this folder and compares them with FINGERPRINTS.txt (MATCH / CHANGED). Usage: python -I verify_posthoc_fingerprints.py"""
import hashlib, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ok = True
for line in open(os.path.join(HERE, 'FINGERPRINTS.txt'), encoding='utf-8'):
    m = re.match(r'^([0-9a-f]{64})\s+(\S+)', line)
    if m:
        h = hashlib.sha256(open(os.path.join(HERE, m.group(2)), 'rb').read()).hexdigest(); print('MATCH' if h == m.group(1) else 'CHANGED', m.group(2)); ok &= (h == m.group(1))
sys.exit(0 if ok else 1)
