"""make_ziggurat_tables.py - writes ziggurat_nor.f64: ki_double (256), wi_double (256), fi_double (256), ziggurat_nor_r, ziggurat_nor_inv_r
as little-endian float64, from numpy's ziggurat_constants.h. These are the numeric constants of numpy's normal generator (numpy is
BSD-3-Clause); the header is read from GitHub into memory (it is not saved), and only the numbers are written. Needed once: numpy_rng.R
reads the small binary file. Run: python -I make_ziggurat_tables.py"""
import os, re, sys, urllib.request
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TAGS = ['v' + np.__version__, 'main']
txt = None
for tag in TAGS:
    url = f'https://raw.githubusercontent.com/numpy/numpy/{tag}/numpy/random/src/distributions/ziggurat_constants.h'
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            txt = r.read().decode('utf-8')
        print('read', url, len(txt), 'characters')
        break
    except Exception as e:
        print('could not read', url, '-', e)
if txt is None:
    sys.exit('no header read')


def table(name):
    m = re.search(r'\b' + name + r'\s*\[\s*\]\s*=\s*\{(.*?)\};', txt, re.S)
    if not m:
        sys.exit('array not found: ' + name)
    out = []
    for tok in m.group(1).split(','):
        tok = re.sub(r'/\*.*?\*/', '', tok, flags=re.S).strip()
        if not tok:
            continue
        tok = re.sub(r'(?i)(ull|ul|u|l)$', '', tok)
        out.append(float(int(tok, 16)) if tok.lower().startswith('0x') else float(tok))
    return out


def const(name):
    m = re.search(r'(?im)^[^\n]*\b' + name + r'\b[ \t]*(?:=[ \t]*)?([-+]?[0-9][-+0-9.eE]*)', txt)
    if not m:
        sys.exit('constant not found: ' + name)
    print('   ', name, '<-', m.group(0).strip()[:100])
    return float(m.group(1))


ki, wi, fi = table('ki_double'), table('wi_double'), table('fi_double')
r = const('ZIGGURAT_NOR_R')
try:
    inv_r = const('ZIGGURAT_NOR_INV_R')
except SystemExit:
    print('    ZIGGURAT_NOR_INV_R is not in this header; showing the lines that mention "inv" or "nor_":')
    for ln in txt.splitlines():
        if re.search(r'(?i)inv|nor_', ln) and len(ln) < 160:
            print('       ', ln.strip())
    inv_r = float('0.27366123732975827203338247596')     # ZIGGURAT_NOR_INV_R of numpy's distributions.c (= 1 / ziggurat_nor_r)
    print('    using the literal', inv_r, ' | 1/r =', 1.0 / r, ' | difference', inv_r - 1.0 / r)
assert len(ki) == len(wi) == len(fi) == 256, (len(ki), len(wi), len(fi))
print('ki[0..2]', ki[:3], ' wi[0..2]', wi[:3], ' fi[0..2]', fi[:3], ' r', r, ' inv_r', inv_r)
assert all(k == int(k) and k < 2 ** 53 for k in ki)
np.array(ki + wi + fi + [r, inv_r], dtype='<f8').tofile(os.path.join(HERE, 'ziggurat_nor.f64'))
print('wrote ziggurat_nor.f64', 770 * 8, 'bytes')
