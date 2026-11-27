"""probe_pyc.py - is the cached bytecode of coexp_checks.py in the project folder the same program as its current source? Compares the code objects."""
import hashlib, marshal, os, sys, types
HERE = os.path.dirname(os.path.abspath(__file__)); P = os.path.normpath(os.path.join(HERE, '..', '..'))
src = os.path.join(P, 'coexp_checks.py'); pyc = os.path.join(P, '__pycache__', 'coexp_checks.cpython-314.pyc')
data = open(pyc, 'rb').read(); flags = int.from_bytes(data[4:8], 'little'); mtime = int.from_bytes(data[8:12], 'little'); size = int.from_bytes(data[12:16], 'little')
st = os.stat(src); text = open(src, 'rb').read()
print('source: size', st.st_size, 'mtime', int(st.st_mtime), '| pyc header: flags', flags, 'mtime', mtime, 'size', size, '| match:', int(st.st_mtime) == mtime and st.st_size == size)
print('source sha256', hashlib.sha256(text).hexdigest()[:16])
old = marshal.loads(data[16:]); print('pyc co_filename:', old.co_filename)
new = compile(text, old.co_filename, 'exec', dont_inherit=True, optimize=0)


def walk(a, b, path='<module>'):
    diffs = []
    if a.co_code != b.co_code: diffs.append(path + ': bytecode differs')
    ca, cb = a.co_consts, b.co_consts
    if len(ca) != len(cb): diffs.append(path + ': number of constants differs'); return diffs
    for x, y in zip(ca, cb):
        if isinstance(x, types.CodeType) and isinstance(y, types.CodeType): diffs += walk(x, y, path + '/' + x.co_name)
        elif x != y and not (x != x and y != y): diffs.append('%s: constant %r vs %r' % (path, str(x)[:50], str(y)[:50]))
    return diffs


d = walk(old, new); print('differences between the cached bytecode and the current source:', len(d))
for x in d[:15]: print('  ', x)
