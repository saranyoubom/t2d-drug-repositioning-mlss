"""make_numpy_reference.py - reference vectors from numpy for the R port of numpy's random generator (R_PORT_PLAN.md, gate G0).
Writes little-endian float64 files (integers are stored exactly as doubles) and manifest.json into ./reference/.
Run: python -I make_numpy_reference.py   (Python 3.14 with numpy; it reads and writes nothing else)."""
import json, os, random
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, 'reference')
os.makedirs(REF, exist_ok=True)


def split32(u64):
    u = np.asarray(u64, dtype=np.uint64)
    return (u & np.uint64(0xFFFFFFFF)).astype(np.float64), (u >> np.uint64(32)).astype(np.float64)


def limbs16(x, n):
    return [float((x >> (16 * i)) & 0xFFFF) for i in range(n)]


SEEDS = [
    {'name': 'int2026', 'kind': 'int', 'value': 2026, 'big': True},
    {'name': 'ss2026_5_1_7', 'kind': 'seedseq', 'value': [2026, 5, 1, 7], 'big': True},
    {'name': 'int0', 'kind': 'int', 'value': 0},
    {'name': 'int1', 'kind': 'int', 'value': 1},
    {'name': 'int7', 'kind': 'int', 'value': 7},
    {'name': 'list2031_1', 'kind': 'list', 'value': [2031, 1]},
    {'name': 'ss2029_3_2_1_9', 'kind': 'seedseq', 'value': [2029, 3, 2, 1, 9]},
    {'name': 'int5003', 'kind': 'int', 'value': 5003},
    {'name': 'int47bit', 'kind': 'int', 'value': 123456789012345},
]

OPS_SMALL = [
    {'op': 'random', 'n': 7}, {'op': 'integers', 'low': 0, 'high': 17, 'n': 5}, {'op': 'standard_normal', 'n': 9},
    {'op': 'permutation', 'n': 12}, {'op': 'choice_noreplace', 'pop': 20, 'k': 6}, {'op': 'random', 'n': 3},
    {'op': 'integers', 'low': 0, 'high': 400, 'n': 4}, {'op': 'normal', 'loc': 0.5, 'scale': 2.0, 'n': 5},
    {'op': 'permutation', 'n': 60}, {'op': 'integers_scalar', 'high': 2147483647}, {'op': 'choice_scalar', 'pop': 30},
    {'op': 'permutation_array', 'x': [3.5, 1.25, 9.0, 7.75]}, {'op': 'integers', 'low': 0, 'high': 2, 'n': 11},
    {'op': 'random_matrix', 'rows': 6, 'cols': 17}, {'op': 'choice_noreplace', 'pop': 3, 'k': 2},
    {'op': 'normal_array', 'loc': [0.0, 1.0, 2.0], 'scale': [1.0, 0.5, 2.0], 'rows': 5},
    {'op': 'choice_replace', 'pop': 25, 'n': 10}, {'op': 'integers', 'low': 0, 'high': 10, 'n': 20},
    {'op': 'permutation', 'n': 1000}, {'op': 'choice_noreplace', 'pop': 15000, 'k': 400},
    {'op': 'choice_noreplace', 'pop': 12000, 'k': 17}, {'op': 'standard_normal', 'n': 5000},
    {'op': 'random', 'n': 5001}, {'op': 'integers', 'low': 0, 'high': 3, 'n': 4000},
    {'op': 'integers', 'low': 0, 'high': 4294967295, 'n': 50}, {'op': 'standard_normal', 'n': 3}, {'op': 'random', 'n': 1},
]
OPS_BIG = [
    {'op': 'standard_normal', 'n': 200000}, {'op': 'random', 'n': 300001}, {'op': 'integers', 'low': 0, 'high': 3, 'n': 100000},
    {'op': 'permutation', 'n': 20000}, {'op': 'standard_normal', 'n': 7}, {'op': 'random', 'n': 2},
]


def make_gen(spec):
    if spec['kind'] == 'seedseq':
        return np.random.default_rng(np.random.SeedSequence(spec['value']))
    return np.random.default_rng(spec['value'])


def run_op(rng, op):
    k = op['op']
    if k == 'random':
        return rng.random(op['n'])
    if k == 'random_matrix':
        return rng.random((op['rows'], op['cols'])).ravel()
    if k == 'integers':
        return rng.integers(op['low'], op['high'], op['n']).astype(np.float64)
    if k == 'integers_scalar':
        return np.array([float(rng.integers(op['high']))])
    if k == 'standard_normal':
        return rng.standard_normal(op['n'])
    if k == 'normal':
        return rng.normal(op['loc'], op['scale'], op['n'])
    if k == 'normal_array':
        return rng.normal(np.array(op['loc']), np.array(op['scale']), (op['rows'], len(op['loc']))).ravel()
    if k == 'permutation':
        return rng.permutation(op['n']).astype(np.float64)
    if k == 'permutation_array':
        return rng.permutation(np.array(op['x'], dtype=np.float64))
    if k == 'choice_noreplace':
        return rng.choice(np.arange(op['pop']) * 10 + 1, op['k'], replace=False).astype(np.float64)
    if k == 'choice_scalar':
        return np.array([float(rng.choice(np.arange(op['pop']) * 3))])
    if k == 'choice_replace':
        return rng.choice(np.arange(op['pop']) * 3, op['n'], replace=True).astype(np.float64)
    raise ValueError(k)


manifest = {'numpy': np.__version__, 'python': os.sys.version.split()[0], 'seeds': [], 'stdlib': []}
for spec in SEEDS:
    name = spec['name']
    ss = np.random.SeedSequence(spec['value'])
    words = ss.generate_state(4, np.uint64)
    lo, hi = split32(words)
    bg = np.random.PCG64(np.random.SeedSequence(spec['value']))
    st = bg.state['state']
    head = list(ss.pool.astype(np.float64)) + [x for pair in zip(lo, hi) for x in pair] + limbs16(st['state'], 8) + limbs16(st['inc'], 8)
    raw = bg.random_raw(3000)
    rlo, rhi = split32(raw)
    rng = make_gen(spec)
    ops = OPS_SMALL + (OPS_BIG if spec.get('big') else [])
    outs, lens = [], []
    for op in ops:
        o = np.asarray(run_op(rng, op), dtype=np.float64).ravel()
        outs.append(o)
        lens.append(int(o.size))
    np.concatenate([np.array(head, dtype=np.float64), rlo, rhi, np.concatenate(outs)]).astype('<f8').tofile(os.path.join(REF, name + '.f64'))
    manifest['seeds'].append({'name': name, 'kind': spec['kind'], 'value': spec['value'], 'file': name + '.f64', 'n_head': len(head),
                              'n_raw': 3000, 'ops': ops, 'op_lengths': lens})

for sd, n in ((2026, 30), (2026, 200), (7, 17)):
    lst = list(range(n))
    r = random.Random(sd)
    r.shuffle(lst)
    r2 = random.Random(sd)
    manifest['stdlib'].append({'seed': sd, 'n': n, 'shuffled': lst, 'first_getrandbits32': [r2.getrandbits(32) for _ in range(5)]})

with open(os.path.join(REF, 'manifest.json'), 'w', encoding='utf-8') as fh:
    json.dump(manifest, fh)
print('numpy', np.__version__, '- wrote', len(manifest['seeds']), 'seed files and manifest.json to', REF)
for s in manifest['seeds']:
    print('  ', s['name'], s['n_head'] + 2 * s['n_raw'] + sum(s['op_lengths']), 'doubles')
