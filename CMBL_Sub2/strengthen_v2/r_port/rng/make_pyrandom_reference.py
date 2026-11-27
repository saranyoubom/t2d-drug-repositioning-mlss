"""Reference values of Python's random.Random for the R port rng/pyrandom.R (gate G0b): getrandbits(32) streams for two seeds (one wider than 32 bits), getrandbits(7), randbelow-driven shuffles of 100 and
5000 elements. Writes pyrandom_reference.csv (columns: case, index, value). Usage: python -I make_pyrandom_reference.py <out.csv>"""
import sys, random
rows = []
r = random.Random(2026); rows += [('bits32_seed2026', i, r.getrandbits(32)) for i in range(2000)]
r = random.Random(2 ** 40 + 5); rows += [('bits32_seed2p40p5', i, r.getrandbits(32)) for i in range(700)]
r = random.Random(2026); rows += [('bits7_seed2026', i, r.getrandbits(7)) for i in range(60)]
r = random.Random(2026); x = list(range(100)); r.shuffle(x); rows += [('shuffle100_seed2026', i, v) for i, v in enumerate(x)]
r = random.Random(2026); x = list(range(5000)); r.shuffle(x); rows += [('shuffle5000_seed2026', i, v) for i, v in enumerate(x)]
r = random.Random(7); x = list(range(37)); r.shuffle(x); rows += [('shuffle37_seed7', i, v) for i, v in enumerate(x)]
with open(sys.argv[1], 'w', encoding='utf8', newline='') as f:
    f.write('case,index,value\n')
    for c, i, v in rows: f.write('%s,%d,%d\n' % (c, i, v))
print('wrote', len(rows), 'values to', sys.argv[1])
