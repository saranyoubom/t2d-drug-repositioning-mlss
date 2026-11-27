"""make_pyformat_reference.py - reference values for test_pyformat.R: Python's float(), repr(), '%.16g', '%.6g', '%.2f' and '%.6f' for hard and random doubles (fixed seed), written to reference/pyformat_reference.tsv.
Run once with CPython: python -I make_pyformat_reference.py"""
import os, random, struct
HERE = os.path.dirname(os.path.abspath(__file__))
hard = ['0.529113', '0.0920667', '6.5753507660701875e-11', '3.2546139159239827e-25', '1.9865758833475414e-31', '996151588624343.25', '96462868410162.625', '2.91e-11', '1.85', '1234565', '2.5', '0.125', '1.005', '0.0005',
        '100', '-3', '0.30000000000000004', '0.3333333333333333', '5e-324', '2.2250738585072014e-308', '1.7976931348623157e+308', '1e15', '1e16', '1e-4', '1e-5', '29400650', '-0.0207317230304731', '0.9920039980009996', '58.25049701789305', '0.0', '-0.0']
random.seed(2026)
texts = list(hard)
for _ in range(150): texts.append(repr(struct.unpack('<d', struct.pack('<Q', random.getrandbits(63)))[0]))      # any positive double, subnormals included
for _ in range(150): texts.append('%.*g' % (random.randint(1, 12), random.uniform(0, 1) * 10 ** random.randint(-12, 12)))      # short decimals as they appear in result files
texts = [t for t in texts if float(t) == float(t) and abs(float(t)) != float('inf')]
with open(os.path.join(HERE, 'reference', 'pyformat_reference.tsv'), 'w', encoding='ascii', newline='\n') as f:
    f.write('text\trepr\tg16\tg6\tf2\tf6\n')
    for t in texts:
        x = float(t); f.write('\t'.join([t, repr(x), '%.16g' % x, '%.6g' % x, '%.2f' % x, '%.6f' % x]) + '\n')
print(len(texts), 'reference rows')
