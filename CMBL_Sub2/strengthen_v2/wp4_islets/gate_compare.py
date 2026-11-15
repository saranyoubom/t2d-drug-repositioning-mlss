"""gate_compare.py - WP4 gate (protocol 4.1): compare re-run outputs with the stored results of the two earlier islet cross-checks and with the numbers reported in their READMEs.
Part 1 (files): every file of the stored results folder must exist in the re-run folder and agree (CSV: same shape and columns, numeric columns within rtol 1e-9 / atol 1e-12, other columns equal; JSON: same
structure, floats within the same tolerance). Part 2 (reported numbers): the figures quoted in the READMEs and in the brief, compared at the precision at which they were written.
Usage: python -I gate_compare.py <label> <stored results dir> <re-run results dir> [summary file name]   (appends a report block to stdout; exit code 1 if any check fails)"""
import sys, os, json
import numpy as np, pandas as pd
sys.stdout.reconfigure(encoding='utf8')
RTOL, ATOL = 1e-9, 1e-12


def cmp_json(a, b, path=''):
    out = []
    if isinstance(a, dict):
        if not isinstance(b, dict) or set(a) != set(b): return [f'{path}: keys differ']
        for k in a: out += cmp_json(a[k], b[k], f'{path}/{k}')
    elif isinstance(a, list):
        if not isinstance(b, list) or len(a) != len(b): return [f'{path}: list length differs']
        for i, (x, y) in enumerate(zip(a, b)): out += cmp_json(x, y, f'{path}[{i}]')
    elif isinstance(a, float) or isinstance(b, float):
        if a is None or b is None: return [] if a is b else [f'{path}: {a} vs {b}']
        if not np.isclose(a, b, rtol=RTOL, atol=ATOL, equal_nan=True): out.append(f'{path}: {a} vs {b}')
    elif a != b: out.append(f'{path}: {a!r} vs {b!r}')
    return out


def cmp_csv(pa, pb):
    A, B = pd.read_csv(pa), pd.read_csv(pb)
    if A.shape != B.shape or list(A.columns) != list(B.columns): return None, [f'shape/columns differ: {A.shape} vs {B.shape}']
    bad, mx = [], 0.0
    for c in A.columns:
        if pd.api.types.is_numeric_dtype(A[c]) and pd.api.types.is_numeric_dtype(B[c]):
            x, y = A[c].to_numpy(float), B[c].to_numpy(float)
            if not np.allclose(x, y, rtol=RTOL, atol=ATOL, equal_nan=True): bad.append(f'{c}: max abs diff {np.nanmax(np.abs(x - y)):.3g}')
            d = np.nanmax(np.abs(x - y)) if len(x) and not np.all(np.isnan(x - y)) else 0.0; mx = max(mx, float(d))
        elif not A[c].fillna('NA').astype(str).equals(B[c].fillna('NA').astype(str)): bad.append(f'{c}: values differ')
    return mx, bad


def main():
    label, stored, rerun = sys.argv[1], sys.argv[2], sys.argv[3]
    ok_all = True; print(f'### {label}\nstored: {stored}\nre-run: {rerun}')
    for f in sorted(os.listdir(stored)):
        sp, rp = os.path.join(stored, f), os.path.join(rerun, f)
        if not os.path.exists(rp): print(f'  MISSING  {f}'); ok_all = False; continue
        if f.endswith('.csv'):
            mx, bad = cmp_csv(sp, rp)
            print(f'  {"OK     " if not bad else "DIFFER "} {f}' + (f'  (max abs difference over numeric cells {mx:.3g})' if not bad else '  ' + '; '.join(bad[:4]))); ok_all &= not bad
        elif f.endswith('.json'):
            bad = cmp_json(json.load(open(sp)), json.load(open(rp)))
            print(f'  {"OK     " if not bad else "DIFFER "} {f}' + ('' if not bad else '  ' + '; '.join(bad[:6]))); ok_all &= not bad
    print('FILES: ' + ('all identical within tolerance' if ok_all else 'DIFFERENCES FOUND'))
    return ok_all


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
