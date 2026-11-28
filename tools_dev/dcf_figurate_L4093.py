#!/usr/bin/env python3
# DCF L4093 - do Gary's grid numbers (star S_n = 6n(n-1)+1, centred hexagonal H_n = 3n(n-1)+1) favour any reciprocal period?
# Periods 25..40 (L4001 sets). Null: P(d in family) = c(d mod m) / (dF/dn at n(d)), c = exact residue counts, m = 360360.
# Gary NOMEREMORTAL Mason. MIT.
import os, sys, math
from math import isqrt
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib.util
sp = importlib.util.spec_from_file_location('PT', os.path.join(HERE, 'dcf_pertri_L4001.py')); PT = importlib.util.module_from_spec(sp)
sys.argv = ['x']; src = open(os.path.join(HERE, 'dcf_pertri_L4001.py')).read().split('print("DCF L4001 tool output')[0]
g = {}; exec(compile(src, 'pertri', 'exec'), g); exact_period_set = g['exact_period_set']
m = 360360
def counts(poly, nmod):
    n = np.arange(nmod, dtype=np.int64); r = poly(n) % m
    return np.bincount(r, minlength=m)
FAM = {
  'tri':  (lambda n: n * (n + 1) // 2, 2 * m, 0.5, lambda d: math.sqrt(2 * d),        lambda d: (lambda s: s * s == 8 * d + 1)(isqrt(8 * d + 1))),
  'star': (lambda n: 6 * n * (n - 1) + 1, m, 1.0, lambda d: 2 * math.sqrt(6 * d),      lambda d: (lambda D, s: s * s == D and (6 + s) % 12 == 0)(24 * d + 12, isqrt(24 * d + 12))),
  'hex':  (lambda n: 3 * n * (n - 1) + 1, m, 1.0, lambda d: 2 * math.sqrt(3 * d),      lambda d: (lambda D, s: s * s == D and (3 + s) % 6 == 0)(12 * d - 3, isqrt(max(12 * d - 3, 0)))),
}
print("DCF L4093 tool output")
C = {k: counts(v[0], v[1]) for k, v in FAM.items()}
for n in list(range(25, 41)) + [16, 24]:
    _, D = exact_period_set(n)
    out = [f"N {n} size {len(D)}"]
    for k, (poly, nmod, scale, deriv, member) in FAM.items():
        O = sorted(d for d in D if member(d)); E = sum(C[k][d % m] * scale / deriv(d) for d in D)
        out.append(f"{k} O {len(O)} E {E:.5f} members {O}")
    print(" | ".join(out))
# gate: brute force membership check on 1..2e6 against each formula
for k, (poly, nmod, scale, deriv, member) in FAM.items():
    true = set(); i = 1
    while poly(np.int64(i)) <= 2_000_000:
        v = int(poly(np.int64(i)));  true.add(v) if v > 0 else None; i += 1
    got = {d for d in range(1, 2_000_001) if member(d)}
    print(f"GATE {k} membership_exact {int(true == got)}")
print("END")
