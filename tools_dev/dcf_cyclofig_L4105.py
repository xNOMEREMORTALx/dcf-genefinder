#!/usr/bin/env python3
# DCF L4105 - which cyclotomic values Phi_m(10^k)/d are figurate numbers? Gary NOMEREMORTAL Mason. MIT.
from math import isqrt
import sympy
x = sympy.Symbol('x')
def sq(D):
    if D < 0: return None
    s = isqrt(D); return s if s * s == D else None
FAM = {
 'triangular':        lambda d: (lambda s: s is not None and (s - 1) % 2 == 0)(sq(8 * d + 1)),
 'square':            lambda d: sq(d) is not None,
 'pentagonal':        lambda d: (lambda s: s is not None and (1 + s) % 6 == 0)(sq(24 * d + 1)),
 'centred_triangular':lambda d: (lambda s: s is not None and (3 + s) % 6 == 0)(sq(24 * d - 15)),
 'centred_square':    lambda d: (lambda s: s is not None and (2 + s) % 4 == 0)(sq(8 * d - 4)),
 'centred_hexagonal': lambda d: (lambda s: s is not None and (3 + s) % 6 == 0)(sq(12 * d - 3)),
 'star':              lambda d: (lambda s: s is not None and (6 + s) % 12 == 0)(sq(24 * d + 12)),
}
# self-check of each membership test against generated members
gen = {'triangular': lambda n: n*(n+1)//2, 'square': lambda n: n*n, 'pentagonal': lambda n: n*(3*n-1)//2,
       'centred_triangular': lambda n: 3*n*(n-1)//2+1, 'centred_square': lambda n: 2*n*(n-1)+1,
       'centred_hexagonal': lambda n: 3*n*(n-1)+1, 'star': lambda n: 6*n*(n-1)+1}
print("DCF L4105 tool output")
for f in FAM:
    mem = {gen[f](n) for n in range(1, 3000)}; lim = max(mem)
    ok = all(FAM[f](v) == (v in mem) for v in range(1, min(lim, 200000)))
    print(f"GATE {f} membership_ok {int(ok)}")
hits = []
for m in range(1, 13):
    P = sympy.cyclotomic_poly(m, x)
    for k in range(1, 21):
        val = int(P.subs(x, 10 ** k))
        for d in (1, 2, 3, 4, 5, 6, 7, 9, 11, 13):
            if val % d: continue
            V = val // d
            for f, t in FAM.items():
                if t(V): hits.append((m, k, d, f, V))
for h in hits:
    m, k, d, f, V = h
    print(f"HIT m {m} k {k} d {d} family {f} value {V if V < 10**12 else str(V)[:8] + '...(' + str(len(str(V))) + ' digits)'} big {int(V > 1000)}")
print("END")
