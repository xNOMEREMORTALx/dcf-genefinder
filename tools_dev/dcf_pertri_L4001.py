#!/usr/bin/env python3
# DCF L4001 - triangular numbers among numbers with exact reciprocal period n (Gary NOMEREMORTAL Mason). MIT. Not for medical use.
# Pre-registration: preregistrations/PREREGISTERED_L4001_period_triangular.txt
# Usage: python3 -I dcf_pertri_L4001.py > L4001_tool_output.txt
import math, sys
from math import gcd, isqrt
import sympy
import mpmath

def lcm(a, b):
    return a // gcd(a, b) * b

def factor_rep(n):
    """Factor 10^n - 1 via cyclotomic pieces Phi_k(10), k | n. Returns {p: e}, all p proved prime."""
    x = sympy.Symbol('x')
    fac = {}
    for k in sympy.divisors(n):
        val = int(sympy.cyclotomic_poly(k, x).subs(x, 10))
        for p, e in sympy.factorint(val).items():
            fac[p] = fac.get(p, 0) + e
    prod = 1
    for p, e in fac.items():
        assert sympy.isprime(p), p
        prod *= p ** e
    assert prod == 10 ** n - 1
    return fac

def exact_period_set(n):
    fac = factor_rep(n)
    items = []  # per prime: list of (p^k, ord) for k = 0..e
    for p, e in sorted(fac.items()):
        opts = [(1, 1)]
        for k in range(1, e + 1):
            q = p ** k
            opts.append((q, int(sympy.n_order(10, q))))
        items.append(opts)
    out = []
    def rec(i, val, o):
        if i == len(items):
            if o == n and val > 1:
                out.append(val)
            return
        for q, oq in items[i]:
            no = lcm(o, oq)
            if n % no:
                continue
            rec(i + 1, val * q, no)
    rec(0, 1, 1)
    return fac, out

PR = (3, 5, 7, 11, 13)
def weight(d):
    w = 1.0 / math.sqrt(2.0 * d)
    for p in PR:
        w *= 2.0 if d % p == 0 else (1 - 2 / p) / (1 - 1 / p)
    return w

def is_tri(d):
    s = isqrt(8 * d + 1)
    return s * s == 8 * d + 1

def period_digits(t, n):
    # integer long division: one period of 1/t, n digits, leading zeros kept
    r = 1
    ds = []
    for _ in range(n):
        r *= 10
        ds.append(r // t)
        r %= t
    assert r == 1
    s1 = ''.join(map(str, ds))
    mpmath.mp.dps = 12000
    s2 = mpmath.nstr(mpmath.mpf(1) / t, 10000, strip_zeros=False).split('.')[1][:n]
    assert s1 == s2, (t, s1, s2)
    return s1

def run(n, show=True):
    fac, D = exact_period_set(n)
    tri = sorted(d for d in D if is_tri(d))
    E = sum(weight(d) for d in D)
    return fac, D, tri, E

print("DCF L4001 tool output - triangular numbers with exact reciprocal period n")
# Gate G0
_, D24, tri24, E24 = run(24)
_, D16, tri16, E16 = run(16)
print(f"GATE n=24 |D|={len(D24)} tri={tri24} E={E24:.4f}")
print(f"GATE n=16 |D|={len(D16)} smallest_tri={tri16[0] if tri16 else None} tri={tri16} E={E16:.4f}")
for n in range(25, 41):
    fac, D, tri, E = run(n)
    fs = ' * '.join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(fac.items()))
    print(f"N {n} size {len(D)} O {len(tri)} E {E:.6f} tri {tri}")
    print(f"  factors 10^{n}-1 = {fs}")
    for t in tri:
        k = (isqrt(8 * t + 1) - 1) // 2
        s = period_digits(t, n)
        print(f"  T {t} = T_{k} digitsum {sum(map(int, s))} period {s}")
print("END")
