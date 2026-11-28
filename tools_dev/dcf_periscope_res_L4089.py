#!/usr/bin/env python3
# DCF L4089 - corner camera resolution: smallest separation of two bright patches that can be told apart. Gary NOMEREMORTAL Mason. MIT.
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib.util
sp = importlib.util.spec_from_file_location('P', os.path.join(HERE, 'dcf_periscope_L4081.py')); P = importlib.util.module_from_spec(sp); sp.loader.exec_module(P)
N, M = P.N, P.M
def recon(A, y, sigma):                  # weak prior: mean 0, variance 1
    lam = sigma ** 2 / 1.0
    return np.linalg.solve(A.T @ A + lam * np.eye(N), A.T @ y)
def resolved(xh, i, j):
    lo, hi = min(i, j), max(i, j)
    pk = []
    for c in (lo, hi):
        w = range(max(0, c - 2), min(N, c + 3)); k = max(w, key=lambda q: xh[q]); pk.append((k, xh[k]))
    if pk[0][0] >= pk[1][0]: return False
    dip = xh[pk[0][0]:pk[1][0] + 1].min()
    return dip <= 0.75 * min(pk[0][1], pk[1][1]) and min(pk[0][1], pk[1][1]) > 0
def smin(occ, noise, rng):
    A = P.transport(occ)
    for s in range(2, 21):
        ok = 0
        for _ in range(10):
            i = int(rng.integers(4, N - 4 - s)); j = i + s
            x = np.zeros(N); x[i] = x[j] = 1.0
            y0 = A @ x; sig = noise * y0.mean(); y = y0 + rng.normal(0, sig, M)
            ok += resolved(recon(A, y, sig), i, j)
        if ok >= 8: return s
    return 99
rng = np.random.default_rng(4089); print("DCF L4089 tool output")
for noise in (0.01, 0.10):
    for occ in (True, False):
        print(f"SMIN noise {noise} occluder {int(occ)} s_min {smin(occ, noise, rng)}")
print("END")
