"""dcf_localalign (L3585): Smith-Waterman local alignment score (BLOSUM62, linear gap 8), numpy row-vectorized.
score(a,b) = best local alignment score (length-robust: only the best-matching region counts)."""
import numpy as np, runpy, os
_G = runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dcf_align.py'))
M, enc = _G['M'], _G['enc']
def score(a, b, gap=8):
    A, B = enc(a), enc(b); m = len(B)
    prev = np.zeros(m + 1, dtype=np.int32); ar = np.arange(1, m + 1, dtype=np.int32) * gap; best = 0
    for i in range(len(A)):
        diag = prev[:-1] + M[A[i], B]; up = prev[1:] - gap
        h = np.maximum(np.maximum(diag, up), 0)
        cand = np.maximum.accumulate(np.concatenate(([0], h + ar)))[1:] - ar
        row = np.maximum(h, cand); row = np.maximum(row, 0)
        best = max(best, int(row.max())); prev = np.concatenate(([0], row)).astype(np.int32)
    return best
