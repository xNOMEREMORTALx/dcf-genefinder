"""dcf_align (L3569): global protein alignment (Needleman-Wunsch, BLOSUM62, linear gap 8), numpy row-vectorized, no installs.
identity(a,b) = identical aligned positions / aligned (non-gap) positions. Author: Gary NOMEREMORTAL Mason (DCF)."""
import numpy as np
_B = """A R N D C Q E G H I L K M F P S T W Y V
A 4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
R -1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
N -2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
D -2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
C 0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
E -1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
G 0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
H -2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3
L -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1
K -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2
M -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1
F -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2
S 1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2
T 0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1
V 0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4"""
_L = _B.split('\n'); AA = _L[0].split(); IDX = {a: i for i, a in enumerate(AA)}
M = np.array([[int(x) for x in l.split()[1:]] for l in _L[1:]], dtype=np.int32)
M = np.pad(M, ((0, 1), (0, 1)), constant_values=-1)  # index 20 = unknown
def enc(s): return np.array([IDX.get(c, 20) for c in s], dtype=np.int64)
def align(a, b, gap=8):
    A, B = enc(a), enc(b); n, m = len(A), len(B)
    H = np.zeros((n + 1, m + 1), dtype=np.int32); H[0] = -gap * np.arange(m + 1); H[:, 0] = -gap * np.arange(n + 1)
    T = np.zeros((n + 1, m + 1), dtype=np.int8)   # 0 diag, 1 up, 2 left
    ar = np.arange(1, m + 1, dtype=np.int32) * gap
    for i in range(1, n + 1):
        diag = H[i - 1, :-1] + M[A[i - 1], B]; up = H[i - 1, 1:] - gap
        best = np.maximum(diag, up); t = np.where(diag >= up, 0, 1).astype(np.int8)
        cand = np.maximum.accumulate(np.concatenate(([H[i, 0] + 0], best + ar)))[1:] - ar
        left = cand > best
        H[i, 1:] = np.maximum(best, cand); t[left] = 2; T[i, 1:] = t
    i, j = n, m; same = cols = 0
    while i > 0 and j > 0:
        k = T[i, j]
        if k == 0:
            cols += 1; same += A[i - 1] == B[j - 1]; i -= 1; j -= 1
        elif k == 1: i -= 1
        else: j -= 1
    return same / max(cols, 1), int(H[n, m])
def identity(a, b): return align(a, b)[0]
