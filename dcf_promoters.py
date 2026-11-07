"""dcf_promoters.py (DCF, L3049) - label-free sigma-70/sigma-A promoter report, run after dcf_genefinder_mobile.py.
Learns the genome's own -35 / -10 boxes (seeded at TTGACA / TATAAT, spacer 15-19) from the 250 nt upstream of predicted operon-first
genes, and reports them only if the self-check passes: the 17-nt spacer share of the top-quarter sites must be > 1.5x that of the same
training on composition-shuffled windows. Usage: python3 dcf_promoters.py genome.gb_or_fasta predictions.gff3
Author: Gary NOMEREMORTAL Mason. License: MIT."""
import sys, re, random
import numpy as np
comp = str.maketrans('ACGT', 'TGCA')
def rc(s): return s.translate(comp)[::-1]
def read_seq(path):
    t = open(path).read()
    if '\nORIGIN' in t: return ''.join(re.findall(r'[acgtn]+', t.split('\nORIGIN')[1])).upper()
    return ''.join(l.strip() for l in t.splitlines() if not l.startswith('>')).upper()
def read_gff(path):
    G = []
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        if len(f) > 6 and f[2] == 'CDS': G.append((int(f[3]), int(f[4]), f[6] == '-'))
    return G
def windows(seq, genes, W=250):
    L = len(seq); out = []
    for strand in (False, True):
        g = sorted((a, b) for a, b, c in genes if c == strand)
        if not strand:
            prev = -10**9
            for a, b in g:
                if a - 1 - prev > W and a > W + 1: out.append(seq[a - W - 1:a - 1])
                prev = max(prev, b)
        else:
            prev = 10**12
            for a, b in sorted(g, key=lambda x: -x[1]):
                if prev - b > W and b + W <= L: out.append(rc(seq[b:b + W]))
                prev = min(prev, a)
    return [w for w in out if set(w) <= set('ACGT')]
IDX = {c: i for i, c in enumerate('ACGT')}
def onehot(w):
    a = np.zeros((len(w), 4)); a[np.arange(len(w)), [IDX[c] for c in w]] = 1; return a
def init(cons, bg):
    M = np.full((6, 4), 0.1)
    for k, c in enumerate(cons): M[k, IDX[c]] += 2.0
    M /= M.sum(1, keepdims=True); return np.log(M / bg)
def hexscores(X, P):
    n = len(X) - 5; return sum(X[k:k + n] @ P[k] for k in range(6))
def best(X, P35, P10):
    a = hexscores(X, P35); b = hexscores(X, P10); n = len(b); bv = -1e9; arg = None
    for s in range(15, 20):
        k = n - (6 + s)
        if k <= 0: continue
        v = a[:k] + b[6 + s:6 + s + k]; i = int(np.argmax(v))
        if v[i] > bv: bv = v[i]; arg = (i, i + 6 + s, s)
    return bv, arg
def train(Wn):
    Xs = [onehot(w) for w in Wn]; c = np.concatenate(Xs).sum(0); bg = c / c.sum()
    P35 = init('TTGACA', bg); P10 = init('TATAAT', bg)
    for _ in range(5):
        res = [best(X, P35, P10) for X in Xs]; top = np.argsort([-r[0] for r in res])[:max(1, len(Xs) // 4)]
        C35 = np.full((6, 4), 0.5); C10 = np.full((6, 4), 0.5)
        for j in top:
            i35, i10, s = res[j][1]; C35 += Xs[j][i35:i35 + 6]; C10 += Xs[j][i10:i10 + 6]
        P35 = np.log(C35 / C35.sum(1, keepdims=True) / bg); P10 = np.log(C10 / C10.sum(1, keepdims=True) / bg)
    res = [best(X, P35, P10) for X in Xs]; top = np.argsort([-r[0] for r in res])[:max(1, len(Xs) // 4)]
    sp = np.bincount([res[j][1][2] for j in top], minlength=20)[15:20]
    cons = lambda P: ''.join('ACGT'[i] for i in P.argmax(1))
    return sp / max(sp.sum(), 1), cons(P35), cons(P10), [(res[j][1], res[j][0]) for j in top]
def main(seqpath, gffpath, seed=7):
    seq = read_seq(seqpath); Wn = windows(seq, read_gff(gffpath))
    if len(Wn) < 100: print("promoter report: too few operon-first genes (%d)" % len(Wn)); return None
    fr, c35, c10, _ = train(Wn); random.seed(seed)
    fs, _, _, _ = train([''.join(random.sample(w, len(w))) for w in Wn])
    E17 = fr[2] / max(fs[2], 1e-9)
    print("promoter report: %d upstream windows; learned -35 %s, -10 %s; spacer 15-19 shares %s; 17-nt enrichment vs shuffled %.2f"
          % (len(Wn), c35, c10, ' '.join('%.2f' % x for x in fr), E17))
    if E17 > 1.5: print("  CLEAR sigma-70-type promoter signal (TTGACA-17-TATAAT family): promoter positions are meaningful in this genome")
    else: print("  no clear consensus-promoter signal: promoters here are too degenerate for a label-free scan (as in E. coli)")
    return E17
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
