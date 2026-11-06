"""dcf_grade.py - grade gene predictions (GFF3) against a GenBank annotation.
found   = annotated CDS covered >= 50% by a same-strand prediction
correct = prediction covered >= 50% by a same-strand annotated CDS
exact start = predicted start equals the annotated start, among annotated CDS whose stop is predicted on the right strand
Usage: python3 dcf_grade.py genome.gb predictions.gff3
Author: Gary NOMEREMORTAL Mason (Digital Cosmos Framework)."""
import re, sys
import numpy as np

def read_gb(path):
    txt = open(path).read()
    seq = ''.join(re.findall(r'[acgtn]+', txt.split('\nORIGIN')[1])).upper() if '\nORIGIN' in txt else ''
    cds = []
    for m in re.finditer(r'\n     CDS             (complement\()?<?(\d+)\.\.>?(\d+)\)?\n', txt):
        a, b = int(m.group(2)), int(m.group(3)); cds.append((a, b, bool(m.group(1))))
    return seq, cds

def grade(pred, ann):
    import bisect
    def cov(x, others):
        a, b, s = x; best = 0
        for (c, d, t) in others:
            if t == s and c <= b and d >= a: best = max(best, min(b, d) - max(a, c) + 1)
        return best >= 0.5 * (b - a + 1)
    # speed: bucket by 10 kb
    from collections import defaultdict
    def bucket(L):
        B = defaultdict(list)
        for x in L:
            for k in range(x[0] // 10000, x[1] // 10000 + 1): B[k].append(x)
        return B
    BA, BP = bucket(ann), bucket(pred)
    near = lambda x, B: {y for k in range(x[0] // 10000, x[1] // 10000 + 1) for y in B[k]}
    f = sum(cov(x, near(x, BP)) for x in ann) / len(ann); c = sum(cov(x, near(x, BA)) for x in pred) / len(pred)
    return f, c

def main(gb, gff):
    seq, cds = read_gb(gb); ann = [(a, b, -1 if c else 1) for a, b, c in cds]
    P = []
    for l in open(gff):
        f = l.rstrip('\n').split('\t')
        if len(f) > 6 and f[2] == 'CDS': P.append((int(f[3]), int(f[4]), 1 if f[6] == '+' else -1))
    fd, cr = grade(P, ann); sk = {((b if s == 1 else a), s): (a, b) for a, b, s in P}; ex = n = 0
    for a, b, c in cds:
        s = -1 if c else 1; k = sk.get(((b if s == 1 else a), s))
        if k: n += 1; ex += (k[0] == a if s == 1 else k[1] == b)
    print("found %.1f%%  correct %.1f%%  exact start %.1f%% (of %d same-stop genes)" % (100*fd, 100*cr, 100*ex/max(n, 1), n))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
