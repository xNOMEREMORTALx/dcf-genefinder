"""dcf_metals.py (DCF, L3125) - label-free metal-use profile of a predicted proteome, run after dcf_genefinder_mobile.py.
Translates the predicted genes (standard code) and reports:
  * Fe-S / redox signature: C-x-x-C motifs per 1000 amino acids and R = observed / expected from the proteome's own cysteine level
    (L3117: methanogens 5.4-10.6, aerobic bacteria 2.6-5.7; densities 1.07 vs 0.32 per 1000 aa).
  * Classic C2H2 zinc fingers: C-x(2,4)-C-x(12)-H-x(3,5)-H per 1000 proteins (L3121: 0-1.7 in 18 prokaryotes; validated on 8 human
    proteins, 8/8 within +-1 of literature finger counts).
Usage: python3 dcf_metals.py genome.gb_or_fasta predictions.gff3       Author: Gary NOMEREMORTAL Mason. License: MIT."""
import sys, re
comp = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(comp)[::-1]
def read_seq(path):
    t = open(path).read()
    if '\nORIGIN' in t: return re.sub(r'[^ACGT]', 'N', ''.join(re.findall(r'[A-Za-z]+', t.split('\nORIGIN')[1].split('\n//')[0])).upper())
    return re.sub(r'[^ACGT]', 'N', ''.join(l.strip() for l in t.splitlines() if not l.startswith('>')).upper())
def read_gff(path):
    G = []
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        if len(f) > 6 and f[2] == 'CDS': G.append((int(f[3]), int(f[4]), f[6]))
    return G
B = 'TCAG'; AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}
ZF = re.compile(r'C.{2,4}C.{12}H.{3,5}H')
def main(seqpath, gffpath):
    seq = read_seq(seqpath); P = []
    for a, b, s in read_gff(gffpath):
        g = seq[a - 1:b] if s == '+' else rc(seq[a - 1:b])
        P.append(''.join(CODE.get(g[i:i + 3], 'X') for i in range(0, len(g) - 3, 3)))
    allp = ''.join(P); n = len(allp); fC = allp.count('C') / max(n, 1)
    cxxc = sum(len(re.findall(r'(?=C..C)', p)) for p in P); exp = sum(max(len(p) - 3, 0) for p in P) * fC * fC
    zf = sum(len(ZF.findall(p)) for p in P); zfp = sum(1 for p in P if ZF.search(p))
    R = cxxc / exp if exp else 0
    print("metal-use profile: %d proteins, %d aa, cysteine %.4f" % (len(P), n, fC))
    print("  Fe-S/redox  C-x-x-C: %.2f per 1000 aa, %.1fx the cysteine-level expectation" % (1000 * cxxc / max(n, 1), R))
    if fC < 0.004:
        print("  note: cysteine below 0.4%, so R is unstable (tiny expected count); use the density, not R (L3189 R-floor)")
    print("  Zn fingers  C2H2: %d motifs in %d proteins (%.2f proteins per 1000)" % (zf, zfp, 1000 * zfp / max(len(P), 1)))
    print("  reading: %s" % ("Fe-S rich (anaerobe-like, e.g. methanogens)" if 1000 * cxxc / max(n, 1) > 0.7 else "Fe-S moderate (aerobe-like)"))
    return R
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
