"""dcf_redox.py (DCF, L3277) - energy-metal profile of a predicted proteome, run after dcf_genefinder_mobile.py.
Author: Gary NOMEREMORTAL Mason (Digital Cosmos Framework). License: MIT.

Reports, per 1000 amino acids of the predicted proteins (standard code):
  * IRON   3-3-4 ferredoxin ladder C-x2-C-x2-C-x3-C (4Fe-4S)          - frozen L3237
  * IRON   plain C-x-x-C                                              - frozen L3125
  * COPPER CuA site C-x-E-x-C-x3-H-x2-M (proteins carrying it)          - frozen L3257
  * ZINC   ribbon C-x2-C-x(9..20)-C-x2-C, ladder cysteines masked        - frozen L3273
  * ZINC   knuckle (CCHC) C-x2-C-x4-H-x4-C                            - L3277, descriptive only (not yet tested)
and the four-tier energy call from the ladder, cuts 0.0045 | 0.023 | 0.060:
  tier 1 lactic fermenter | tier 2 oxygen respirer | tier 3 H2-making fermenter | tier 4 ferredoxin-based respirer/autotroph.
Blind record (L3269/L3273): 14 of 16 fresh genomes in the right tier, Spearman 0.96. Boundary cases: microaerophiles with
ferredoxin enzymes (Arcobacter, Sulfurihydrogenibium) read as tier 3. NOT for medical or clinical use.

Usage: python3 dcf_redox.py genome.gb_or_fasta predictions.gff3"""
import sys, re

CODE = {}
_b = 'TCAG'; _aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, x in enumerate(_b):
    for j, y in enumerate(_b):
        for k, z in enumerate(_b): CODE[x + y + z] = _aa[16 * i + 4 * j + k]
comp = str.maketrans('ACGTN', 'TGCAN')
CUTS = (0.0045, 0.023, 0.060)
TIERS = {1: "lactic fermenter (no ferredoxin energy step)", 2: "oxygen respirer",
         3: "H2-making fermenter (ferredoxin feeds a hydrogenase)", 4: "ferredoxin-based respirer / autotroph"}

def read_seq(path):
    t = open(path, errors='replace').read()
    if '\nORIGIN' in t:
        return re.sub(r'[^ACGT]', 'N', ''.join(re.findall(r'[A-Za-z]+', t.split('\nORIGIN')[1].split('\n//')[0])).upper())
    return re.sub(r'[^ACGT]', 'N', ''.join(l.strip() for l in t.splitlines() if not l.startswith('>')).upper())

def proteins(seq, gff):
    for l in open(gff):
        f = l.split('\t')
        if len(f) > 6 and f[2] == 'CDS':
            a, e, s = int(f[3]), int(f[4]), f[6]
            g = seq[a - 1:e] if s == '+' else seq[a - 1:e].translate(comp)[::-1]
            yield ''.join(CODE.get(g[i:i + 3], 'X') for i in range(0, len(g) - 3, 3))

def tier(ladder):
    return 1 if ladder < CUTS[0] else 2 if ladder < CUTS[1] else 3 if ladder < CUTS[2] else 4

def profile(genome, gff):
    seq = read_seq(genome); n = lad = cx = rib = knk = cua = 0
    for p in proteins(seq, gff):
        n += len(p)
        L = [m.start() for m in re.finditer(r'(?=C..C..C...C)', p)]; lad += len(L)
        cx += len(re.findall(r'(?=C..C)', p))
        q = list(p)
        for st in L:
            for k in range(st, st + 11):
                if q[k] == 'C': q[k] = 'z'
        rib += len(re.findall(r'(?=C..C.{9,20}C..C)', ''.join(q)))
        knk += len(re.findall(r'(?=C..C....H....C)', p))
        cua += bool(re.search(r'C.E.C...H..M', p))
    k = 1000.0 / max(n, 1)
    return dict(aa=n, ladder=lad * k, cxxc=cx * k, ribbon=rib * k, knuckle=knk * k, cua=cua, tier=tier(lad * k))

def bar(x, top, w=40):
    return '#' * max(0, min(w, int(round(w * x / top))))

def report(d, name=""):
    out = ["DCF energy-metal profile %s (%d aa of predicted protein)" % (name, d['aa']),
           "  IRON   ferredoxin ladder %.4f /1000 aa  %s" % (d['ladder'], bar(d['ladder'], 0.15)),
           "  IRON   C-x-x-C           %.4f /1000 aa" % d['cxxc'],
           "  COPPER CuA proteins      %d" % d['cua'],
           "  ZINC   ribbon            %.4f /1000 aa" % d['ribbon'],
           "  ZINC   knuckle (CCHC)    %.4f /1000 aa  (descriptive, untested)" % d['knuckle'],
           "  ENERGY TIER %d: %s" % (d['tier'], TIERS[d['tier']]),
           "  scale  |1 lactic <0.0045 |2 O2 <0.023 |3 H2 ferm <0.060 |4 Fd-based|"]
    if CUTS[1] - 0.002 <= d['ladder'] <= CUTS[1] + 0.002 or CUTS[2] - 0.003 <= d['ladder'] <= CUTS[2] + 0.003:
        out.append("  note: reading is within a boundary band; the neighbouring tier is also plausible")
    if d['tier'] == 3 and d['cua'] > 0:
        out.append("  note: tier 3 with CuA - possibly a microaerophile with ferredoxin enzymes (boundary case)")
    return "\n".join(out)

if __name__ == "__main__":
    if len(sys.argv) < 3: sys.exit(__doc__)
    print(report(profile(sys.argv[1], sys.argv[2]), sys.argv[1]))
