"""dcf_neurometal v2 (L3425 amendment A1): label-free sorting of voltage-gated channels by the METAL ION each passes.
Author: Gary NOMEREMORTAL Mason (DCF). Not for medical use.
Input: RefSeq mRNA FASTA (.fna.gz); one transcript per gene (lowest NM_ number), longest ATG ORF.
Step 1 voltage sensor S4 (frozen L3293, passed blind) must be present.
Step 2 shape and pore:
  >= 1500 aa (4 domains): count the calcium ring  [TS]-x-[ED]-x-W  (one E/D per pore loop)
       >= 3  -> Ca (calcium channel, EEEE/EEDD ring)
       <= 1  AND sodium lysine loop [TS]-[FY]-[KR]-[GA]-W (the K of DEKA) -> Na
       otherwise -> undecided   [A1: v1 called any large S4 protein with ring <=1 Na; Brwd3 false call in mouse 1-2]
       == 2  -> undecided
  < 1500 aa: C-I-G-Y-G -> HCN (K/Na pacemaker); else [TSC]-[VITL]-G-[YF]-G (k_filter2) -> K; else other
Gene key = LAST parenthesis in the header [A1: v1 took the first, merging all mouse Kcnh genes under '(eag-related)'].
Mode 'tool': prints accession | aa | S4 | ring | class, NO gene names.  Mode 'grade': adds gene names.
Usage: python3 -I dcf_neurometal_L3425.py tool|grade file1.rna.fna.gz [...]"""
import gzip, re, sys
CODE = {}; b = 'TCAG'; aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, x in enumerate(b):
    for j, y in enumerate(b):
        for k, z in enumerate(b): CODE[x + y + z] = aa[16 * i + 4 * j + k]
def orf(s):
    best = ''
    for f in range(3):
        p = ''.join(CODE.get(s[i:i + 3], 'X') for i in range(f, len(s) - 2, 3))
        for seg in p.split('*'):
            k = seg.find('M')
            if k >= 0 and len(seg) - k > len(best): best = seg[k:]
    return best
H = '[AVILMFWYGSTC]'
S4 = r'(?=[RK]' + H + '{2}[RK]' + H + '{2}[RK]' + H + '{2}[RK]' + H + '{2}[RK])'
RING = r'(?=[TS].[ED].W)'
def classify(p):
    s4 = len(re.findall(S4, p))
    if not s4: return s4, 0, None
    r = len(re.findall(RING, p))
    if len(p) >= 1500: return s4, r, 'Ca' if r >= 3 else 'Na' if r <= 1 and re.search(r'[TS][FY][KR][GA]W', p) else 'undecided'
    if 'CIGYG' in p: return s4, r, 'HCN'
    if re.search(r'[TSC][VITL]G[YF]G', p): return s4, r, 'K'
    return s4, r, 'other'
if __name__ == '__main__':
    mode = sys.argv[1]; genes = {}
    for path in sys.argv[2:]:
        for blk in gzip.open(path, 'rt').read().split('>')[1:]:
            h, _, s = blk.partition('\n'); acc = h.split()[0]
            if not acc.startswith('NM_'): continue
            m = re.findall(r'\(([A-Za-z0-9\-\.]+)\)', h); g = m[-1] if m else acc
            num = int(re.sub(r'\D', '', acc.split('.')[0]))
            if g not in genes or num < genes[g][0]: genes[g] = (num, acc, s.replace('\n', '').upper())
    print('genes:', len(genes))
    rows = []
    for g, (_, acc, s) in genes.items():
        p = orf(s); s4, r, c = classify(p)
        if c: rows.append((c, acc, len(p), s4, r, g))
    for c, acc, n, s4, r, g in sorted(rows):
        print('%-9s %-14s %5d aa  S4 x%d  ring %d' % (c, acc, n, s4, r) + ('  ' + g if mode == 'grade' else ''))
