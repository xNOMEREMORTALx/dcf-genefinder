"""dcf_k2p (L3425 amendment A2): two-pore-domain potassium (K2P 'leak') channels, label-free.
K2P channels set the resting potential and are targets of general anaesthetics (TREK/TASK). Rule:
 protein < 700 aa, NO voltage sensor S4 (frozen L3293), and >= 2 potassium filters T-[VIT]-G-[YFL]-G.
Gene key = last parenthesis; one transcript per gene (lowest NM_), longest ATG ORF. Not for medical use.
Usage: python3 -I dcf_k2p_L3425.py tool|grade file.rna.fna.gz [...]   Author: Gary NOMEREMORTAL Mason (DCF)."""
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
F = r'T[VIT]G[YFL]G'
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
    for g, (_, acc, s) in sorted(genes.items(), key=lambda x: x[1][1]):
        p = orf(s)
        if len(p) < 700 and not re.search(S4, p) and len(re.findall(F, p)) >= 2:
            print('K2P  %-14s %4d aa  %s' % (acc, len(p), ' '.join(re.findall(F, p))) + ('  ' + g if mode == 'grade' else ''))
