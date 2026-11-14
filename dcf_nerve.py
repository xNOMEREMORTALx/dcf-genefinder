"""dcf_nerve.py (DCF phone companion, L3433): nerve channels sorted by the metal ion they pass, from an mRNA FASTA.
Author: Gary NOMEREMORTAL Mason (Digital Cosmos Framework). License: MIT. Not for medical use.
Rules are exactly those that passed blind on mouse.3 (L3429, 10 of 10):
  dcf_neurometal v2 (L3425 A1): voltage sensor S4; >= 1500 aa: calcium ring [TS]x[ED]xW >= 3 -> Ca,
      <= 1 with the sodium lysine loop [TS][FY][KR][GA]W -> Na, else undecided; < 1500 aa: CIGYG -> HCN,
      [TSC][VITL]G[YF]G -> K, else other.
  dcf_k2p (L3425 A2): < 700 aa, no S4, >= 2 filters T[VIT]G[YFL]G -> K2P leak channel (anaesthetic targets TREK/TASK).
Input: RefSeq mRNA FASTA (.fna or .fna.gz); one transcript per gene (lowest NM_), longest ATG reading frame."""
import gzip, re
CODE = {}; _b = 'TCAG'; _aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, x in enumerate(_b):
    for j, y in enumerate(_b):
        for k, z in enumerate(_b): CODE[x + y + z] = _aa[16 * i + 4 * j + k]
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
    """returns (S4 count, ring count, class or None) - v2 rule"""
    s4 = len(re.findall(S4, p))
    if not s4: return s4, 0, None
    r = len(re.findall(RING, p))
    if len(p) >= 1500: return s4, r, 'Ca' if r >= 3 else 'Na' if r <= 1 and re.search(r'[TS][FY][KR][GA]W', p) else 'undecided'
    if 'CIGYG' in p: return s4, r, 'HCN'
    if re.search(r'[TSC][VITL]G[YF]G', p): return s4, r, 'K'
    return s4, r, 'other'
def is_k2p(p): return len(p) < 700 and not re.search(S4, p) and len(re.findall(r'T[VIT]G[YFL]G', p)) >= 2
def genes_from(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt', errors='replace') as h: txt = h.read()
    genes = {}
    for blk in txt.split('>')[1:]:
        h, _, s = blk.partition('\n'); acc = h.split()[0] if h.split() else ''
        if not acc.startswith('NM_'): continue
        m = re.findall(r'\(([A-Za-z0-9\-\.]+)\)', h); g = m[-1] if m else acc
        num = int(re.sub(r'\D', '', acc.split('.')[0]) or 0)
        if g not in genes or num < genes[g][0]: genes[g] = (num, acc, re.sub(r'\s', '', s).upper())
    return genes
NAMES = {'Ca': 'CALCIUM channels (transmitter release, heart, rhythm)', 'Na': 'SODIUM channels (the nerve spike)',
         'K': 'POTASSIUM channels (reset after each spike)', 'HCN': 'HCN pacemaker channels (K/Na)',
         'K2P': 'K2P potassium LEAK channels (resting potential; anaesthetic targets)', 'undecided': 'undecided 4-domain',
         'other': 'voltage sensor, no filter (not called)'}
def report(path):
    genes = genes_from(path); out = {}
    for g, (_, acc, s) in genes.items():
        p = orf(s); _, _, c = classify(p)
        if c is None and is_k2p(p): c = 'K2P'
        if c: out.setdefault(c, []).append("%s (%d aa)" % (g, len(p)))
    lines = ["DCF nerve channels by metal ion: %s  (%d genes)" % (path.replace('\\', '/').rsplit('/', 1)[-1], len(genes))]
    for c in ['Na', 'K', 'Ca', 'HCN', 'K2P', 'undecided', 'other']:
        if c in out: lines.append("  %-9s %3d  %s\n            %s" % (c, len(out[c]), NAMES[c], ", ".join(sorted(out[c]))))
    lines.append("  Blind record: mouse.3 (L3429) 30 calls, 0 wrong; 23/28 family genes; K2P 7/7. Not for medical use.")
    return "\n".join(lines)
if __name__ == '__main__':
    import sys; print(report(sys.argv[1]))
