"""dcf_nerve.py v2 (DCF phone/genefinder companion, L3549): nerve chemistry by metal, the nerve 'battery', and the
excitability (arousal/rhythm) parts list, read label-free from sequence.
Author: Gary NOMEREMORTAL Mason (Digital Cosmos Framework). License: MIT. Not for medical use.
Input: RefSeq mRNA FASTA (NM_/XM_, longest ATG ORF per gene), protein FASTA (e.g. Ensembl pep.all), or a UniProt TSV
(Entry + Sequence columns). .gz accepted.
Rules = dcf_neurochem6 v8 (blind record: mouse, zebrafish, chicken, fly, nematode, all Swiss-Prot animals, 24,796 TrEMBL Na/Ca
channels, early animals, comb jellies) + T-type calcium split (development L3549: 14/15 T-type, 0/30 false):
  channels: voltage sensor S4 (or weak sensor + K filter); 4-domain: sodium lysine T-F-K-G-W -> Na (DEKA), T-F-N-G-W -> Na sensor
            (Nax, DENA), calcium ring [TS]x[ED]xW >= 3 -> Ca (>= 2 ring sites with D -> T-type, else L/N/P/Q/R-type);
            1-domain: C-I-G-Y-G -> HCN, [TSC][VITL]G[YF]G -> K; K2P leak; glutamate receptor Y-T-A-N-L-A-A;
            iron hydroxylase [CVI]-H-E-[LI]-[LM]-G-H (dopamine/noradrenaline/serotonin synthesis);
  pumps (battery chargers): D-K-T-G-T-[LIVM]-[TI] + loop [NY]-V-P-E-G-L Na/K, [AS]-[IV]-P-E-G-L Ca, C-P-C + 2 x GM[TSHN]CxxC Cu."""
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
RING = r'(?=[TS].[ED].W)'; LYS = r'[TS][FY][KR][GA]W'; NAX = r'[TS]F[NQ]GW'
def weak(p):
    return bool(re.search(r'(?=[RK]..[RK]..[RK]..[RK])', p)) or any(sum(x in 'RK' for x in m.groups()) == 4
               for m in re.finditer(r'(?=([RKQ])..([RKQ])..([RKQ])..([RKQ])..([RKQ]))', p))
def channel(p):
    s4 = bool(re.search(S4, p)); r = len(re.findall(RING, p)); lys = bool(re.search(LYS, p)); nax = bool(re.search(NAX, p))
    kv7 = len(p) < 1500 and re.search(r'(?=R..Q..R..[RH])', p) and re.search(r'[TS][VI]GYGD', p)
    if not (s4 or kv7) and re.search(r'[TS][VI]GYGD|SVGFGN|CIGYG', p): kv7 = weak(p)
    if not (s4 or kv7): return None
    if len(p) >= 1500 and s4 and (r >= 3 or lys or nax):
        if lys: return 'Na'
        if nax: return 'Nax'
        if r >= 3:
            hits = [p[m.start():m.start() + 5] for m in re.finditer(RING, p)]
            return 'CaT' if sum(h[2] == 'D' for h in hits) >= 2 else 'Ca'
        return 'undecided'
    if 'CIGYG' in p: return 'HCN'
    if re.search(r'[TSC][VITL]G[YF]G', p): return 'K'
    return None
def pump_class(p):
    if not re.search(r'DKTGT[LIVM][TI]', p): return None
    if re.search(r'[NY]VPEGL', p): return 'pump_NaK'
    if re.search(r'[AS][IV]PEGL', p): return 'pump_Ca'
    if 'CPC' in p and len(re.findall(r'GM[TSHN]C..C', p)) >= 2: return 'pump_Cu'
    return 'pump_other'
def is_k2p(p): return len(p) < 1100 and not re.search(S4, p) and len(re.findall(r'[TS][TSV][IVT]G[YFL]G', p)) >= 2
def classify(p):
    """kept for older callers: (S4 count, ring count, channel class)"""
    return len(re.findall(S4, p)), len(re.findall(RING, p)), channel(p)
def tags(p):
    t = []; c = channel(p)
    if c: t.append(c)
    if is_k2p(p): t.append('K2P')
    if 'YTANLAA' in p: t.append('iGluR')
    if re.search(r'[CVI]HE[LI][LM]GH', p): t.append('FeHyd')
    pc = pump_class(p)
    if pc: t.append(pc)
    return t
def proteins(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt', errors='replace') as h: txt = h.read()
    out = {}
    if not txt.startswith('>'):
        L = txt.splitlines(); hd = L[0].split('\t'); ie = hd.index('Entry'); iq = hd.index('Sequence')
        for l in L[1:]:
            f = l.split('\t')
            if len(f) > iq: out[f[ie]] = f[iq]
        return out
    for blk in txt.split('>')[1:]:
        h, _, s = blk.partition('\n'); seq = re.sub(r'\s', '', s).upper().rstrip('*'); acc = h.split()[0] if h.split() else ''
        if acc.startswith(('NR_', 'XR_')): continue   # non-coding RNA records
        if acc.startswith(('NM_', 'XM_')):
            m = re.findall(r'\(([A-Za-z0-9\-\.:]+)\)', h); g = m[-1] if m else acc; p = orf(seq)
        else:
            m = re.search(r'gene_symbol:(\S+)|gene:(\S+)', h); g = (m.group(1) or m.group(2)) if m else acc
            p = seq if set(seq[:200]) - set('ACGTN') else orf(seq)
        if g not in out or len(p) > len(out[g]): out[g] = p
    return out
NAMES = {'Na': 'SODIUM channels, DEKA (the nerve spike)', 'Nax': 'SODIUM SENSOR Nax (DENA)', 'K': 'POTASSIUM channels (reset after each spike)',
         'Ca': 'CALCIUM channels, high-voltage (transmitter release)', 'CaT': 'CALCIUM channels, T-type (rhythmic bursts, sleep spindles)',
         'HCN': 'HCN pacemakers (rhythm drive)', 'K2P': 'K2P LEAK channels (resting level; anaesthetic targets)',
         'iGluR': 'glutamate receptors (excitation, learning)', 'FeHyd': 'IRON hydroxylases (dopamine, noradrenaline, serotonin)',
         'pump_NaK': 'Na/K PUMPS (charge the battery)', 'pump_Ca': 'CALCIUM pumps', 'pump_Cu': 'COPPER pumps', 'pump_other': 'P-type ATPases, no metal',
         'undecided': 'undecided 4-domain'}
ORDER = ['Na', 'Nax', 'K', 'Ca', 'CaT', 'HCN', 'K2P', 'iGluR', 'FeHyd', 'pump_NaK', 'pump_Ca', 'pump_Cu', 'pump_other', 'undecided']
def report(path, show=40):
    P = proteins(path); out = {}
    for g, p in P.items():
        for t in tags(p): out.setdefault(t, []).append(g)
    n = lambda k: len(out.get(k, []))
    L = ["DCF nerve chemistry by metal: %s  (%d genes/proteins)" % (path.replace('\\', '/').rsplit('/', 1)[-1], len(P))]
    for c in ORDER:
        if c in out:
            names = sorted(out[c]); L.append("  %-10s %3d  %s\n             %s%s" % (c, len(names), NAMES[c], ", ".join(names[:show]), ' ...' if len(names) > show else ''))
    ch = sum(n(k) for k in ('Na', 'Nax', 'K', 'Ca', 'CaT', 'HCN', 'K2P')); pu = sum(n(k) for k in ('pump_NaK', 'pump_Ca', 'pump_Cu'))
    L.append("  BATTERY: %d dischargers (channels) : %d chargers (metal pumps)%s" % (ch, pu, ' = %.1f per pump' % (ch / pu) if pu else ''))
    L.append("  EXCITABILITY KIT (arousal and rhythm parts; not a measure of consciousness): pacemaker HCN %d, T-type Ca %d, K2P leak %d,"
             " monoamine iron enzymes %d, glutamate receptors %d" % (n('HCN'), n('CaT'), n('K2P'), n('FeHyd'), n('iGluR')))
    L.append("  Sodium spike type: %s" % ('DEKA present (bilaterian-type fast spike)' if n('Na') else 'no DEKA sodium channel (early-animal or ambulacrarian type, or none)'))
    L.append("  Not for medical use.")
    return "\n".join(L)
if __name__ == '__main__':
    import sys; print(report(sys.argv[1]))
