"""dcf_motif_scan.py (DCF, L3297) - scan proteins for DCF protein-family motifs, label-free.
Author: Gary NOMEREMORTAL Mason (Digital Cosmos Framework). License: MIT.

Input, any of:
  genome.gb + predictions.gff3      (proteins from dcf_genefinder_mobile.py calls)
  mRNA FASTA (.fna/.fa, may be .gz) (RefSeq NM_ only, one transcript per gene = lowest NM_ number, longest ATG reading frame)
  protein FASTA (.faa)
Presets (status = how far each motif has been tested; 'dev' = development only, not yet blind):
"""
import gzip, re, sys

PRESETS = {
  "hedgehog":  (r"H[YF]EGRA[VI]D.{40,80}GCF", "Hedgehog signal protein: zinc site + G-C-F self-cutting site (cholesterol)", "dev L3297: IHH only of 26,088 human transcripts; 0 in 70 microbes"),
  "fd_ladder": (r"C..C..C...C", "4Fe-4S ferredoxin 3-3-4 ladder", "blind: 14/16 tiers (L3269), photosynthesis P1/P2 (L3293)"),
  "cuA":       (r"C.E.C...H..M", "copper CuA site (oxygen respiration)", "100% precision, 12 genomes (L3257)"),
  "zn_ribbon": (r"C..C.{9,20}C..C", "zinc ribbon knuckle pair (structural Zn; unmasked here, dcf_redox masks ladder cysteines)", "blind L3273/L3285 (masked version)"),
  "rubredoxin":(r"C..CG.{24,32}C..CG", "rubredoxin iron knuckle pair", "blind RB1/RB2 (L3285)"),
  "c2h2_zf":   (r"C.{2,4}C.{12}H.{3,5}H", "C2H2 zinc finger", "L3261 (zinc-finger rhythm)"),
  "cchc":      (r"C..C....H....C", "CCHC zinc knuckle (RNA binding)", "dev L3281: 11/12 known human genes"),
  "volt_s4":   (r"(?=[RK][AVILMFWYGSTC]{2}[RK][AVILMFWYGSTC]{2}[RK][AVILMFWYGSTC]{2}[RK][AVILMFWYGSTC]{2}[RK])", "voltage sensor S4 (nerve, muscle channels)", "test L3293: 11/13 recall, 11/11 precision"),
  "k_filter":  (r"[TSC][VIT]G[YF]G", "potassium-channel selectivity filter", "test L3293: recall 11/13, precision 0.42 (failed bar 0.45)"),
  "k_filter2": (r"[TSC][VITL]G[YF]G", "potassium filter incl. Kv3/Kv4 T-T-L-G-Y-G (correction found L3301)", "UNTESTED correction of L3293 N1"),
  "cys_loop":  (r"C.{4,6}[FY]P[FYWM]D.{3}C", "Cys-loop neurotransmitter receptor (corrected for GABA-A F/Y-P-M-D)", "UNTESTED correction of L3293 N3 (failed)"),
}
__doc__ += "\n".join("  %-11s %s  [%s]" % (k, v[1], v[2]) for k, v in PRESETS.items())

CODE = {}
_b = 'TCAG'; _aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, x in enumerate(_b):
    for j, y in enumerate(_b):
        for k, z in enumerate(_b): CODE[x + y + z] = _aa[16 * i + 4 * j + k]
comp = str.maketrans('ACGTN', 'TGCAN')

def longest_orf(s):
    best = ''
    for f in range(3):
        p = ''.join(CODE.get(s[i:i + 3], 'X') for i in range(f, len(s) - 2, 3))
        for seg in p.split('*'):
            k = seg.find('M')
            if k >= 0 and len(seg) - k > len(best): best = seg[k:]
    return best

def _open(path): return gzip.open(path, 'rt', errors='replace') if path.endswith('.gz') else open(path, errors='replace')

def proteins_from(path, gff=None):
    """yields (name, protein)"""
    if gff:
        sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
        import dcf_redox as R
        seq = R.read_seq(path)
        for n, p in enumerate(R.proteins(seq, gff), 1): yield ("gene_%d" % n, p)
        return
    txt = _open(path).read()
    if txt.lstrip().startswith('LOCUS'): sys.exit("GenBank genome: give the .gff3 from the gene finder as the second argument")
    blocks = [b.partition('\n') for b in txt.split('>')[1:]]
    seqs = [(h, re.sub(r'\s', '', s).upper()) for h, _, s in blocks]
    is_dna = sum(set(s[:300]) <= set('ACGTN') for _, s in seqs[:50]) > len(seqs[:50]) / 2
    if not is_dna:
        for h, s in seqs: yield (h.split()[0], s)
        return
    best = {}
    for h, s in seqs:
        acc = h.split()[0]
        if not acc.startswith('NM_') and any(x[0].startswith('NM_') for x in seqs[:5]): continue
        m = re.search(r'\(([A-Za-z0-9\-\.]+)\)', h); g = m.group(1) if m else acc
        num = int(re.sub(r'\D', '', acc.split('.')[0]) or 0)
        if g not in best or num < best[g][0]: best[g] = (num, s)
    for g, (_, s) in best.items(): yield (g, longest_orf(s))

def scan(path, gff=None, names=None, custom=None):
    pats = {k: PRESETS[k][0] for k in (names or PRESETS)}
    if custom: pats["custom"] = custom
    hits = {k: [] for k in pats}; n = aa = 0
    for name, p in proteins_from(path, gff):
        n += 1; aa += len(p)
        for k, pat in pats.items():
            c = len(re.findall(pat, p))
            if c:
                m = re.search(pat, p)
                hits[k].append((name, c, m.start() + 1))
    return n, aa, hits

def report(path, gff=None, names=None, custom=None, show=12):
    n, aa, hits = scan(path, gff, names, custom)
    out = ["DCF motif scan: %s  (%d proteins, %d aa)" % (path.rsplit('/', 1)[-1], n, aa)]
    for k, h in hits.items():
        tot = sum(c for _, c, _ in h)
        desc = PRESETS[k][1] if k in PRESETS else custom
        out.append("  %-11s %5d proteins, %.4f per 1000 aa  - %s" % (k, len(h), 1000.0 * tot / max(aa, 1), desc))
        if 0 < len(h) <= show or k in ("hedgehog", "custom"):
            out.append("              " + ", ".join("%s (x%d, first at aa %d)" % x for x in h[:show]))
    out.append("  status of each motif: python3 dcf_motif_scan.py (no arguments). Not for medical use.")
    return "\n".join(out)

if __name__ == "__main__":
    a = sys.argv[1:]
    if not a: sys.exit(__doc__)
    custom = None
    if "--motif" in a: i = a.index("--motif"); custom = a[i + 1]; a = a[:i] + a[i + 2:]
    sel = None
    if "--only" in a: i = a.index("--only"); sel = a[i + 1].split(","); a = a[:i] + a[i + 2:]
    print(report(a[0], a[1] if len(a) > 1 else None, sel, custom))
