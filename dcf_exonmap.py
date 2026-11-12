"""dcf_exonmap.py (DCF, L3297) - place an mRNA on a stretch of chromosome: exons, introns, strand, and chromosome coordinates.
Author: Gary NOMEREMORTAL Mason (Digital Cosmos Framework). License: MIT. Exact-match seed-and-extend; no external tools.

Usage: python3 dcf_exonmap.py region.fasta mrna.fasta [region_start]
  region.fasta  a chromosome slice, e.g. NCBI efetch ...&id=NC_000007.14&seq_start=155794980&seq_stop=155817463&rettype=fasta
  mrna.fasta    the RefSeq mRNA (NM_...)
  region_start  chromosome coordinate of the slice's first base (read from the FASTA header ':155794980-155817463' if omitted)
Prints each exon (mRNA positions -> chromosome positions), intron lengths with their GT..AG splice signals, strand,
the gene span, and what fraction of the mRNA was placed (the poly-A tail is not in the genome)."""
import re, sys

comp = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(comp)[::-1]

def read_fasta(path):
    head, seq = '', []
    for l in open(path, errors='replace'):
        if l.startswith('>'):
            if head: break
            head = l[1:].strip()
        else: seq.append(l.strip())
    return head, re.sub(r'[^ACGT]', 'N', ''.join(seq).upper())

def exons_on(genome, mrna, seed=25, min_exon=12):
    """greedy left-to-right placement of mRNA blocks in increasing genome order. Returns [(m_start, m_end, g_start, g_end)] 0-based half-open."""
    out, i, gpos = [], 0, 0
    while i < len(mrna) - seed:
        k = genome.find(mrna[i:i + seed], gpos)
        if k < 0: i += 1; continue
        j = seed
        while i + j < len(mrna) and k + j < len(genome) and mrna[i + j] == genome[k + j]: j += 1
        b = 0                                   # extend left into a previous skipped stretch, if it matches
        while i - b - 1 >= (out[-1][1] if out else 0) and k - b - 1 >= gpos and mrna[i - b - 1] == genome[k - b - 1]: b += 1
        if j + b >= min_exon:
            out.append((i - b, i + j, k - b, k + j)); gpos = k + j
        i += j
    return out

def splice_adjust(G, mrna, ex, w=6):
    """an exon boundary is ambiguous when the bases next to it match on both sides; slide it (within the ambiguous
    stretch only) to the position that gives the canonical GT..AG intron, if one exists."""
    ex = [list(e) for e in ex]
    for n in range(len(ex) - 1):
        a, b = ex[n], ex[n + 1]
        if a[1] != b[0]: continue                  # mRNA not contiguous across this junction
        best = None
        for d in range(-w, w + 1):
            if d > 0 and (mrna[a[1]:a[1] + d] != G[a[3]:a[3] + d] or b[1] - b[0] <= d): continue
            if d < 0 and (mrna[a[1] + d:a[1]] != G[b[2] + d:b[2]] or a[1] - a[0] <= -d): continue
            intr = G[a[3] + d:b[2] + d]
            if intr[:2] == 'GT' and intr[-2:] == 'AG' and (best is None or abs(d) < abs(best)): best = d
        if best:
            a[1] += best; a[3] += best; b[0] += best; b[2] += best
    return [tuple(e) for e in ex]

def best_strand(genome, mrna):
    plus = splice_adjust(genome, mrna, exons_on(genome, mrna)); minus = splice_adjust(rc(genome), mrna, exons_on(rc(genome), mrna))
    cov = lambda ex: sum(e[1] - e[0] for e in ex)
    return ('+', plus) if cov(plus) >= cov(minus) else ('-', minus)

def report(region_path, mrna_path, region_start=None):
    gh, g = read_fasta(region_path); mh, m = read_fasta(mrna_path)
    if region_start is None:
        r = re.search(r':c?(\d+)-(\d+)', gh); region_start = int(r.group(1)) if r else 1
    strand, ex = best_strand(g, m)
    L = len(g); out = ["DCF exon map: %s  ->  %s" % (mh[:70], gh[:70]),
                       "  region %d bp starting at chromosome position %d; mRNA %d nt; strand %s" % (L, region_start, len(m), strand)]
    def chrom(gs, ge):          # 0-based half-open on the strand used -> 1-based inclusive chromosome coordinates
        if strand == '+': return region_start + gs, region_start + ge - 1
        return region_start + (L - ge), region_start + (L - gs) - 1
    G = g if strand == '+' else rc(g)
    for n, (ms, me, gs, ge) in enumerate(ex, 1):
        a, b = chrom(gs, ge)
        out.append("  exon %d: mRNA %d-%d (%d nt) -> chromosome %d-%d" % (n, ms + 1, me, me - ms, min(a, b), max(a, b)))
        if n < len(ex):
            ig = G[ge:ex[n][2]]
            out.append("     intron %d: %d bp, starts %s ... ends %s%s" % (n, len(ig), ig[:2], ig[-2:], "  (GT..AG splice signals)" if ig[:2] == 'GT' and ig[-2:] == 'AG' else ""))
    placed = sum(e[1] - e[0] for e in ex)
    if ex:
        c = [chrom(e[2], e[3]) for e in ex]; lo = min(min(x) for x in c); hi = max(max(x) for x in c)
        out.append("  gene span on the chromosome: %d-%d (%d bp), %d exons; %.1f%% of the mRNA placed" % (lo, hi, hi - lo + 1, len(ex), 100.0 * placed / len(m)))
    else: out.append("  mRNA not found in this region")
    return "\n".join(out), dict(strand=strand, exons=ex, placed=placed, mrna_len=len(m), span=(lo, hi) if ex else None)

if __name__ == "__main__":
    if len(sys.argv) < 3: sys.exit(__doc__)
    print(report(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else None)[0])
