#!/usr/bin/env python3
# DCF L4041 - gene finder on error-rich (long-read-style) assemblies. Gary NOMEREMORTAL Mason. MIT.
# Usage: python3 -I dcf_longread_L4041.py GENOME.gb GENEFINDER.py OUTDIR
import os, re, sys, subprocess
import numpy as np
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
import dcf_grade as G

RATES = [0.0, 1e-4, 1e-3, 5e-3]
SEEDS = [4041, 4042, 4043]

def split_gb(txt):
    head, rest = txt.split('\nFEATURES', 1)
    feat, origin = rest.split('\nORIGIN', 1)
    seq = ''.join(re.findall(r'[A-Za-z]+', origin.split('\n//')[0])).upper()
    return head, 'FEATURES' + feat, seq

def mutate(seq, e, rng):
    L = len(seq)
    if e == 0:
        return seq, np.arange(L + 1), []
    n = rng.poisson(e * L)
    n_indel = int(round(0.8 * n)); n_sub = n - n_indel
    runs = [(m.start(), m.end()) for m in re.finditer(r'(A+|C+|G+|T+)', seq) if m.end() - m.start() >= 3]
    ops = {}
    for i in rng.choice(len(runs), size=n_indel, replace=True):
        a, b = runs[i]
        ops[a] = 'ins' if rng.random() < 0.5 else 'del'          # act on the first base of the run
    for p in rng.choice(L, size=n_sub, replace=False):
        if p not in ops: ops[int(p)] = 'sub'
    out, newidx, indel_new = [], np.zeros(L + 1, dtype=np.int64), []
    for i, ch in enumerate(seq):
        newidx[i] = len(out)
        op = ops.get(i)
        if op == 'del':
            indel_new.append((len(out), i)); continue
        if op == 'ins':
            indel_new.append((len(out), i)); out.append(ch); out.append(ch); continue
        if op == 'sub':
            out.append('ACGT'[('ACGT'.index(ch) + 1 + int(rng.integers(3))) % 4] if ch in 'ACGT' else ch); continue
        out.append(ch)
    newidx[L] = len(out)
    return ''.join(out), newidx, indel_new

def lift_features(feat, newidx):
    lines, out, inloc = feat.split('\n'), [], False
    for ln in lines:
        if len(ln) > 21 and ln[:5] == '     ' and ln[5] != ' ':
            inloc = True
        elif ln.strip().startswith('/'):
            inloc = False
        if inloc and len(ln) > 21:
            loc = re.sub(r'\d+', lambda m: str(int(newidx[min(int(m.group()) - 1, len(newidx) - 2)]) + 1), ln[21:])
            ln = ln[:21] + loc
        out.append(ln)
    return '\n'.join(out)

def write_gb(path, name, head, feat, seq):
    head = re.sub(r'^LOCUS\s+\S+\s+\d+ bp', 'LOCUS       %s %d bp' % (name, len(seq)), head)
    with open(path, 'w') as fh:
        fh.write(head + '\n' + feat + '\nORIGIN\n')
        s = seq.lower()
        for i in range(0, len(s), 60):
            fh.write('%9d %s\n' % (i + 1, ' '.join(s[i + j:i + j + 10] for j in range(0, 60, 10))))
        fh.write('//\n')

def grade(gb, gff):
    seq, cds = G.read_gb(gb); ann = [(a, b, -1 if c else 1) for a, b, c in cds]
    P = []
    for l in open(gff):
        f = l.rstrip('\n').split('\t')
        if len(f) > 6 and f[2] == 'CDS': P.append((int(f[3]), int(f[4]), 1 if f[6] == '+' else -1))
    fd, cr = G.grade(P, ann)
    sk = {((b if s == 1 else a), s): (a, b) for a, b, s in P}; ex = n = 0
    for a, b, c in cds:
        s = -1 if c else 1; k = sk.get(((b if s == 1 else a), s))
        if k: n += 1; ex += (k[0] == a if s == 1 else k[1] == b)
    return 100 * fd, 100 * cr, 100 * ex / max(n, 1), len(P), ann

def main(gbpath, gf, outdir):
    os.makedirs(outdir, exist_ok=True)
    txt = open(gbpath).read()
    head, feat, seq = split_gb(txt)
    print("DCF L4041 tool output - %s, %d bp" % (txt.split()[1], len(seq)))
    for e in RATES:
        for seed in (SEEDS[:1] if e == 0 else SEEDS):
            rng = np.random.default_rng(seed)
            mseq, newidx, indels = mutate(seq, e, rng)
            name = 'm_e%g_s%d' % (e, seed)
            gb = os.path.join(outdir, name + '.gb'); gff = os.path.join(outdir, name + '.gff3')
            write_gb(gb, name, head, lift_features(feat, newidx), mseq)
            r = subprocess.run([sys.executable, '-I', gf, gb, '--table', 'auto', '--out', gff], capture_output=True, text=True)
            log = r.stdout + r.stderr
            open(os.path.join(outdir, name + '.log'), 'w').write(log)
            tab = re.search(r'genetic code table (\S+) -> predicted (\d+) genes', log)
            fs = re.search(r'(\d+) possible frameshifts .*?self-estimated indel sensitivity (\d+)%, 95% interval (\d+)-(\d+)%', log)
            f, c, x, npred, ann = grade(gb, gff)
            # indels inside annotated (graded) CDS, in new coordinates
            inside = [p for p, _ in indels if any(a - 1 <= p <= b - 1 for a, b, s in ann)]
            hits = 0
            fsfile = gff.replace('.gff3', '') + '_frameshifts.txt'
            pos = []
            if os.path.exists(fsfile):
                for l in open(fsfile):
                    if l.startswith('#') or not l.strip(): continue
                    parts = l.split()
                    try: pos.append(int(float(parts[-1])))
                    except ValueError: pass
            pos = np.array(sorted(pos))
            for p in inside:
                if len(pos) and np.min(np.abs(pos - (p + 1))) <= 30: hits += 1
            print("RUN e %g seed %d table %s npred %s found %.2f correct %.2f exact %.2f errors %d indels %d indels_in_cds %d "
                  "fs_reported %s fs_hits %d fs_est %s fs_lo %s fs_hi %s"
                  % (e, seed, tab.group(1) if tab else 'NA', tab.group(2) if tab else npred, f, c, x,
                     len(indels) + int(round(0)), len(indels), len(inside), fs.group(1) if fs else 'NA', hits,
                     fs.group(2) if fs else 'NA', fs.group(3) if fs else 'NA', fs.group(4) if fs else 'NA'))
    print("END")

if __name__ == '__main__':
    main(*sys.argv[1:4])
