#!/usr/bin/env python3
"""dcf_anchor.py - does a gene family sit at chromatin-loop anchors more (or less) than matched genes?
Gary NOMEREMORTAL Mason, Digital Cosmos Framework. MIT licence. Not for medical use. numpy only.

The DCF matched-gene, dependence-aware test (Manuscript 3; rounds L3725, L3825a, L3853-L3869, L3909):
  gene table  GENCODE GTF, protein-coding genes on chr1-22 and X; TSS by strand; length = end - start + 1
  at anchor   the gene's TSS lies within 5 kb of a loop anchor (anchor interval widened by +-5 kb)
  density     number of other protein-coding TSSs within +-1 Mb
  matching    each member gets its 50 nearest pool genes in (z log10 length, z log10(1 + density)), Euclidean, stable ties
  null        10,000 sets drawing one candidate per member; all maps scored together (mean share over maps), so genes and
              anchors shared between maps keep their dependence
  output      observed share, null mean and sd, excess, p_upper = (1 + #null >= obs)/(1 + N), p_lower likewise

Usage
  python3 -I dcf_anchor.py --gtf gencode.v49.annotation.gtf.gz --genes KCNMA1,KCNJ5,...   (or --regex '^KCN' or --genes-file f)
        --loops map1.bedpe.gz [map2 ...] [--format top|4dn|all] [--ref ref1 ref2 --mode both|none]
        [--exclude REGEX] [--seed 4037] [--sets 10000]
  python3 -I dcf_anchor.py --selftest            synthetic planted-signal and calibration check (no data needed)
Loop formats: top = ENCODE loops, top 28,270 by column 7 (else the first 28,270); 4dn = only loops whose column 7 lists
CTCF-ChIAPET; all = every loop. 'chr' is added to bare chromosome names. --ref files (4DN union) give reference CTCF anchors;
--mode both keeps loops with both anchors on a reference anchor, none keeps loops with neither.
"""
import sys, re, gzip, argparse, collections
import numpy as np

def opener(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def gene_table(gtf):
    G = []
    for l in opener(gtf):
        if l.startswith('#'): continue
        x = l.split('\t')
        if len(x) < 9 or x[2] != 'gene' or 'gene_type "protein_coding"' not in x[8] or not re.match(r'^chr([0-9]+|X)$', x[0]): continue
        n = re.search(r'gene_name "([^"]+)"', x[8]).group(1)
        a, b = int(x[3]), int(x[4])
        G.append((x[0], a if x[6] == '+' else b, b - a + 1, n))
    return G

def covariates(G):
    chrom = np.array([g[0] for g in G]); tss = np.array([g[1] for g in G]); length = np.array([g[2] for g in G], float)
    dens = np.zeros(len(G))
    for c in np.unique(chrom):
        idx = np.where(chrom == c)[0]; T = tss[idx]; o = np.sort(T)
        dens[idx] = np.searchsorted(o, T + 1_000_000, side='right') - np.searchsorted(o, T - 1_000_000, side='left') - 1
    ll = np.log10(length); ld = np.log10(1 + dens)
    return (ll - ll.mean()) / ll.std(), (ld - ld.mean()) / ld.std()

def read_loops(p, fmt):
    L = []
    for l in opener(p):
        if l.startswith('#'): continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 6 or not f[1].isdigit(): continue
        if fmt == '4dn' and (len(f) < 7 or 'CTCF-ChIAPET' not in f[6].split(',')): continue
        L.append(f)
    if fmt == 'top' and len(L) > 28270:
        def sc(f):
            try: return float(f[6])
            except (IndexError, ValueError): return None
        if all(sc(f) is not None for f in L): L = [f for _, f in sorted(enumerate(L), key=lambda t: (-sc(t[1]), t[0]))][:28270]
        else: L = L[:28270]
    fix = lambda c: c if c.startswith('chr') else 'chr' + c
    return [(fix(f[0]), int(f[1]), int(f[2]), fix(f[3]), int(f[4]), int(f[5])) for f in L]

def ref_overlap(refs):
    CI = collections.defaultdict(list)
    for p in refs:
        for c1, s1, e1, c2, s2, e2 in read_loops(p, '4dn'):
            CI[c1].append((s1, e1)); CI[c2].append((s2, e2))
    CS = {c: (np.array(sorted(v))[:, 0], np.array(sorted(v))[:, 1]) for c, v in CI.items()}
    MW = max(int((CS[c][1] - CS[c][0]).max()) for c in CS)
    def ov(c, s, e):
        if c not in CS: return False
        st, en = CS[c]; i0 = np.searchsorted(st, s - MW); i1 = np.searchsorted(st, e)
        return bool((en[i0:i1] > s).any())
    return ov

def anchors_of(loops, ov=None, mode=None):
    A = collections.defaultdict(list)
    for c1, s1, e1, c2, s2, e2 in loops:
        if ov is not None:
            k = ov(c1, s1, e1) + ov(c2, s2, e2)
            if (mode == 'both' and k < 2) or (mode == 'none' and k > 0): continue
        A[c1].append((s1, e1)); A[c2].append((s2, e2))
    return A

def hits(G, A, pad=5000):
    H = np.zeros(len(G), bool); byc = collections.defaultdict(list)
    for i, g in enumerate(G): byc[g[0]].append(i)
    for c, idx in byc.items():
        if c not in A: continue
        a = np.array(A[c]); o = np.argsort(a[:, 0], kind='stable'); s = a[o, 0] - pad; e = a[o, 1] + pad
        me = np.maximum.accumulate(e); T = np.array([G[i][1] for i in idx]); k = np.searchsorted(s, T, side='right') - 1
        H[idx] = (k >= 0) & (me[np.maximum(k, 0)] >= T)
    return H

def matched_test(Hm, members, pool, zl, zd, nsets=10000, seed=4037, k=50):
    C = np.array([pool[np.argsort(np.hypot(zl[pool] - zl[i], zd[pool] - zd[i]), kind='stable')[:k]] for i in members])
    rng = np.random.default_rng(seed); obs = Hm[members].mean()
    null = np.array([Hm[C[np.arange(len(C)), rng.integers(0, k, len(C))]].mean() for _ in range(nsets)])
    bal = (float(zl[members].mean() - zl[C].mean()), float(zd[members].mean() - zd[C].mean()))   # residual covariate gap
    # L4053: members that share candidates share the sampling error of the control pool's own mean. Effective number of
    # distinct controls D = 1 / sum_j w_j^2 (w_j = chance a random null draw picks gene j); variance inflation 1 + m/D.
    u, cnt = np.unique(C.ravel(), return_counts=True); w = cnt / cnt.sum(); D = 1.0 / np.sum(w ** 2)
    sd_adj = null.std() * np.sqrt(1 + len(members) / D)
    from math import erf, sqrt
    p_adj = 0.5 * (1 - erf(((obs - null.mean()) / sd_adj if sd_adj > 0 else 0.0) / sqrt(2)))
    return dict(n=len(members), obs=obs, mean=null.mean(), sd=null.std(), excess=obs - null.mean(), balance=bal, D_eff=D, p_upper_adj=p_adj,
                p_upper=(1 + (null >= obs).sum()) / (nsets + 1), p_lower=(1 + (null <= obs).sum()) / (nsets + 1))

def naive_test(Hm, members, pool, nsets=10000, seed=4037):
    rng = np.random.default_rng(seed); obs = Hm[members].mean()
    null = np.array([Hm[rng.choice(pool, len(members), replace=False)].mean() for _ in range(nsets)])
    return (1 + (null >= obs).sum()) / (nsets + 1)

# ---------------- synthetic self-test (L4037) ----------------
def synthetic(rng, n_genes=6000, n_chr=10, b_len=0.9):
    chrom = rng.integers(0, n_chr, n_genes)
    pos = rng.random(n_genes) * 1.5e8
    length = 10 ** rng.normal(4.3, 0.5, n_genes)
    G = [('chr%d' % (c + 1), int(p), int(L), 'G%d' % i) for i, (c, p, L) in enumerate(zip(chrom, pos, length))]
    zl, zd = covariates(G)
    base = 1 / (1 + np.exp(-(-1.2 + b_len * zl - 0.4 * zd)))     # anchor probability vs length (the confounder; L4069: b_len < 0 = real direction)
    return G, zl, zd, base

def selftest(reps_power=20, reps_null=100, nsets=2000):
    rng = np.random.default_rng(4037)
    power = cal = naive_fp = 0
    for r in range(reps_power + reps_null):
        G, zl, zd, base = synthetic(rng)
        long_idx = np.where(zl > 1.0)[0]
        members = rng.choice(long_idx, 60, replace=False)
        p = base.copy()
        planted = r < reps_power
        if planted: p[members] = np.minimum(1, p[members] + 0.15)
        H = (rng.random(len(G)) < p).astype(float)
        pool = np.setdiff1d(np.arange(len(G)), members)
        res = matched_test(H, members, pool, zl, zd, nsets=nsets, seed=int(rng.integers(1 << 30)))
        if planted:
            power += res['p_upper'] < 0.01
        else:
            cal += res['p_upper'] < 0.05
            naive_fp += naive_test(H, members, pool, nsets=nsets, seed=int(rng.integers(1 << 30))) < 0.05
    print("SELFTEST planted_detected %d/%d matched_false_pos %d/%d naive_false_pos %d/%d"
          % (power, reps_power, cal, reps_null, naive_fp, reps_null))
    return power, cal, naive_fp

def selftest2(seed=4049, reps_power=20, reps_null=100, nsets=2000, ks=(50, 10)):
    """L4049: same synthetic design as L4037, fresh seed; compares neighbour counts k on identical replicates and reports
    the residual covariate gap (members minus controls, in SD units of z log length)."""
    rng = np.random.default_rng(seed)
    pw = {k: 0 for k in ks}; fp = {k: 0 for k in ks}; gap = {k: [] for k in ks}
    for r in range(reps_power + reps_null):
        G, zl, zd, base = synthetic(rng)
        members = rng.choice(np.where(zl > 1.0)[0], 60, replace=False)
        p = base.copy(); planted = r < reps_power
        if planted: p[members] = np.minimum(1, p[members] + 0.15)
        H = (rng.random(len(G)) < p).astype(float); pool = np.setdiff1d(np.arange(len(G)), members)
        s0 = int(rng.integers(1 << 30))
        for k in ks:
            res = matched_test(H, members, pool, zl, zd, nsets=nsets, seed=s0, k=k); gap[k].append(res['balance'][0])
            if planted: pw[k] += res['p_upper'] < 0.01
            else: fp[k] += res['p_upper'] < 0.05
    for k in ks:
        print("SELFTEST2 k %d planted_detected %d/%d matched_false_pos %d/%d mean_zlength_gap %.4f"
              % (k, pw[k], reps_power, fp[k], reps_null, float(np.mean(gap[k]))))

def selftest3(seed=4053, reps_power=20, reps_null=100, nsets=2000, k=50):
    """L4053: same design, fresh seed; empirical p_upper vs the variance-inflated p_upper_adj on identical replicates."""
    rng = np.random.default_rng(seed); pw = [0, 0]; fp = [0, 0]; D = []
    for r in range(reps_power + reps_null):
        G, zl, zd, base = synthetic(rng)
        members = rng.choice(np.where(zl > 1.0)[0], 60, replace=False)
        p = base.copy(); planted = r < reps_power
        if planted: p[members] = np.minimum(1, p[members] + 0.15)
        H = (rng.random(len(G)) < p).astype(float); pool = np.setdiff1d(np.arange(len(G)), members)
        res = matched_test(H, members, pool, zl, zd, nsets=nsets, seed=int(rng.integers(1 << 30)), k=k); D.append(res['D_eff'])
        for j, pv in enumerate((res['p_upper'], res['p_upper_adj'])):
            if planted: pw[j] += pv < 0.01
            else: fp[j] += pv < 0.05
    print("SELFTEST3 k %d D_eff_mean %.1f empirical: planted %d/%d false_pos %d/%d | adjusted: planted %d/%d false_pos %d/%d"
          % (k, float(np.mean(D)), pw[0], reps_power, fp[0], reps_null, pw[1], reps_power, fp[1], reps_null))

def selftest4(seed=4069, reps_null=200, nsets=2000, k=50):
    """L4069: realistic confounder direction (longer genes sit at anchors LESS, as L3853 found); null replicates only."""
    rng = np.random.default_rng(seed); fu = fl = 0; gaps = []
    for r in range(reps_null):
        G, zl, zd, base = synthetic(rng, b_len=-0.9)
        members = rng.choice(np.where(zl > 1.0)[0], 60, replace=False)
        H = (rng.random(len(G)) < base).astype(float); pool = np.setdiff1d(np.arange(len(G)), members)
        res = matched_test(H, members, pool, zl, zd, nsets=nsets, seed=int(rng.integers(1 << 30)), k=k)
        fu += res['p_upper'] < 0.05; fl += res['p_lower'] < 0.05; gaps.append(res['balance'][0])
    print("SELFTEST4 k %d null_reps %d false_pos_upper %d false_pos_lower %d mean_zlength_gap %.4f" % (k, reps_null, fu, fl, float(np.mean(gaps))))

def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--selftest', action='store_true'); ap.add_argument('--selftest2', action='store_true'); ap.add_argument('--selftest3', action='store_true'); ap.add_argument('--selftest4', action='store_true'); ap.add_argument('--k', type=int, default=50)
    ap.add_argument('--gtf'); ap.add_argument('--genes'); ap.add_argument('--genes-file'); ap.add_argument('--regex')
    ap.add_argument('--loops', nargs='+'); ap.add_argument('--format', default='top', choices=['top', '4dn', 'all'])
    ap.add_argument('--ref', nargs='*'); ap.add_argument('--mode', choices=['both', 'none'])
    ap.add_argument('--exclude', default=None); ap.add_argument('--seed', type=int, default=4037); ap.add_argument('--sets', type=int, default=10000)
    a = ap.parse_args()
    if a.selftest:
        selftest(); return
    if a.selftest2:
        selftest2(); return
    if a.selftest3:
        selftest3(); return
    if a.selftest4:
        selftest4(); return
    G = gene_table(a.gtf); names = [g[3] for g in G]; zl, zd = covariates(G)
    if a.regex: want = {n for n in names if re.match(a.regex, n)}
    else:
        want = set((a.genes or '').split(',')) if a.genes else {l.strip() for l in open(a.genes_file) if l.strip()}
    members = np.array([i for i, n in enumerate(names) if n in want])
    excl = re.compile(a.exclude) if a.exclude else None
    pool = np.array([i for i, n in enumerate(names) if n not in want and not (excl and excl.match(n))])
    ov = ref_overlap(a.ref) if a.ref else None
    Hs = []
    for p in a.loops:
        A = anchors_of(read_loops(p, a.format), ov, a.mode); Hs.append(hits(G, A).astype(float))
        print('map %s anchor share (all genes) %.4f' % (p, Hs[-1].mean()))
    r = matched_test(np.mean(Hs, axis=0), members, pool, zl, zd, a.sets, a.seed, a.k)
    print('genes %d (missing from GTF: %s)' % (r['n'], ','.join(sorted(want - set(names))) or 'none'))
    print('T_obs %.4f null mean %.4f sd %.4f excess %+.4f p_upper %.4f p_lower %.4f'
          % (r['obs'], r['mean'], r['sd'], r['excess'], r['p_upper'], r['p_lower']))
    print('balance (members minus matched controls, SD units): z log length %+.3f, z log density %+.3f' % r['balance'])
    print('effective distinct controls %.1f for %d genes; variance-inflated p_upper_adj %.4f (L4053)' % (r['D_eff'], r['n'], r['p_upper_adj']))

if __name__ == '__main__':
    main()
