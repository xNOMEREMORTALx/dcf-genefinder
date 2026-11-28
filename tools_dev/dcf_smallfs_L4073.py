#!/usr/bin/env python3
# DCF L4073 - frameshift detector on a small genome (dev candidate v2.33dev_smallfs: 500 kb limit lowered to 100 kb).
# Reuses the L4041 mutation pipeline. Gary NOMEREMORTAL Mason. MIT.
import os, sys, re, subprocess
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dcf_longread_L4041 as LR

def main(gbpath, gf, outdir):
    os.makedirs(outdir, exist_ok=True)
    txt = open(gbpath).read(); head, feat, seq = LR.split_gb(txt)
    print("DCF L4073 tool output -", txt.split()[1])
    for e, seeds in ((0.0, [4073]), (1e-3, [4073, 4074, 4075])):
        for seed in seeds:
            rng = np.random.default_rng(seed); mseq, newidx, indels = LR.mutate(seq, e, rng)
            name = 'fs_e%g_s%d' % (e, seed); gb = os.path.join(outdir, name + '.gb'); gff = os.path.join(outdir, name + '.gff3')
            LR.write_gb(gb, name, head, LR.lift_features(feat, newidx), mseq)
            log = subprocess.run([sys.executable, '-I', gf, gb, '--table', 'auto', '--out', gff], capture_output=True, text=True).stdout
            fs = re.search(r'(\d+) possible frameshifts .*?self-estimated indel sensitivity (\d+)%, 95% interval (\d+)-(\d+)%', log)
            _, _, _, _, ann = LR.grade(gb, gff)
            inside = [p + 1 for p, _ in indels if any(a - 1 <= p <= b - 1 for a, b, s in ann)]
            allpos = np.array([p + 1 for p, _ in indels])
            pos = []
            fsf = gff.replace('.gff3', '') + '_frameshifts.txt'
            if os.path.exists(fsf):
                pos = [int(l.split()[-1]) for l in open(fsf) if l.strip() and not l.startswith('#')]
            pos = np.array(pos)
            hits = sum(1 for p in inside if len(pos) and np.min(np.abs(pos - p)) <= 150)
            true_rep = sum(1 for q in pos if len(allpos) and np.min(np.abs(allpos - q)) <= 150)
            print("FS e %g seed %d ran %d reported %d inside_cds %d hits150 %d reported_true150 %d est %s lo %s hi %s"
                  % (e, seed, int(fs is not None), len(pos), len(inside), hits, true_rep,
                     fs.group(2) if fs else 'NA', fs.group(3) if fs else 'NA', fs.group(4) if fs else 'NA'))
    print("END")

if __name__ == '__main__':
    main(*sys.argv[1:4])
