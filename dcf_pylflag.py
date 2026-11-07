"""dcf_pylflag.py (DCF, L3073) - label-free stop-codon read-through flag, run after dcf_genefinder_mobile.py.
Finds consecutive same-strand predicted genes joined in frame across a single stop codon (TAA/TAG/TGA) and reports the TAG excess over the
genome's own TAA/TGA joins (a pyrrolysine / amber-read-through signature). Usage: python3 dcf_pylflag.py genome.gb predictions.gff3
Author: Gary NOMEREMORTAL Mason. License: MIT."""
import sys, re
comp=str.maketrans('ACGTN','TGCAN')
def rc(s): return s.translate(comp)[::-1]
def read_seq(path):
    t=open(path).read()
    if '\nORIGIN' in t: return re.sub(r'[^ACGT]','N',''.join(re.findall(r'[A-Za-z]+',t.split('\nORIGIN')[1].split('\n//')[0])).upper())
    return re.sub(r'[^ACGT]','N',''.join(l.strip() for l in t.splitlines() if not l.startswith('>')).upper())
def read_gff(path):
    G=[]
    for l in open(path):
        f=l.rstrip('\n').split('\t')
        if len(f)>6 and f[2]=='CDS': G.append((int(f[3]),int(f[4]),f[6]))
    return G
STOPS=('TAA','TAG','TGA')
def joins(seq,genes,maxgap=600,minpart=150):
    L=len(seq); R=rc(seq); out=[]; use={s:0 for s in STOPS}
    for strand in '+-':
        S=seq if strand=='+' else R
        # convert to strand coordinates (0-based start, end exclusive) in transcription order
        g=sorted(((a-1,b) if strand=='+' else (L-b,L-a+1)) for a,b,s in genes if s==strand)
        for a,b in g:
            st=S[b-3:b]
            if st in use: use[st]+=1
        for (a1,b1),(a2,b2) in zip(g,g[1:]):
            st=S[b1-3:b1]
            if st not in STOPS or b1-a1<minpart or b2-a2<minpart: continue
            if a2<b1 or a2-b1>maxgap or (a2-b1)%3: continue
            if any(S[i:i+3] in STOPS for i in range(b1,a2,3)): continue
            out.append((strand,a1,b1,a2,b2,st))
    return out,use
def main(seqpath,gffpath):
    seq=read_seq(seqpath); J,use=joins(seq,read_gff(gffpath))
    n={s:sum(1 for j in J if j[5]==s) for s in STOPS}
    exp=(n['TAA']+n['TGA'])*use['TAG']/max(use['TAA']+use['TGA'],1)
    ex=n['TAG']/exp if exp>0 else float('inf') if n['TAG'] else 0.0
    print("read-through joins: TAA %d, TAG %d, TGA %d; stop use TAA %d TAG %d TGA %d; TAG excess %.2f"%(n['TAA'],n['TAG'],n['TGA'],use['TAA'],use['TAG'],use['TGA'],ex))
    if ex>=2: print("  POSSIBLE PYRROLYSINE / amber read-through: TAG joins are %.1fx the genome's TAA/TGA rate"%ex)
    return J,ex,len(seq)
if __name__=='__main__':
    main(sys.argv[1],sys.argv[2])
