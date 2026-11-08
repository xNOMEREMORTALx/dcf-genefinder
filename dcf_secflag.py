"""dcf_secflag.py (DCF, L3161) - label-free selenocysteine (TGA read-through) candidate finder, run after dcf_genefinder_mobile.py.
For each predicted gene: (a) look upstream through an in-frame TGA to the most upstream start codon (upstream part >= 60 nt);
(b) look downstream through a TGA stop (continuation >= 120 nt). Developed in L3157 (recall 45/69 on 8 genomes); usable only on
AT-rich genomes (GC < 0.45), where false candidates are few. No dependencies. Usage: python3 dcf_secflag.py genome.gb pred.gff3
Author: Gary NOMEREMORTAL Mason. License: MIT."""
import sys, re
comp=str.maketrans('ACGTN','TGCAN')
def rc(s): return s.translate(comp)[::-1]
def read_seq(path):
    t=open(path).read()
    if '\nORIGIN' in t: return re.sub(r'[^ACGT]','N',''.join(re.findall(r'[A-Za-z]+',t.split('\nORIGIN')[1].split('\n//')[0])).upper())
    return re.sub(r'[^ACGT]','N',''.join(l.strip() for l in t.splitlines() if not l.startswith('>')).upper())
def read_gff(path):
    return [(int(f[3]),int(f[4]),f[6]) for f in (l.rstrip('\n').split('\t') for l in open(path)) if len(f)>6 and f[2]=='CDS']
STOP={'TAA','TAG','TGA'}; START={'ATG','GTG','TTG'}
def cands(seq,G,minup=60,mincont=120):
    L=len(seq); R=rc(seq); out={}
    for strand in '+-':
        S=seq if strand=='+' else R
        for a,b in sorted(((a-1,b) if strand=='+' else (L-b,L-a+1)) for a,b,s in G if s==strand):
            i=a-3
            while i>=0 and S[i:i+3] not in STOP: i-=3
            if i>=0 and S[i:i+3]=='TGA':
                j=i-3
                while j>=0 and S[j:j+3] not in STOP: j-=3
                us=[k for k in range(j+3,i,3) if S[k:k+3] in START]
                if us and i-us[0]>=minup: out[(strand,i)]='upstream'
            if S[b-3:b]=='TGA':
                k=b
                while k<L-3 and S[k:k+3] not in STOP: k+=3
                if k-b>=mincont: out.setdefault((strand,b-3),'downstream')
    return sorted((s,(i+1 if s=='+' else L-(i+3)+1),t) for (s,i),t in out.items())
if __name__=='__main__':
    seq=read_seq(sys.argv[1]); C=cands(seq,read_gff(sys.argv[2])); gc=(seq.count('G')+seq.count('C'))/len(seq)
    print("selenocysteine candidates: %d (GC %.2f%s)"%(len(C),gc,"" if gc<0.45 else "; GC >= 0.45: many false candidates expected, list not reliable"))
    for s,p,t in C: print("  TGA at %d (%s strand), %s read-through"%(p,s,t))
