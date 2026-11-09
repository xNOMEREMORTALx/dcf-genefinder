#!/usr/bin/env python3
"""dcf_genefinder_mobile.py v2.32 — DCF label-free spectral gene finder, NO Biopython needed (works in Pydroid 3 on Android).
Gary NOMEREMORTAL Mason, Digital Cosmos Framework (L2723; algorithm = dcf_genefinder v2: period-3 self-training,
shadow removal, GC frame-plot seed for GC-rich genomes, anti-Shine-Dalgarno starts from the genome's own 16S rRNA).
Needs only numpy (Pydroid: Pip menu -> install numpy).

Run with arguments:   python3 dcf_genefinder_mobile.py GENOME.gb [--table auto|11|4|2|3|5|9] [--out pred.gff3]
Run with no arguments (Pydroid 'play' button): it asks for the file name and the genetic-code table.
"""
import sys, os, re
try:
    import numpy as np
except ImportError:
    sys.exit("numpy is required: in Pydroid open the menu -> Pip -> install 'numpy'")

STOPS={'11':{'TAA','TAG','TGA'},'4':{'TAA','TAG'},'3':{'TAA','TAG'},'5':{'TAA','TAG'},'9':{'TAA','TAG'},'2':{'TAA','TAG','AGA','AGG'}}
# L2810 (v2.17): TGA read-through threshold. Dev gap (25 genomes): standard code <= +0.042, TGA=Trp >= +0.084; midpoint of the 95% CI gap.
TGA_THR=0.065
SEED='gc3'      # L2819: GC-rich seed; 'gc2' (second codon position GC-poorest) is used only by the v2.19 rescue
LAST={}
COMP=str.maketrans('ACGTRYKMBDHVN','TGCAYRMKVHDBN')
def revcomp(s): return s.translate(COMP)[::-1]

# ---------- minimal GenBank / FASTA reader (replaces Biopython) ----------
def parse_location(loc):
    """returns (list of (start0,end), strand) for simple/complement/join locations"""
    strand=-1 if loc.startswith('complement(') else 1
    nums=re.findall(r'<?(\d+)\.\.>?(\d+)',loc)
    if not nums:
        m=re.findall(r'(\d+)',loc); nums=[(m[0],m[0])] if m else []
    parts=[(int(a)-1,int(b)) for a,b in nums]
    return parts,strand
def read_genbank(path):
    txt=open(path,errors='replace').read()
    if txt.lstrip().startswith('>'):                      # FASTA
        lines=txt.strip().splitlines(); rid=lines[0][1:].split()[0]
        return rid,''.join(l.strip() for l in lines[1:]).upper(),[]
    rid=(re.search(r'^VERSION\s+(\S+)',txt,re.M) or re.search(r'^LOCUS\s+(\S+)',txt,re.M)).group(1)
    seq=''.join(re.findall(r'[A-Za-z]+',txt.split('\nORIGIN',1)[1].split('\n//')[0])).upper() if '\nORIGIN' in txt else ''
    feats=[]
    if 'FEATURES' in txt:
        block=txt.split('\nFEATURES',1)[1].split('\nORIGIN',1)[0].splitlines()[1:]
        cur=None
        for line in block:
            if len(line)>21 and line[5]!=' ':                      # new feature
                if cur: feats.append(cur)
                cur={'type':line[5:21].strip(),'loc':line[21:].strip(),'q':{}}; lastq=None
            elif cur is not None:
                s=line[21:].strip()
                if s.startswith('/'):
                    k,_,v=s[1:].partition('='); cur['q'][k]=v.strip('"'); lastq=k
                elif lastq is None: cur['loc']+=s
                else: cur['q'][lastq]+=' '+s.strip('"')
        if cur: feats.append(cur)
    for f in feats: f['parts'],f['strand']=parse_location(f['loc'])
    return rid,seq,feats
def extract(seq,f):
    s=''.join(seq[a:b] for a,b in f['parts'])
    return s if f['strand']==1 else revcomp(s)

# ---------- the gene finder ----------
code={c:i for i,c in enumerate('ACGT')}
def arr(s): return np.array([code.get(c,0) for c in s],np.int8)
def selftrain(F,R,AF,AR,stops,minl,seedl,mask=None,iters=3,best_start=False):
    """period-3 self-training; returns kept ORFs [(a,b,strand) in strand coords], LC, LN"""
    L=len(F)
    def orfs(S,s):
        o=[]
        for f in range(3):
            last=f
            for i in range(f,len(S)-2,3):
                if S[i:i+3] in stops:
                    if i-last>=minl: o.append((last,i+3,s))
                    last=i+3
        return o
    O=orfs(F,1)+orfs(R,-1)
    seg=lambda o:(AF if o[2]==1 else AR)[o[0]:o[1]-3]
    gc=lambda o:(o[0],o[1]) if o[2]==1 else (L-o[1],L-o[0])
    if mask is not None:
        O=[o for o in O if mask[gc(o)[0]:gc(o)[1]].mean()<=0.5]
    def llr(x,LC,LN):
        i=np.arange(1,len(x)); return float((LC[i%3,x[i-1],x[i]]-LN[x[i-1],x[i]]).sum())
    sel=[o for o in O if o[1]-o[0]>=seedl]
    gcg=(F.count('G')+F.count('C'))/max(L,1)
    if gcg>0.55:
        def g3(o):
            x=seg(o); n=len(x)//3*3; y=x[:n].reshape(-1,3); g=((y==1)|(y==2)).mean(0); return (g[1]<g[0] and g[1]<g[2]) if SEED=='gc2' else (g[2]>g[0] and g[2]>g[1])
        sel=[o for o in sel if g3(o)]
    kept=[]; LC=LN=None
    for it in range(iters):
        C=np.ones((3,4,4)); N=np.ones((4,4)); m=np.ones(L,bool)
        for o in sel:
            x=seg(o); i=np.arange(1,len(x)); np.add.at(C,(i%3,x[i-1],x[i]),1); a,b=gc(o); m[a:b]=False
        idx=np.flatnonzero(m[1:]&m[:-1])+1; np.add.at(N,(AF[idx-1],AF[idx]),1)
        C/=C.sum(2,keepdims=True); N/=N.sum(1,keepdims=True); LC,LN=np.log(C),np.log(N)
        if best_start:   # L2973 (v2.27): best-start scoring for the final gene calls of every non-mitochondrial genome;
                         # code detection keeps whole-ORF scoring (v2.26 used it there too and lost Hodgkinia's code call)
            # L2969 (v2.26): GC-rich genomes have long stop-to-stop ORFs whose non-coding upstream part drags the score
            # below zero; score each ORF from its best start codon (max suffix sum) and trim it to that start.
            sc=[]
            for o in O:
                x=seg(o); i=np.arange(1,len(x)); c=LC[i%3,x[i-1],x[i]]-LN[x[i-1],x[i]]
                suf=np.concatenate([np.cumsum(c[::-1])[::-1],[0.0]]); S_=F if o[2]==1 else R; best=None
                for j in range(0,o[1]-o[0]-minl+1,3):
                    if S_[o[0]+j:o[0]+j+3] in ('ATG','GTG','TTG') and (best is None or suf[j]>best[0]): best=(suf[j],j)
                sc.append((float(suf[0]),o) if best is None else (float(best[0]),(o[0]+best[1],o[1],o[2])))
            sc.sort(key=lambda t:-t[0])
        else: sc=sorted(((llr(seg(o),LC,LN),o) for o in O),key=lambda t:-t[0])
        kept=[]; occ=np.zeros(L,np.int8)
        for v,o in sc:
            if v<=0: break
            a,b=gc(o)
            if (occ[a:b]==-o[2]).mean()>0.5: continue          # antisense shadow of a better call
            if (occ[a:b]==o[2]).mean()>0.5: continue           # L2724: same-strand duplicate frame of a better call
            kept.append(o); occ[a:b]=np.where(occ[a:b]==0,o[2],occ[a:b])
        sel=kept
    return kept,LC,LN,llr
def agr_is_stop(F,R,AF,AR,mask=None):
    """L2722: AGA/AGG resolver. If AGR are sense codons, treating them as stops fragments the genes:
    median length of the 10 longest self-trained genes drops to ~0.3-0.4x. Ratio > 0.6 => AGR are stops."""
    def top(stops):
        kept,LC,LN,llr=selftrain(F,R,AF,AR,stops,150,450,mask)
        l=sorted([o[1]-o[0] for o in kept],reverse=True)[:10]; return np.median(l) if l else 0
    r=top({'TAA','TAG','AGA','AGG'})/max(1,top({'TAA','TAG'}))
    return r>0.6, r
def detect_code(F,R,AF,AR):
    """L2725: genetic-code detector v2 (read-through test; permissive fallback; never silent)."""
    L=len(F); small=L<200000; ALL={'TAA','TAG','TGA','AGA','AGG'}
    kept,LC,LN,llr=selftrain(F,R,AF,AR,ALL,150 if small else 300,450 if small else 600)
    by={c:[] for c in ALL}
    for a,b,s in kept:
        S=F if s==1 else R; cont=S[b:b+150]
        if len(cont)<150: continue
        x=arr(cont); i=np.arange(1,len(x)); by[S[b-3:b]].append(float((LC[i%3,x[i-1],x[i]]-LN[x[i-1],x[i]]).mean()))
    nT=len(by['TAA'])+len(by['TAG'])
    if len(by['TGA'])>=5 and nT>=3:
        base=np.mean(by['TAA']+by['TAG']); tga=(np.mean(by['TGA'])-base)>TGA_THR
        # L2780: bootstrap 95% interval of the TGA read-through margin (genes resampled); reported, verdict unchanged
        _r=np.random.default_rng(7); _t=np.array(by['TGA']); _b=np.array(by['TAA']+by['TAG'])
        _bs=[_r.choice(_t,len(_t)).mean()-_r.choice(_b,len(_b)).mean() for _ in range(500)]; _lo,_hi=np.percentile(_bs,[2.5,97.5])
        _note=' [TGA margin %+.3f, 95%% interval %+.3f to %+.3f%s]'%(np.mean(_t)-base,_lo,_hi,'; UNCERTAIN: interval crosses the %.3f threshold'%TGA_THR+'' if _lo<TGA_THR<_hi else '')
        if _lo<TGA_THR<_hi and not small:
            # L2784 escalation: only when uncertain, redo the read-through with ORFs >=150 nt (more TGA-ended genes); a clear result decides
            k2,LC2,LN2,_l=selftrain(F,R,AF,AR,ALL,150,450); by2={c:[] for c in ALL}
            for a,b,s in k2:
                S=F if s==1 else R; cont=S[b:b+150]
                if len(cont)<150: continue
                x=arr(cont); i=np.arange(1,len(x)); by2[S[b-3:b]].append(float((LC2[i%3,x[i-1],x[i]]-LN2[x[i-1],x[i]]).mean()))
            _t2=np.array(by2['TGA']); _b2=np.array(by2['TAA']+by2['TAG'])
            if len(_t2)>=5 and len(_b2)>=3:
                _bs2=[_r.choice(_t2,len(_t2)).mean()-_r.choice(_b2,len(_b2)).mean() for _ in range(500)]; _lo2,_hi2=np.percentile(_bs2,[2.5,97.5])
                m2=np.mean(_t2)-np.mean(_b2)
                if not (_lo2<TGA_THR<_hi2):
                    tga=m2>TGA_THR
                    _note+=' -> escalated with ORFs>=150 nt: margin %+.3f, 95%% interval %+.3f to %+.3f (%d TGA-ended), now clear'%(m2,_lo2,_hi2,len(_t2))
                else: _note+=' -> escalation (ORFs>=150 nt) still uncertain: %+.3f to %+.3f; treat the code call with caution'%(_lo2,_hi2)
        ag_sense=(len(by['AGA'])+len(by['AGG']))>=max(3,0.3*nT)
        if not tga: return '11','TGA is a stop (standard bacterial/archaeal code)'+_note
        if small:
            stop,r=agr_is_stop(F,R,AF,AR); ag_sense=not stop
        return ('4','TGA = Trp, AGA/AGG sense (tables 4/3/5/9 share these stops)'+_note) if ag_sense else ('2','TGA = Trp, AGA/AGG = stop (vertebrate mitochondrial code)'+_note)
    kept,LC,LN,llr=selftrain(F,R,AF,AR,{'TAA','TAG'},90,300)
    if len(kept)<5: return None,'UNDETERMINED (too few genes) - please choose a table'
    cnt=[0,0]
    for a,b,s in kept:
        S=F if s==1 else R; n=b-a; lo=a+(n//4)//3*3; hi=a+(3*n//4)//3*3
        for i in range(lo,hi,3):
            for off in (0,1,2):
                if S[i+off:i+off+3]=='TGA': cnt[0 if off==0 else 1]+=1
    ratio=cnt[0]/max(1,cnt[1]/2)
    if ratio<=0.3: return '11','TGA is a stop'
    stop,r=agr_is_stop(F,R,AF,AR)
    return ('2','TGA = Trp, AGA/AGG = stop (fragmentation ratio %.2f)'%r) if stop else ('5','TGA = Trp, AGA/AGG = sense (fragmentation ratio %.2f; tables 5/9/3 share these stops)'%r)
def fs_all(S,stops):
    """L2736-L2739: frameshift candidate scores from one self-training pass.
    returns (A, B, P, pred): A = v2.7 rule [pair test + first switch inside a call]; B = L2739 within-call switches
    (any persistent frame change touching the call's frame, old frame left open >=60 nt); P = pair test (shared); pred = predicted genes"""
    L=len(S); R=revcomp(S); AF,AR=arr(S),arr(R)
    kept,LC,LN,_=selftrain(S,R,AF,AR,stops,300,600); A=[]; B=[]; P=[]
    for strand,Sq,Aq in ((1,S,AF),(-1,R,AR)):
        n=len(Aq); prev=np.roll(Aq,1); T=np.zeros((3,n+1))
        for f in range(3):
            kk=(np.arange(n)-f)%3; t=LC[kk,prev,Aq]-LN[prev,Aq]; t[0]=0; T[f,1:]=np.cumsum(t)
        conv=lambda c:c if strand==1 else L-c
        calls=sorted((a,b) for a,b,s2 in kept if s2==strand)
        for (a1,b1),(a2,b2) in zip(calls,calls[1:]):
            if a1%3==a2%3 or a2>=b1-15 or b1-a2>600: continue
            fA,fB=a1%3,a2%3
            d,c=max(((T[fB][b1]-T[fB][c])-(T[fA][b1]-T[fA][c]),c) for c in range(a2,b1-14,3))
            P.append((d,conv(c)))
        for a,b in calls:
            if b-a<600: continue
            f0=a%3; los=np.arange(a,b-120,15)
            if len(los)<17: continue
            sc=np.stack([T[f][los+120]-T[f][los] for f in range(3)]); bst=sc.argmax(0); top=sc.max(0); srt=np.sort(sc,0); marg=srt[2]-srt[1]
            if b-a>=900:
                for k in range(4,len(los)-8):
                    if (bst[k:k+8]!=f0).all() and len(set(bst[k:k+8]))==1 and (top[k:k+8]>3).all():
                        A.append((float(np.mean(top[k:k+8]-sc[f0][k:k+8]))+1,conv(int(los[k])+60))); break
            for k in range(8,len(los)-8):
                X=bst[k-8:k]; Y=bst[k:k+8]
                if len(set(X))==1 and len(set(Y))==1 and X[0]!=Y[0] and f0 in (X[0],Y[0]):
                    g=int(los[k])+60; x=int(X[0]); j=g+((x-g)%3); d=0
                    while j+3<=n and Sq[j:j+3] not in stops and d<3000: j+=3; d+=3
                    if d>=60: B.append((float(min(marg[k-8:k].mean(),marg[k:k+8].mean())),conv(g)))
    pred=[((a,b) if s==1 else (L-b,L-a)) for a,b,s in kept]
    return A+P,B,P,pred
def fs_rules(A,B,P,L,thr=None):
    """apply both rules; thresholds from the genome itself (rate caps) unless given"""
    if thr is None:
        top=lambda v,cap,floor:max(floor,(sorted(v,reverse=True)[int(cap*L/1e6)] if int(cap*L/1e6)<len(v) else floor))
        thr=(top([d for d,p in A],15,5.0), top([d for d,p in P],7.5,5.0), top([m for m,p in B],7.5,0.0))
    r27=[(d,p) for d,p in A if d>thr[0]]
    r28=[(d,p) for d,p in P if d>thr[1]]+[(m,p) for m,p in B if m>thr[2]]
    return r27,r28,thr
def fs_inject(F0,regions,rng,n):
    regions=[r for r in regions if r[1]-r[0]>=900]
    sites=sorted(rng.sample([rng.randrange(a+300,b-300) for a,b in regions],min(n,len(regions))))
    S=list(F0); shift=0; tp=[]
    for p in sites:
        q=p+shift
        if rng.random()<0.5: S.insert(q,S[q]); shift+=1
        else: del S[q]; shift-=1
        tp.append(q)
    return ''.join(S),np.array(tp)
def ssr(S,umax=5):
    """label-free simple-sequence-repeat finder (L2740, L2742).
    For repeat unit u (1-4 nt) a tract is a maximal stretch with S[i]==S[i+u]. The genome's own tract-length histogram decays
    geometrically in its bulk; phase-variable / variable-surface tracts form an excess tail. The last bulk decay ratio (lengths
    with >=100 tracts) is extrapolated, and l* = first length whose observed tail (>= l) is >=10x the extrapolation with >=3 tracts.
    L2742: homopolymers are split into A/T and G/C classes, each with its own tail (C. jejuni poly-G hid under the A/T tail)."""
    A=np.frombuffer(S.encode(),np.uint8); L=len(A); out=[]; thr={}
    for u in range(1,umax+1):
        classes=(('AT',b'AT'),('GC',b'GC')) if u==1 else (('',None),)
        for cname,cls in classes:
            m=(A[:-u]==A[u:])
            if cls is not None: m=m&np.isin(A[:-u],np.frombuffer(cls,np.uint8))
            d=np.diff(np.concatenate([[0],m.astype(np.int8),[0]])); st=np.flatnonzero(d==1); en=np.flatnonzero(d==-1)
            h=np.bincount(en-st+u,minlength=200).astype(float); bulk=[l for l in range(u+1,200) if h[l]>=100]
            if not bulk: thr[str(u)+cname]=999; continue
            lb=max(bulk); r=h[lb]/h[lb-1]; pred=np.zeros(200); pred[:lb+1]=h[:lb+1]
            for l in range(lb+1,200): pred[l]=pred[l-1]*r
            co=np.cumsum(h[::-1])[::-1]; cp=np.cumsum(pred[::-1])[::-1]
            lstar=next((l for l in range(lb+1,199) if co[l]>=3 and co[l]>=10*cp[l]),999); thr[str(u)+cname]=lstar
            for i,j in zip(st,en):
                l=j-i+u
                if l<lstar: continue
                unit=S[i:i+u]
                if u==2 and unit[0]==unit[1]: continue                         # really a homopolymer
                if u==3 and len(set(unit))==1: continue
                if u==4 and (unit[:2]==unit[2:] or len(set(unit))==1): continue  # really a di-/mononucleotide repeat
                if u==5 and len(set(unit))==1: continue                         # L2743: pentanucleotide unit (Neisseria Opa CTCTT)
                out.append((int(i),int(l),unit,u))
        if u==1: thr[1]=min(thr['1AT'],thr['1GC'])
        else: thr[u]=thr[str(u)]
    return out,thr
def _run(path,table='11',out='pred.gff3',mask_rrna=True):
    rid,F,feats=read_genbank(path); L=len(F); R=revcomp(F); AF,AR=arr(F),arr(R)
    if table=='auto':
        table,why=detect_code(F,R,AF,AR); print("auto genetic-code detection:",why)
        if table is None: return
        if table in ('11','4'):
            # L2827 (v2.20): second opinion = gene-length ratio. Treating a Trp codon as a stop chops genes: median self-trained ORF
            # length under table 4 / table 11 is 0.98-1.38 on 23 standard-code genomes and 1.64-1.98 on 6 TGA=Trp genomes (dev).
            _m=None; _rr=[f for f in feats if f['type']=='rRNA']
            if _rr:
                _m=np.zeros(L,bool)
                for f in _rr:
                    for a,b in f['parts']: _m[a:b]=True
            _sm=L<200000; _md={}
            for _t in ('11','4'):
                _k,_,_,_=selftrain(F,R,AF,AR,STOPS[_t],150 if _sm else 300,450 if _sm else 600,_m); _md[_t]=float(np.median([o[1]-o[0] for o in _k])) if _k else 1.0
            _r=_md['4']/max(_md['11'],1.0); _rt='4' if _r>=1.5 else '11'; LAST['ratio']=_r
            if _rt==table: print("  length-ratio test agrees: table 4/table 11 median ORF ratio %.2f (cut 1.5)"%_r)
            elif L<300000 and (F.count('G')+F.count('C'))/max(L,1)<0.25 and (_r>1.45)!=(table=='4'):
                # L2957 (v2.25): tiny AT-rich genomes (< 300 kb, GC < 0.25): the read-through test was wrong in both blind symbionts
                # (Sulcia, Nasuia; L2953) and the length ratio right; here trust the ratio with cut 1.45 and no uncertain band.
                _rt='4' if _r>1.45 else '11'
                print("  tiny AT-rich genome: length ratio %.2f overrides the read-through call (table %s -> %s; L2957 rule)"%(_r,table,_rt))
                table=_rt; LAST['disagree']=True
            elif _rt=='4' and table=='11' and L>2500000:
                # L3149 (v2.31): a genome-size prior. Known table-4 bacteria (Mycoplasma, Spiroplasma, Hodgkinia, ...) have genomes
                # <= ~2.2 Mb (Spiroplasma 0.7-2.2 Mb); development table-4 calls were all <= 1.21 Mb. Above 2.5 Mb the length ratio may not
                # override a read-through 'stop' call (blind miss L3145: Gluconobacter, 2.70 Mb, ratio 1.69).
                print("  DISAGREEMENT: length ratio %.2f says table 4, but the genome is %.2f Mb (> 2.5 Mb, larger than any known table-4 bacterium); keeping table 11 (L3149 rule)"%(_r,L/1e6))
                LAST['disagree']=True
            elif _r<1.4 or _r>1.6:
                print("  DISAGREEMENT: read-through says table %s, length ratio %.2f says table %s; using the length ratio (L2827: 29/29 on development)"%(table,_r,_rt))
                table=_rt; LAST['disagree']=True
            else: print("  length ratio %.2f is near the 1.5 cut; keeping the read-through call, UNCERTAIN"%_r)
    stops=STOPS[table]; mito=table in ('2','3','5','9') or (L<200000)
    mask=None
    if mask_rrna:
        rr=[f for f in feats if f['type']=='rRNA']
        if rr:
            mask=np.zeros(L,bool)
            for f in rr:
                for a,b in f['parts']: mask[a:b]=True
    _gc0=(F.count('G')+F.count('C'))/max(L,1)   # L2969 (v2.26): GC-rich genomes keep ORFs >= 180 nt (short-gene misses)
    kept,LC,LN,llr=selftrain(F,R,AF,AR,stops,150 if mito else (180 if _gc0>0.55 else 240),450 if mito else 600,mask,best_start=not mito)
    gcg=(F.count('G')+F.count('C'))/max(L,1)
    ideal=None
    tails=[extract(F,f)[-14:] for f in feats if f['type']=='rRNA' and '16S' in f['q'].get('product','')]
    if tails and table in ('11','4') and L>=200000:
        from collections import Counter; ideal=revcomp(Counter(tails).most_common(1)[0][0])
    def pair(segm):
        if not ideal: return 0
        best=0
        for j in range(len(segm)):
            for k in range(len(ideal)):
                run=0
                for t in range(min(len(segm)-j,len(ideal)-k)):
                    if segm[j+t]==ideal[k+t]: run+=1; best=max(best,run)
                    else: run=0
        return best
    START={'ATG':0,'GTG':-1.5,'TTG':-3} if not mito else {'ATG':0,'ATA':-1,'ATT':-2,'GTG':-2}
    preds=[]
    _ext=[]
    for a,b,s in kept:
        S=F if s==1 else R
        if not mito:   # L2993 (v2.28): search starts back to the previous in-frame stop; v2.26-v2.27 searched only downstream of the
            while a>=3 and S[a-3:a] not in stops: a-=3     # best-start trim and lost exact starts (E. coli 81.4% -> 63.0%)
        cands=[i for i in range(a,b-3,3) if S[i:i+3] in START]
        st=a if not cands else max(cands,key=lambda i:2*pair(S[max(0,i-17):i-4])+START[S[i:i+3]]+llr(arr(S[i:b-3]),LC,LN))
        preds.append((st,b,s)); _ext.append((a,b,s,cands))
    import os as _os
    # L3209/L3217 (v2.32): in the learned start re-pick, coding weight 0.75, fixed penalties GTG -4 / TTG -5, and a re-picked start must
    # leave the minimum gene length. Blind: exact starts +2.2 (L3213, 13 genomes) and +2.3 (L3229, 6 genomes); genes found unchanged or better.
    _MINL=180 if _gc0>0.55 else 240; _LW=0.75; _PEN={'GTG':-4.0,'TTG':-5.0}
    if not mito and preds:   # L3001 (v2.29): on by default
        # L3001 (G11c/d): learn the genome's own start signal: log-odds of upstream bases (positions -20..-1) and of the start codon,
        # chosen starts of confident genes vs the alternative in-frame start codons of the same ORFs; then re-pick every start.
        W=30; MW=1.0; NP=2; SDW=2.0   # L3009 (v2.30): window 30 nt reaches archaeal promoter elements of leaderless genes
        B='ACGT'
        def win(S,i): return S[i-W:i] if i>=W else None
        gl={}
        for (st,b,s),(a,b2,s2,cands) in zip(preds,_ext):
            S=F if s==1 else R; gl[(b,s)]=llr(arr(S[st:b-3]),LC,LN)
        thr=np.median(list(gl.values()))
        cur=list(preds)
        for _p in range(NP):
            fg=np.ones((W,4)); bg=np.ones((W,4)); cf={}; cb={}
            for (st,b,s),(a,b2,s2,cands) in zip(cur,_ext):
                if gl[(b,s)]<thr or len(cands)<2: continue
                S=F if s==1 else R
                for i in cands:
                    w=win(S,i)
                    if w is None or any(c not in B for c in w): continue
                    M=fg if i==st else bg; D=cf if i==st else cb
                    for k,c in enumerate(w): M[k,B.index(c)]+=1
                    D[S[i:i+3]]=D.get(S[i:i+3],0)+1
            PW=np.log(fg/fg.sum(1,keepdims=True))-np.log(bg/bg.sum(1,keepdims=True))
            tf=sum(cf.values())+3; tb=sum(cb.values())+3
            CL={c:np.log((cf.get(c,0)+1)/tf)-np.log((cb.get(c,0)+1)/tb) for c in START}
            def ms(S,i):
                w=win(S,i)
                if w is None: return 0.0
                return float(sum(PW[k,B.index(c)] for k,c in enumerate(w) if c in B))
            new=[]
            for (st,b,s),(a,b2,s2,cands) in zip(cur,_ext):
                S=F if s==1 else R
                if not cands: new.append((st,b,s)); continue
                _ok=[i for i in cands if b-i>=_MINL] or cands
                new.append((max(_ok,key=lambda i:MW*ms(S,i)+CL[S[i:i+3]]+_PEN.get(S[i:i+3],0)+SDW*pair(S[max(0,i-17):i-4])+_LW*llr(arr(S[i:b-3]),LC,LN)),b,s))
            cur=new
        preds=cur
    # L2731: per-gene fit to the genome's own period-3 model (mean log-likelihood ratio per nt)
    fits=[]
    for st,b,s in preds:
        x=(AF if s==1 else AR)[st:b-3]; i=np.arange(1,len(x)); fits.append(float((LC[i%3,x[i-1],x[i]]-LN[x[i-1],x[i]]).mean()) if len(x)>3 else 0.0)
    lo5=np.percentile(fits,5) if fits else 0; hi5=np.percentile(fits,95) if fits else 0
    with open(out,'w') as fh:
        fh.write("##gff-version 3\n")
        for k,(st,b,s) in enumerate(preds):
            ga,gb=(st,b) if s==1 else (L-b,L-st)
            tag='atypical_possible_foreign' if fits[k]<=lo5 else ('strong_codon_bias_possible_high_expression' if fits[k]>=hi5 else 'typical')
            fh.write("%s\tDCF\tCDS\t%d\t%d\t%.3f\t%s\t0\tID=dcf%04d;fit=%.3f;class=%s\n"%(rid,ga+1,gb,fits[k],'+' if s==1 else '-',k+1,fits[k],tag))
    # L2736-L2739: frameshift candidates (assembly indels / frameshifted pseudogenes / phase-variable genes).
    #   Two rules, thresholds from rate caps on this genome; each rule's sensitivity is self-estimated by planting 100 synthetic
    #   1-nt indels in this genome's own predicted genes; the rule with the higher self-estimate is reported (L2739, pre-registered pass).
    if len(kept)>=50 and L>=500000:
        import random as _r
        A,B,P,pg=fs_all(F,stops); r27,r28,thr=fs_rules(A,B,P,L)
        S1,t1=fs_inject(F,pg,_r.Random(11),100); A1,B1,P1,_=fs_all(S1,stops); s27,s28,_=fs_rules(A1,B1,P1,len(S1),thr)
        hit=lambda fl:float(np.mean([len(fl)>0 and np.min(np.abs(np.array([p for d,p in fl])-t))<=150 for t in t1]))
        e27,e28=hit(s27),hit(s28)
        use,est,name=(r28,e28,'within-call rule (GC-rich mode)') if e28>e27 else (r27,e27,'pair rule (standard mode)')
        fsc=sorted(use,reverse=True); fname=out.replace('.gff3','')+'_frameshifts.txt'
        with open(fname,'w') as fh:
            fh.write("# possible frameshifts (assembly indel / frameshifted pseudogene / phase-variable gene): score position\n")
            nT=len(t1); half=1.96*(est*(1-est)/max(nT,1))**0.5     # L2772: binomial 95% interval (19/20 graded genomes fall inside)
            fh.write("# %s chosen; self-estimated sensitivity for 1-nt indels: %.0f%% (95%% interval %.0f-%.0f%%, from %d planted indels; other rule: %.0f%%)\n"%(name,100*est,100*max(0,est-half),100*min(1,est+half),nT,100*min(e27,e28)))
            for d,p in fsc: fh.write("%.1f\t%d\n"%(d,p+1))
        print("  %d possible frameshifts -> %s (%s; self-estimated indel sensitivity %.0f%%, 95%% interval %.0f-%.0f%%)"%(len(fsc),fname,name,100*est,100*max(0,est-half),100*min(1,est+half)))
    # L2740: phase-variation (contingency-gene) scan: simple-sequence repeats in excess of this genome's own tract-length tail
    if L>=500000:
        tr,thr_ssr=ssr(F)
        near=[(a+1 if s2==1 else L-b+1,b if s2==1 else L-a,'+' if s2==1 else '-') for a,b,s2 in preds]
        rows=[]
        for i,l,unit,u in sorted(tr):
            hit=[g for g in near if (g[0]-101<=i<=g[1] if g[2]=='+' else g[0]<=i<=g[1]+100)]
            rows.append((i,unit,l,'%d-%d,%s'%hit[0] if hit else 'intergenic'))
        # L2753: a repeat class (>=5 tracts) is called switch-like only if >=40% of its tracts lie in or <=100 nt upstream of a
        # predicted gene; decaying genomes (M. leprae) accumulate repeats in dead DNA (0-22% genic) vs 50-83% in true switchers
        cname=lambda unit:(('poly-G/C' if unit in 'GC' else 'poly-A/T') if len(unit)==1 else 'unit-%d'%len(unit))
        cl={}
        for i,unit,l,gn in rows: c=cl.setdefault(cname(unit),[0,0]); c[1]+=1; c[0]+=(gn!='intergenic')
        verdict=['%s %d tracts, %.0f%% genic -> %s'%(k,v[1],100*v[0]/v[1],'SWITCH-LIKE' if v[0]>=0.4*v[1] else 'non-genic (decay-type)') for k,v in sorted(cl.items()) if v[1]>=5]
        with open(out.replace('.gff3','')+'_phasevar.txt','w') as fh:
            fh.write("# excess simple-sequence-repeat tracts (possible phase-variable genes); self-derived minimum lengths for poly-A/T, poly-G/C, unit 2, 3, 4, 5: %s\n"%('/'.join(str(thr_ssr[k]) if thr_ssr[k]<999 else 'none' for k in ('1AT','1GC',2,3,4,5))))
            fh.write("# repeat classes (>=5 tracts): %s\n"%('; '.join(verdict) if verdict else 'none'))
            fh.write("# tract_start unit length nearest_predicted_gene(start-end,strand; tract in body or <=100 nt upstream)\n")
            for i,unit,l,gn in rows: fh.write("%d\t(%s)n\t%d\t%s\n"%(i+1,unit,l,gn))
        print("  %d excess repeat tracts -> %s; classes: %s"%(len(tr),out.replace('.gff3','')+'_phasevar.txt','; '.join(verdict) if verdict else 'none'))
        # L2807 (v2.16): annotation cross-check. The verdict above stays label-free; if the file carries annotation, also report each
        # class's genic fraction against the annotation's INTACT CDS and the share of tracts inside annotated pseudogenes.
        aint=[(f['parts'][0][0],f['parts'][-1][1],f['strand']) for f in feats if f['type']=='CDS' and 'pseudo' not in f['q']]
        apse=[(f['parts'][0][0],f['parts'][-1][1]) for f in feats if f['type'] in ('CDS','gene') and 'pseudo' in f['q']]
        if aint and verdict:
            ca={}
            for i,l,unit,u in tr:
                k=cname(unit)
                if cl.get(k,[0,0])[1]<5: continue
                c=ca.setdefault(k,[0,0,0]); c[2]+=1
                c[0]+=any((a-100<=i<=b) if s2==1 else (a<=i<=b+100) for a,b,s2 in aint)
                c[1]+=any(a<=i<=b for a,b in apse)
            print("  annotation cross-check (verdict unchanged): "+'; '.join('%s %.0f%% at intact annotated CDS, %.0f%% inside pseudogenes'%(k,100*v[0]/v[2],100*v[1]/v[2]) for k,v in sorted(ca.items())))
    # L2735: genomic-island map = runs of >=5 consecutive genes whose smoothed fit z-score < -1.5 (window 8 genes)
    if len(preds)>=50:
        order=sorted(range(len(preds)),key=lambda k:((preds[k][0],preds[k][1]) if preds[k][2]==1 else (L-preds[k][1],L-preds[k][0])))
        fz=np.array([fits[k] for k in order]); fz=(fz-fz.mean())/(fz.std()+1e-12); sm=np.convolve(fz,np.ones(8)/8,mode='same')
        coords=[((preds[k][0],preds[k][1]) if preds[k][2]==1 else (L-preds[k][1],L-preds[k][0])) for k in order]
        def call(thr):
            isl=[]; i=0
            while i<len(sm):
                if sm[i]<thr:
                    j=i
                    while j<len(sm) and sm[j]<thr: j+=1
                    if j-i>=5: isl.append((coords[i][0],coords[j-1][1],j-i,float(sm[i:j].mean())))
                    i=j
                else: i+=1
            return isl
        thr=-1.5; isl=call(thr)
        if not isl: thr=-1.0; isl=call(thr)   # weaker-contrast genomes (e.g. E. coli): relaxed threshold, reported in the file
        with open(out.replace('.gff3','')+'_islands.txt','w') as fh:
            fh.write("# candidate genomic islands / prophages (threshold z<%.1f): start end n_genes mean_z (label-free, from per-gene fit)\n"%thr)
            for a,b,n,m in isl: fh.write("%d\t%d\t%d\t%.2f\n"%(a+1,b,n,m))
        print("  %d candidate genomic islands (possible prophages / horizontally acquired regions) -> %s"%(len(isl),out.replace('.gff3','')+'_islands.txt'))
    # L2678 / L2852 (v2.21): stem-loop scan with a dinucleotide (first-order Markov) null, for small genomes (mitochondria,
    # plastids, plasmids). Exact hairpins: 8-bp perfect reverse-complement stem, loop 3-15 nt. Human mtDNA: 5 distinct vs
    # 1.6 expected, including the light-strand origin O_L near 5732 (rCRS).
    if L <= 100000:
        import random as _rnd
        _cp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        def _stems(sq, k=8, lmin=3, lmax=15):
            h = []
            for i in range(len(sq) - 2 * k - lmin):
                a = sq[i:i + k]
                if 'N' in a: continue
                rc = ''.join(_cp.get(c, 'N') for c in reversed(a))
                for lp in range(lmin, lmax + 1):
                    j = i + k + lp
                    if sq[j:j + k] == rc: h.append((i + 1, lp)); break
            out = []
            for p_, l_ in h:
                if not out or p_ - out[-1][0] > 5: out.append((p_, l_))
            return out
        _S = ''.join(c if c in 'ACGT' else 'A' for c in F); _obs = _stems(_S)
        _tr = {}
        for x_, y_ in zip(_S, _S[1:]): _tr.setdefault(x_, {}).setdefault(y_, 0); _tr[x_][y_] += 1
        _rg = _rnd.Random(2852)
        def _mk(n):
            z = [_rg.choice(_S)]
            for _ in range(n - 1):
                d_ = _tr[z[-1]]; z.append(_rg.choices(list(d_), weights=list(d_.values()))[0])
            return ''.join(z)
        _null = [len(_stems(_mk(L))) for _ in range(20)]
        _ann = [(f['parts'][0][0], f['parts'][-1][1], f['q'].get('gene', f['q'].get('product', f['type']))[:12], f['type']) for f in feats if f['type'] in ('tRNA', 'rRNA', 'CDS', 'D-loop', 'misc_feature')]
        with open(out.replace('.gff3', '') + '_stems.txt', 'w') as fh:
            fh.write("# exact stem-loops (8-bp stem, loop 3-15), forward strand; Markov null mean %.2f over 20 simulations\n" % (sum(_null) / 20))
            for p_, l_ in _obs:
                w_ = [n_ + '(' + t_ + ')' for a_, b_, n_, t_ in _ann if a_ < p_ <= b_]
                fh.write("%d\t%d\t%s\n" % (p_, l_, ','.join(w_) or 'intergenic'))
        _ex = sum(_null) / 20.0
        import math as _m
        _pv = 1.0 - sum(_m.exp(-_ex) * _ex ** i_ / _m.factorial(i_) for i_ in range(len(_obs))) if _ex > 0 else 0.0
        print("  %d stem-loops (8-bp stem) vs %.1f expected under a dinucleotide null (Poisson P = %.3f%s) -> %s" % (len(_obs), _ex, _pv,
              ', enriched' if _pv < 0.05 else '', out.replace('.gff3', '') + '_stems.txt'))
    if table in ('2', '3', '5', '9'):
        print("  note (L2852): mitochondrial genes often end in an incomplete stop (T or TA, completed by polyadenylation), so 'exact stop codon' scores run low")
    # L2866/L2871 (v2.22): codon-phase strength S (median over predicted genes >= 300 nt; random DNA ~3.5).
    # Development result on 29 genomes: accuracy falls with GC, but at a given GC a higher S means higher accuracy
    # (correct% ~ -46*GC + 14.5*ln S, ln S coefficient 3 s.e. above zero). Informational only.
    try:
        _w = np.exp(2j * np.pi / 3); _Sv = []
        _cm = str.maketrans('ACGT', 'TGCA')
        for st_, b_, s_ in preds:
            x_ = F[st_:b_] if s_ == 1 else F[len(F) - b_:len(F) - st_].translate(_cm)[::-1]   # minus-strand coords are on the reverse complement (fixed L2921)
            m_ = len(x_) // 3
            if m_ < 100: continue
            a_ = np.frombuffer(x_[:3 * m_].encode(), dtype=np.uint8).reshape(m_, 3); S_ = 0.0
            for bb_ in b'ACGT':
                c_ = (a_ == bb_).sum(axis=0); f_ = c_.sum() / (3.0 * m_)
                if 0 < f_ < 1: S_ += abs(c_[0] + c_[1] * _w + c_[2] * _w * _w) ** 2 / (3 * m_ * f_ * (1 - f_))
            _Sv.append(S_)
        if _Sv:
            _gc = (F.count('G') + F.count('C')) / max(len(F), 1); _Sm = float(np.median(_Sv))
            _exp = 55.06 - 45.93 * _gc + 14.52 * np.log(_Sm)
            print("  codon-phase strength S = %.1f (random DNA ~3.5; development genomes 16-111) at GC %.2f -> rough expected precision %.0f%% (L2871 regression; guide only)" % (_Sm, _gc, min(99.0, _exp)))
    except Exception:
        pass
    # L2921 (v2.24): dispersed-repeat flag. Predicted genes whose exact sequence occurs >= 3 times (insertion elements, rRNA-linked
    # paralogs). On B. pertussis: 157/160 annotated multi-copy genes caught, 157/158 flags correct (pre-registered E1; shipped code re-checked below).
    try:
        _cm2 = str.maketrans('ACGT', 'TGCA'); _sq = []
        for st_, b_, s_ in preds:
            _sq.append(F[st_:b_] if s_ == 1 else F[len(F) - b_:len(F) - st_].translate(_cm2)[::-1])
        from collections import Counter as _Ctr
        _cn = _Ctr(_sq); _rep = [((st_ + 1, b_) if s_ == 1 else (len(F) - b_ + 1, len(F) - st_)) + (s_, _cn[q_]) for (st_, b_, s_), q_ in zip(preds, _sq) if _cn[q_] >= 3]
        if _rep:
            with open(out.replace('.gff3', '') + '_repeats.txt', 'w') as fh:
                fh.write("# predicted genes whose exact sequence occurs >= 3 times (possible insertion elements / mobile DNA): start end strand copies\n")
                for r_ in _rep: fh.write("%d\t%d\t%s\t%d\n" % (r_[0], r_[1], '+' if r_[2] == 1 else '-', r_[3]))
            print("  %d predicted genes belong to %d repeated families (>= 3 identical copies; possible mobile DNA) -> %s" % (len(_rep), len({q_ for q_ in _sq if _cn[q_] >= 3}), out.replace('.gff3', '') + '_repeats.txt'))
    except Exception:
        pass
    # L2913 (v2.23): the DNA helical rhythm. AA/TT dinucleotides recur every ~10-11 bp in most bacterial chromosomes
    # (31 genomes: median 11.05 bp, 26/31 above a dinucleotide-Markov null). Informational only.
    try:
        _s = F[:1000000]; _a = np.frombuffer(_s.encode(), dtype=np.uint8)
        _x = (((_a[:-1] == 65) & (_a[1:] == 65)) | ((_a[:-1] == 84) & (_a[1:] == 84))).astype(float)
        if len(_x) > 5000:
            _x -= _x.mean(); _C = np.array([np.dot(_x[:-d], _x[d:]) / (len(_x) - d) for d in range(1, 41)])
            _R = (_C - np.convolve(_C, np.ones(7) / 7, mode='same'))[4:37]; _d = np.arange(5, 38)
            _Ps = np.arange(8, 14.001, 0.05); _A = np.array([abs(np.sum(_R * np.exp(2j * np.pi * _d / P))) for P in _Ps])
            print("  helical rhythm: AA/TT dinucleotides repeat every %.2f bp (amplitude %.1e; bacteria typically 10.5-11.5 bp, L2913)" % (_Ps[_A.argmax()], _A.max()))
    except Exception:
        pass
    print("  per-gene fit written (column 6 / fit=); %d atypical genes flagged (possible horizontal transfer), %d strong-codon-bias genes flagged"%(sum(f<=lo5 for f in fits),sum(f>=hi5 for f in fits)))
    # L2818 (v2.18): reliability self-check. Median length of the tool's own predicted genes: 543-927 nt on 27 working genomes,
    # 246 nt on Hodgkinia (failed model: frame never found, genes chopped). Warning only; calls unchanged.
    if preds:
        _med=float(np.median([b-st for st,b,s in preds]))
        LAST['med']=_med; LAST['gc']=(F.count('G')+F.count('C'))/max(L,1)
        if _med<400: print("  RELIABILITY WARNING: median predicted gene only %d nt (< 400): gene model likely failed; code call and gene calls not trustworthy"%_med)
    print("%s: %d bp, GC %.2f, genetic code table %s -> predicted %d genes, written to %s (16S anti-SD anchor: %s)"%(rid,L,gcg,table,len(preds),out,ideal or 'none'))
    ann=[(f['parts'][0][0],f['parts'][0][1],f['strand']) for f in feats if f['type']=='CDS' and 'pseudo' not in f['q'] and len(f['parts'])==1]
    if ann:
        cov=np.zeros((2,L),bool); gm=np.zeros((2,L),bool)
        for st,b,s in preds:
            a2,b2=(st,b) if s==1 else (L-b,L-st); cov[0 if s==1 else 1,a2:b2]=True
        for a,b,s in ann: gm[0 if s==1 else 1,a:b]=True
        found=np.mean([cov[0 if s==1 else 1,a:b].mean()>=0.5 for a,b,s in ann])
        right=np.mean([gm[0 if s==1 else 1,((st,b) if s==1 else (L-b,L-st))[0]:((st,b) if s==1 else (L-b,L-st))[1]].mean()>=0.5 for st,b,s in preds]) if preds else 0
        pstop={(b if s==1 else L-b,s) for st,b,s in preds}; astop={(b if s==1 else a,s) for a,b,s in ann}
        print("  compared with the file's own annotation (%d CDS): genes found %.1f%% | predictions correct %.1f%% | exact stop codon %.1f%%"%(
            len(ann),100*found,100*right,100*len(pstop&astop)/len(astop)))

def main(path,table='11',out='pred.gff3',mask_rrna=True):
    """v2.19: normal run; if the L2818 reliability warning fires on a GC-rich genome (GC > 0.55), retry with the second-position
    seed (L2819) under tables 11 and 4 and adopt the one whose median predicted gene is >= 400 nt, if exactly one is."""
    global SEED
    import io, contextlib
    SEED='gc3'; LAST.clear(); _run(path,table,out,mask_rrna)
    if (LAST.get('med',999)>=400 and not LAST.get('disagree')) or LAST.get('gc',0)<=0.55: return
    print("  RESCUE (L2819): GC-rich genome, model failed or the two code tests disagreed; retrying with the second-codon-position seed under tables 11 and 4 ...")
    meds={}
    try:
        for t in ('11','4'):
            SEED='gc2'; LAST.clear()
            with contextlib.redirect_stdout(io.StringIO()): _run(path,t,out+'.rescue_'+t,mask_rrna)
            meds[t]=LAST.get('med',0)
        good=[t for t in meds if meds[t]>=400]
        print("  rescue medians: table 11 -> %d nt, table 4 -> %d nt"%(meds['11'],meds['4']))
        if len(good)==1:
            print("  RESCUE ADOPTED: genetic code table %s (TGA %s); rerunning and overwriting %s:"%(good[0],'= Trp' if good[0]=='4' else 'is a stop',out))
            SEED='gc2'; _run(path,good[0],out,mask_rrna)
        else: print("  rescue inconclusive: keep the warning; results above not trustworthy")
    finally:
        SEED='gc3'
        for t in ('11','4'):
            for suf in ('','_frameshifts.txt','_phasevar.txt','_islands.txt'):
                p=(out+'.rescue_'+t).replace('.gff3','')+suf if suf else out+'.rescue_'+t
                if os.path.exists(p): os.remove(p)

if __name__=="__main__":
    a=[x for x in sys.argv[1:] if x]
    if not a:                                            # Pydroid 'play' button: no arguments
        here=os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
        cands=[f for f in os.listdir(here) if f.lower().endswith(('.gb','.gbk','.gbff','.txt','.fna','.fasta','.fa'))]
        print("Genome files found here:",cands)
        path=input("Genome file name (or full path): ").strip() or (cands[0] if cands else '')
        if not os.path.isabs(path): path=os.path.join(here,path)
        t=input("Genetic code table [auto, 11 bacteria/archaea, 4 Mycoplasma, 2 vertebrate mito, 3 yeast mito, 5 insect mito, 9 echinoderm mito] (Enter = auto): ").strip() or 'auto'
        main(path,t,os.path.join(os.path.dirname(path),'pred.gff3'))
    else:
        t=a[a.index('--table')+1] if '--table' in a else 'auto'; o=a[a.index('--out')+1] if '--out' in a else 'pred.gff3'
        main(a[0],t,o)
