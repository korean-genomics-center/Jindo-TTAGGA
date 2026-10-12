import glob, statistics as st
J='${JINDO_ROOT}'
W=1_000_000
for fn in sorted(glob.glob('depth_out/*_win1Mb_v3.tsv')):
    h,t=fn.split('/')[-1].split('_win1Mb')[0].split('_',1)
    CEN={}
    for l in open(f'{J}/Results/Manuscript_Figures/data/centromere/{h}.centromere.bed'):
        f=l.split()
        if len(f)>=3 and f[1].isdigit(): CEN[f[0]]=(int(f[1]),int(f[2]))
    m={}
    for l in open(fn):
        c,w,d,n=l.split(); 
        if int(n)>=W//2: m[(c,int(w))]=float(d)
    v=sorted(m.values()); med=st.median(v); q1,q3=v[len(v)//4],v[3*len(v)//4]
    def lab(k):
        c,w=k; cs=CEN.get(c); a=w*W; b=a+W
        if cs and a<cs[1] and b>cs[0]: return 'cen'
        if c=='chr19' and 21<=w<=24: return 'chr19-SD'
        return 'other'
    out=[k for k,x in m.items() if x<0.5*med or x>1.5*med]
    within=len(v)-len(out)
    cnt={x:sum(lab(k)==x for k in out) for x in ('cen','chr19-SD','other')}
    print(f"{h} {t:8s}: {len(v)} windows | median {med:.1f}x | IQR {q1:.1f}-{q3:.1f} | within 0.5-1.5x {within} ({100*within/len(v):.1f}%) | outside {len(out)} {cnt}")
    for k in sorted(out):
        if lab(k)=='other': print(f"      {k[0]}:{k[1]}-{k[1]+1}Mb  {m[k]:.1f}x")
