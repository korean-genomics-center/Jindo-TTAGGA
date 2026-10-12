import glob, collections
RA='${JINDO_ROOT}/Results/Gene_Annotation/RefAbsent_Functional'
R=f'{RA}/reverse_c4'
def merge(iv):
    iv=sorted(iv); tot=0; s,e=iv[0]
    for a,b in iv[1:]:
        if a<=e+1: e=max(e,b)
        else: tot+=e-s+1; s,e=a,b
    return tot+e-s+1
meta={}
for l in open(f'{R}/missing_geneids_v1.tsv'):
    if l.startswith('geneid'): continue
    f=l.rstrip('\n').split('\t'); meta[f[0]]={'sym':f[1],'nref':int(f[2])}
seq={}; cur=None
for l in open(f'{R}/missing_query_v1.faa'):
    if l.startswith('>'): cur=l[1:].strip(); seq[cur]=''
    else: seq[cur]+=l.strip()
q2g={q:q.split('|')[2] for q in seq}
plen={q:len(s.rstrip('*')) for q,s in seq.items()}
M=collections.Counter(); cur=None
for l in open(f'{R}/seg_missing_v1.txt'):
    if l.startswith('>'): cur=l[1:].strip()
    elif '-' in l: a,b=l.split('-'); M[cur]+=int(b)-int(a)+1
segf={q:M[q]/max(1,plen[q]) for q in seq}
def pres(fn):
    s=set()
    for line in open(fn):
        f=line.rstrip('\n').split('\t')
        if float(f[2])>=50 and 100*int(f[3])/int(f[4])>=50: s.add(f[0].split('|')[2])
    return s
mts=pres(f'{R}/mts_check/rev_hap1_mts50.tsv')|pres(f'{R}/mts_check/rev_hap2_mts50.tsv')
hsp=collections.defaultdict(list); ql={}
for fn in glob.glob(f'{R}/results_tbn/*_m*.tsv'):
    for line in open(fn):
        p=line.rstrip('\n').split('\t')
        if len(p)<10: continue
        ql[p[0]]=int(p[4])
        if float(p[2])>=50:
            a,b=sorted((int(p[6]),int(p[7]))); hsp[(p[0],p[1],fn.split('/')[-1][:4])].append((a,b))
gcov=collections.defaultdict(float)
for (q,s,h),iv in hsp.items(): gcov[q]=max(gcov[q],100*merge(iv)/ql[q])
Q=[q for q in seq if q2g[q] not in mts]
def classify(q): return 'seq-present (Jindo annotation miss)' if gcov[q]>=50 else 'seq-absent from Jindo genome'
def report(name,qs):
    c=collections.Counter(classify(q) for q in qs)
    c4=collections.Counter(classify(q) for q in qs if meta[q2g[q]]['nref']==4)
    n4=sum(1 for q in qs if meta[q2g[q]]['nref']==4)
    print(f"\n== {name}: {len(qs)} genes (in all 4 refs: {n4}) ==")
    for k in ('seq-present (Jindo annotation miss)','seq-absent from Jindo genome'):
        print(f"   {k:38s} {c[k]:5d} ({100*c[k]/max(1,len(qs)):.1f}%)   | in all 4 refs: {c4[k]:5d}")
print(f"missing after max_target_seqs correction: {len(Q)} (removed {len(seq)-len(Q)})")
for t in (50,100):
    print(f"  protein <{t} aa: {sum(plen[q]<t for q in Q)}   SEG>=50%: {sum(segf[q]>=0.5 for q in Q)}")
report('R0 no filter',Q)
F1=[q for q in Q if plen[q]>=50]; report('R1 protein >=50 aa',F1)
F2=[q for q in F1 if segf[q]<0.5]; report('R2 + SEG <50%',F2)
sa=[q for q in F2 if classify(q).startswith('seq-absent')]
print("\nR2 seq-absent examples (symbol, n_refs, len):")
for q in sorted(sa,key=lambda x:-meta[q2g[x]]['nref'])[:25]: print(f"   {meta[q2g[q]]['sym']:14s} {meta[q2g[q]]['nref']}  {plen[q]} aa  cov {gcov[q]:.0f}%")
with open(f'{R}/reverse_class_v1.tsv','w') as o:
    o.write('geneid\tsymbol\tn_refs\tprot_len\tseg_frac\tjindo_genome_cov\tclass\tpass_R2\n')
    for q in Q:
        g=q2g[q]; o.write(f"{g}\t{meta[g]['sym']}\t{meta[g]['nref']}\t{plen[q]}\t{segf[q]:.2f}\t{gcov[q]:.1f}\t{classify(q)}\t{int(q in set(F2))}\n")
print("\nwritten reverse_class_v1.tsv")
