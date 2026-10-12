import glob, collections, statistics as st, openpyxl
RA='${JINDO_ROOT}/Results/Gene_Annotation/RefAbsent_Functional'
J='${JINDO_ROOT}'
C=f'{RA}/tblastn_crossspecies/catalogue_v1'
F=[l.strip() for l in open(f'{C}/F2_ids.txt') if l.strip()]; Fs=set(F)
flag={}
for i,l in enumerate(open(f'{C}/catalogue_flags_v1.tsv')):
    f=l.rstrip('\n').split('\t')
    if i==0: H=f; continue
    flag[f[0]]=dict(zip(H,f))
span={}
for h in ('hap1','hap2'):
    for l in open(f'{J}/Results/Gene_Annotation/BRAKER3/{h}/braker.gff3'):
        f=l.split('\t')
        if len(f)>8 and f[2]=='mRNA':
            k=f"{h}_"+f[8].split('ID=')[1].split(';')[0]
            if k in Fs: span[k]=int(f[4])-int(f[3])+1
def best(fn):
    b={}
    for l in open(fn):
        p=l.split('\t')
        if p[0] not in b or float(p[7])>b[p[0]][1]: b[p[0]]=(p[1],float(p[7]),float(p[2]))
    return b
a=best(f'{C}/rbh_1to2.tsv'); b=best(f'{C}/rbh_2to1.tsv')
pairs={(q,s) for q,(s,_,pid) in a.items() if s in b and b[s][0]==q and pid>=90}
paired={x for p in pairs for x in p}
h1=sum(1 for g in F if g.startswith('hap1_')); h2=len(F)-h1
loci=len(pairs)+(len(F)-len(paired))
print(f"F2 gene models {len(F)} (hap1 {h1} / hap2 {h2})")
print(f"  reciprocal hap1-hap2 pairs (identity>=90%): {len(pairs)} | hap1-only {h1-len(pairs)} | hap2-only {h2-len(pairs)}")
print(f"  >> independent loci: {loci}")
sp=[span[g] for g in F if g in span]; pl=[int(flag[g]['prot_len']) for g in F]
print(f"\nlength: median genomic span {st.median(sp):.0f} bp | >=1 kb {sum(x>=1000 for x in sp)} ({100*sum(x>=1000 for x in sp)/len(sp):.1f}%) | <300 bp {sum(x<300 for x in sp)} ({100*sum(x<300 for x in sp)/len(sp):.1f}%)")
print(f"        median protein {st.median(pl):.0f} aa | >=100 aa {sum(x>=100 for x in pl)} ({100*sum(x>=100 for x in pl)/len(pl):.1f}%)")
grp=collections.Counter(flag[g]['group'] for g in F); print("group:",dict(grp))
dom=collections.defaultdict(set)
for h in ('hap1','hap2'):
    for fn in glob.glob(f'{J}/Results/Gene_Annotation/FunctionalAnnotation/InterProScan/{h}/*.tsv'):
        for l in open(fn):
            p=l.rstrip('\n').split('\t')
            if len(p)>12 and p[11] and p[11]!='-':
                k=f"{h}_{p[0]}"
                if k in Fs: dom[k].add(p[12])
nd=sum(1 for g in F if dom.get(g)); print(f"\nInterPro domain-bearing: {nd} ({100*nd/len(F):.1f}%)")
for gname in ('family','novel'):
    gs=[g for g in F if flag[g]['group']==gname]; n=sum(1 for g in gs if dom.get(g)); print(f"   {gname}: {n} / {len(gs)}")
top=collections.Counter(d for g in F for d in dom.get(g,()))
for d,n in top.most_common(12): print(f"   {n:4d}  {d[:70]}")
try:
    wb=openpyxl.load_workbook(f'{J}/Results/Manuscript_Figures/Supplementary/data_xlsx/Jindo_Supplementary_Data.xlsx',read_only=True)
    ws=wb['S-Data 1 RefAbsent genes']; rows=list(ws.iter_rows(values_only=True)); h=list(rows[0])
    gi=h.index('gene_id'); ei=h.index('expressed')
    ex={r[gi] for r in rows[1:] if r[gi] and str(r[ei]).strip().upper()=='Y'}
    ne=sum(1 for g in F if g in ex); print(f"\nblood expression (S-Data 1 'expressed'): {ne} ({100*ne/len(F):.1f}%)")
except Exception as e: print("expression: FAILED",e)
t=collections.Counter(flag[g]['tier50'] for g in F); print("\ncross-species (tblastn 50/50):",dict(t))
nr={}
for i,l in enumerate(open(f'{RA}/tblastn_crossspecies/crossspecies_per_gene_v1.tsv')):
    f=l.rstrip('\n').split('\t')
    if i: nr[f[0]]=f[-1]
nn=sum(1 for g in F if nr.get(g)); print(f"nr non-canid present: {nn} ({100*nn/len(F):.1f}%)")
