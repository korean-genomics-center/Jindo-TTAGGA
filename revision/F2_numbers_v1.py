import glob, re, collections, statistics as st, openpyxl
RA='${JINDO_ROOT}/Results/Gene_Annotation/RefAbsent_Functional'; J='${JINDO_ROOT}'
C=f'{RA}/tblastn_crossspecies/catalogue_v1'
F=[l.strip() for l in open(f'{C}/F2_ids.txt') if l.strip()]; Fs=set(F)
fl={}
for i,l in enumerate(open(f'{C}/catalogue_flags_v1.tsv')):
    f=l.rstrip('\n').split('\t')
    if i==0: H=f; continue
    fl[f[0]]=dict(zip(H,f))
fam=[g for g in F if fl[g]['group']=='family']; nov=[g for g in F if fl[g]['group']=='novel']
print(f"[1] F2 {len(F)} | hap1 {sum(g.startswith('hap1_') for g in F)} hap2 {sum(g.startswith('hap2_') for g in F)} | family {len(fam)} novel {len(nov)}")
def best(fn):
    b={}
    for l in open(fn):
        p=l.split('\t')
        if float(p[6])>1e-5: continue
        if p[0] not in b or float(p[7])>b[p[0]][3]: b[p[0]]=(p[1],float(p[2]),100*int(p[3])/int(p[4]),float(p[7]))
    return b
a=best(f'{C}/rbh_1to2.tsv'); b=best(f'{C}/rbh_2to1.tsv')
ok=lambda x: x[1]>=90 and x[2]>=80
pairs={(q,v[0]) for q,v in a.items() if ok(v) and v[0] in b and b[v[0]][0]==q and ok(b[v[0]])}
h1=sum(g.startswith('hap1_') for g in F); h2=len(F)-h1
print(f"[2] RBH (>=90% id, >=80% cov) pairs {len(pairs)} | hap1-only {h1-len(pairs)} | hap2-only {h2-len(pairs)} | loci {len(pairs)+len(F)-2*len(pairs)}")
dom=collections.defaultdict(set)
for h in ('hap1','hap2'):
    for fn in glob.glob(f'{J}/Results/Gene_Annotation/FunctionalAnnotation/InterProScan/{h}/*.tsv'):
        for l in open(fn):
            p=l.rstrip('\n').split('\t')
            if len(p)>12 and p[11] not in ('','-'):
                k=f'{h}_{p[0]}'
                if k in Fs: dom[k].add(p[12])
imm=re.compile(r'(Immunoglobulin.*(V-set|variable|Variable))|T[- ]cell receptor',re.I)
immg=[g for g in F if any(imm.search(d) for d in dom.get(g,()))]
print(f"[3] domain: family {sum(1 for g in fam if dom.get(g))} with / {sum(1 for g in fam if not dom.get(g))} without | novel with domain {sum(1 for g in nov if dom.get(g))}")
print(f"    immune V-region (Ig V-set/variable or TCR): total {len(immg)} | in family {sum(1 for g in immg if g in set(fam))}")
chrom={}
for h in ('hap1','hap2'):
    for l in open(f'{J}/Results/Gene_Annotation/BRAKER3/{h}/braker.gff3'):
        f=l.split('\t')
        if len(f)>8 and f[2]=='mRNA':
            k=f"{h}_"+f[8].split('ID=')[1].split(';')[0]
            if k in Fs: chrom[k]=f[0]
print("    immune by chromosome:",collections.Counter(chrom.get(g,'?') for g in immg).most_common(8))
wb=openpyxl.load_workbook(f'{J}/Results/Manuscript_Figures/Supplementary/data_xlsx/Jindo_Supplementary_Data.xlsx',read_only=True)
rows=list(wb['S-Data 1 RefAbsent genes'].iter_rows(values_only=True)); hd=list(rows[0]); gi=hd.index('gene_id'); ei=hd.index('expressed')
ex={r[gi] for r in rows[1:] if r[gi] and str(r[ei]).strip().upper()=='Y'}
print(f"[4] expressed: {sum(g in ex for g in F)} / {len(F)} | immune expressed {sum(g in ex for g in immg)} / {len(immg)}")
cov={}
for i,l in enumerate(open(f'{RA}/tblastn_genome/tblastn_classification_v2.tsv')):
    f=l.rstrip('\n').split('\t')
    if i: cov[f[0]]=[float(x) for x in f[2:6]]
aa=[g for g in F if fl[g]['class']=='annotation-absent']
print(f"[5] annotation-absent {len(aa)} | present (>=50%) in all 4 genomes {sum(1 for g in aa if min(cov[g])>=50)} | SA:present {sum(1 for g in F if fl[g]['class']=='sequence-absent')}")
for g in ('hap1_g2334.t1','hap2_g2255.t1'):
    print(f"[6] {g}: in F2 {g in Fs} | class {fl.get(g,{}).get('class')} | len {fl.get(g,{}).get('prot_len')} | seg {fl.get(g,{}).get('seg_frac')} | tier {fl.get(g,{}).get('tier50')} | dog-genome cov {cov.get(g)}")
