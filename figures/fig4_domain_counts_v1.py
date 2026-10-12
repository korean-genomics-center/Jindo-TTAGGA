"""Domain categories for the refined reference-absent catalogue (F2, n=2,371).
Rules (InterPro description): Immune = Ig V-set/variable or T-cell receptor;
ZGK = zinc finger / G protein-coupled receptor / kinase; Other functional = any other
entry excluding generic immunoglobulin terms; No annotated domain = no InterPro entry."""
import glob, re, collections
J='${JINDO_ROOT}'
C=f'{J}/Results/Gene_Annotation/RefAbsent_Functional/tblastn_crossspecies/catalogue_v1'
F=[l.strip() for l in open(f'{C}/F2_ids.txt') if l.strip()]; Fs=set(F)
grp={}
for i,l in enumerate(open(f'{C}/catalogue_flags_v1.tsv')):
    f=l.rstrip('\n').split('\t')
    if i==0: H=f; continue
    grp[f[0]]=dict(zip(H,f))['group']
dom=collections.defaultdict(set)
for h in ('hap1','hap2'):
    for fn in glob.glob(f'{J}/Results/Gene_Annotation/FunctionalAnnotation/InterProScan/{h}/*.tsv'):
        for l in open(fn):
            p=l.rstrip('\n').split('\t')
            if len(p)>12 and p[11] not in ('','-') and f'{h}_{p[0]}' in Fs: dom[f'{h}_{p[0]}'].add(p[12])
IMM=re.compile(r'(Immunoglobulin.*(V-set|variable))|T[- ]cell receptor',re.I)
IGG=re.compile(r'Immunoglobulin|Ig-like|T[- ]cell receptor',re.I)
ZGK=re.compile(r'zinc finger|G protein-coupled|G-protein coupled|GPCR|kinase',re.I)
fam=[g for g in F if grp[g]=='family']; nov=[g for g in F if grp[g]=='novel']
rows=[]; cnt=collections.Counter(); multi=0
for g in fam:
    d=dom.get(g,set())
    cats=set()
    if any(IMM.search(x) for x in d): cats.add('Immune (Ig/TCR)')
    if any(ZGK.search(x) for x in d): cats.add('Other (ZF/GPCR/kinase)')
    if any(not IGG.search(x) and not ZGK.search(x) for x in d): cats.add('Other functional')
    if not d: cats.add('No annotated domain')
    if len(cats)>1: multi+=1
    for c in cats: cnt[c]+=1
    rows.append((g,';'.join(sorted(cats))))
withdom=sum(1 for g in fam if dom.get(g))
out=f'{J}/Results/Manuscript_Figures/data/figdata_domain_v2.tsv'
with open(out,'w') as o:
    o.write('category\tcount\n')
    for c in ('No annotated domain','Other functional','Immune (Ig/TCR)','Other (ZF/GPCR/kinase)'): o.write(f'{c}\t{cnt[c]}\n')
    o.write(f'__family\t{len(fam)}\n__novel\t{len(nov)}\n__with_domain\t{withdom}\n__multi_class\t{multi}\n')
print(dict(cnt),'| family',len(fam),'novel',len(nov),'| with domain',withdom,'| multi-class',multi,'| sum',sum(cnt.values()))
