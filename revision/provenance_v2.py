import hashlib
J='${JINDO_ROOT}'
def get(fn, names):
    out={}; name=None; buf=[]
    for l in open(fn):
        if l.startswith('>'):
            if name in names: out[name]=''.join(buf).upper()
            if names<=set(out): break
            name=l[1:].split()[0]; buf=[]
        elif name in names: buf.append(l.strip())
    else:
        if name in names: out[name]=''.join(buf).upper()
    return out
md5=lambda s: hashlib.md5(s.encode()).hexdigest()
F1=get(f'{J}/Resources/SUBMISSION_FINAL/assembly/Jindo1-G-TTAGGA.hap1.fasta',{'chr1','chrX'})
F2=get(f'{J}/Resources/SUBMISSION_FINAL/assembly/Jindo1-G-TTAGGA.hap2.fasta',{'chr1','chrY'})
print("[1] Hap2_chrY.fa (TSPY / chrY analyses)")
Y=''.join(l.strip() for l in open(f'{J}/Results/Manuscript_Figures/data/centromere/fasta/Hap2_chrY.fa') if not l.startswith('>')).upper()
print(f"    length {len(Y):,} (final {len(F2['chrY']):,}) | identical to final hap2 chrY: {md5(Y)==md5(F2['chrY'])}")
print("[4] combined reference (cohort SNPs)")
C=get(f'{J}/Analysis/Population_20Jindo/Reanalysis_Combined/REF/jindo_male_combined.fasta',{'chr1','chrX','chrY'})
for k,ref in (('chr1',None),('chrX',F1['chrX']),('chrY',F2['chrY'])):
    if k=='chr1':
        print(f"    chr1 length {len(C['chr1']):,} | = final hap1 chr1: {md5(C['chr1'])==md5(F1['chr1'])} | = final hap2 chr1: {md5(C['chr1'])==md5(F2['chr1'])}")
    else:
        print(f"    {k} length {len(C[k]):,} | identical to final: {md5(C[k])==md5(ref)}")
