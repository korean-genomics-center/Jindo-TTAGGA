J='${JINDO_ROOT}'
MQD=J+'/Results/Genome_Assembly/Assembly_Evaluation/Merqury_Evaluation/FinalASM.J495799_Child.merqury.k21'
CEN=J+'/Results/Manuscript_Figures/data/centromere'
FAI=J+'/Resources/SUBMISSION_FINAL/assembly'

def qv(f):
    d={}
    for l in open(f):
        p=l.split()
        if len(p)>=4 and p[0].startswith('chr'):
            d[p[0]] = '>Q60*' if p[3]=='+inf' else 'Q%.1f'%float(p[3])
    return d
def cen(f):
    d={}
    for l in open(f):
        if l.startswith('#'): continue
        p=l.split()
        if p and p[0].startswith('chr'):
            d[p[0]]='%.2f-%.2f'%(int(p[1])/1e6,int(p[2])/1e6)
    return d
def length(f):
    d={}
    for l in open(f):
        p=l.split()
        if p and p[0].startswith('chr'): d[p[0]]=int(p[1])/1e6
    return d

Q1=qv(MQD+'/FinalASM.J495799_Child.merqury.k21.FinalASM.J495799_Child.Trio_Specific.hap1.qv')
Q2=qv(MQD+'/FinalASM.J495799_Child.merqury.k21.FinalASM.J495799_Child.Trio_Specific.hap2.qv')
C1=cen(CEN+'/Hap1.centromere.bed'); C2=cen(CEN+'/Hap2.centromere.bed')
import os
L1=length(FAI+'/Jindo1-G-TTAGGA.hap1.fasta.fai') if os.path.exists(FAI+'/Jindo1-G-TTAGGA.hap1.fasta.fai') else {}
L2=length(FAI+'/Jindo1-G-TTAGGA.hap2.fasta.fai')

order=['chr%d'%i for i in range(1,39)]+['chrX','chrY']
print('Chromosome\tLen_M\tCen_M\tQV_M\tLen_P\tCen_P\tQV_P')
for c in order:
    r=[c]
    for L,C,Q in ((L1,C1,Q1),(L2,C2,Q2)):
        r += ['%.2f'%L[c] if c in L else '-', C.get(c,'-'), Q.get(c,'-')]
    print('\t'.join(r))

# 검산: chrY가 최저 QV인가
import re
def num(v):
    m=re.match(r'Q([\d.]+)',v or ''); return float(m.group(1)) if m else None
vals=[(c,'M',num(Q1.get(c))) for c in order]+[(c,'P',num(Q2.get(c))) for c in order]
vals=[(c,h,v) for c,h,v in vals if v is not None]
vals.sort(key=lambda x:x[2])
print('\n=== QV 최저 10개 ===')
for c,h,v in vals[:10]: print(f'  {c:6s} {h}  Q{v:.1f}')
print(f'\nchrY(P) QV = {num(Q2.get("chrY"))}')
print('chrY가 전체 최저인가 :', vals[0][0]=='chrY')
