import sys, bisect, itertools, numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
CEN=(7269076,10826466)
SD=[]
for l in open('chrY_SD_v2.bed'):
    f=l.split()
    if len(f)>=3 and f[1].isdigit() and int(f[2])>int(f[1]): SD.append((int(f[1]),int(f[2])))
if not SD: sys.exit("ERROR: chrY_SD_v2.bed has no intervals - check chrY_self_v2.paf")
SD.sort(); st=[a for a,b in SD]
def in_sd(p):
    i=bisect.bisect_right(st,p)-1
    return i>=0 and SD[i][0]<=p<SD[i][1]
S=[l.strip() for l in open('samples.txt') if l.strip()]
L=[l.rstrip('\n').split('\t') for l in open('gt_matrix.tsv')]
pos=np.array([int(x[1]) for x in L])
G0=np.array([[g.replace('|','/').split('/')[0] for g in x[4:4+len(S)]] for x in L],dtype=object)
G0[(G0=='')|(G0=='.')]='.'
inSD=np.array([in_sd(p) for p in pos]); inC=(pos>=CEN[0])&(pos<=CEN[1])
print(f"SD intervals {len(SD)}, {sum(b-a for a,b in SD):,} bp | raw sites {len(pos)} | in SD {inSD.sum()} | in centromere {inC.sum()}")
def run(mask,excl,tag):
    keep=np.array([s!=excl for s in S]); ids=np.array(S)[keep]; G=G0[mask][:,keep]; m=G=='.'
    k=m.sum(1)<=1; G,m=G[k],m[k]; n=len(ids); D=np.zeros((n,n))
    for i,j in itertools.combinations(range(n),2):
        b=~m[:,i]&~m[:,j]; D[i,j]=D[j,i]=np.sum(G[b,i]!=G[b,j])
    lab=fcluster(linkage(squareform(D,checks=False),'average'),2,'maxclust')
    iu=np.triu_indices(n,1); same=lab[iu[0]]==lab[iu[1]]
    w,bt=np.median(D[iu][same]),np.median(D[iu][~same])
    print(f"{tag:34s} sites {G.shape[0]:5d} | split {sum(lab==1)}:{sum(lab==2)} | within {w:.0f} between {bt:.0f} | ratio {bt/w:.1f}")
for excl in (None,'J73478'):
    t='20' if excl is None else '19'
    run(np.ones(len(pos),bool),excl,f'all sites (n={t})')
    run(~inSD&~inC,excl,f'non-SD, non-centromere (n={t})')
