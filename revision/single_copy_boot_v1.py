import bisect, itertools, numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
CEN=(7269076,10826466)
SD=sorted((int(l.split()[1]),int(l.split()[2])) for l in open('chrY_SD_v2.bed') if len(l.split())>=3)
st=[a for a,b in SD]
def in_sd(p):
    i=bisect.bisect_right(st,p)-1; return i>=0 and SD[i][0]<=p<SD[i][1]
S=[l.strip() for l in open('samples.txt') if l.strip()]
L=[l.rstrip('\n').split('\t') for l in open('gt_matrix.tsv')]
pos=np.array([int(x[1]) for x in L])
G0=np.array([[g.replace('|','/').split('/')[0] for g in x[4:4+len(S)]] for x in L],dtype=object); G0[(G0=='')|(G0=='.')]='.'
mask=np.array([not in_sd(p) for p in pos])&~((pos>=CEN[0])&(pos<=CEN[1]))
def part(G,m):
    n=G.shape[1]; D=np.zeros((n,n))
    for i,j in itertools.combinations(range(n),2):
        b=~m[:,i]&~m[:,j]; D[i,j]=D[j,i]=np.sum(G[b,i]!=G[b,j])
    return fcluster(linkage(squareform(D,checks=False),'average'),2,'maxclust')
rng=np.random.default_rng(0)
for excl in (None,'J73478'):
    keep=np.array([s!=excl for s in S]); ids=np.array(S)[keep]
    G=G0[mask][:,keep]; m=G=='.'; k=m.sum(1)<=1; G,m=G[k],m[k]
    lab=part(G,m); ref={frozenset(ids[lab==g]) for g in (1,2)}; hit=0
    for _ in range(100):
        ix=rng.integers(0,G.shape[0],G.shape[0]); l2=part(G[ix],m[ix])
        hit+={frozenset(ids[l2==g]) for g in (1,2)}==ref
    print(f"non-SD non-cen n={len(ids)} sites {G.shape[0]} split {sum(lab==1)}:{sum(lab==2)} bootstrap {hit}/100")
