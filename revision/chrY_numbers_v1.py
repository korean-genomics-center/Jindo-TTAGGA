import numpy as np, itertools, statistics as st
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
def load(gt,smp):
    S=[l.strip() for l in open(smp) if l.strip()]
    L=[l.rstrip('\n').split('\t') for l in open(gt)]
    pos=np.array([int(x[1]) for x in L]); alts=np.array([x[3] for x in L])
    G=np.array([[g.replace('|','/').split('/')[0] for g in x[4:4+len(S)]] for x in L],dtype=object); G[(G=='')|(G=='.')]='.'
    return S,pos,alts,G
def dist(G,m):
    n=G.shape[1]; D=np.zeros((n,n))
    for i,j in itertools.combinations(range(n),2):
        b=~m[:,i]&~m[:,j]; D[i,j]=D[j,i]=np.sum(G[b,i]!=G[b,j])
    return D
def clus(D): return fcluster(linkage(squareform(D,checks=False),'average'),2,'maxclust')
def run(tag,S,pos,alts,G,excl=None,boot=100,extra=None):
    keep=np.array([s!=excl for s in S]); ids=np.array(S)[keep]; G=G[:,keep]; m=G=='.'
    k=m.sum(1)<=1
    if extra is not None: k&=extra
    G,m,P=G[k],m[k],pos[k]; D=dist(G,m); lab=clus(D)
    a,b=[np.where(lab==g)[0] for g in (1,2)]
    if 'J76290' in ids and np.where(ids=='J76290')[0][0] in b: a,b=b,a
    wa=[D[i,j] for i,j in itertools.combinations(a,2)]; wb=[D[i,j] for i,j in itertools.combinations(b,2)]
    bt=[D[i,j] for i in a for j in b]; w=wa+wb
    fx=[]
    for r in range(G.shape[0]):
        va={G[r,i] for i in a if not m[r,i]}; vb={G[r,i] for i in b if not m[r,i]}
        if len(va)==1 and len(vb)==1 and va!=vb: fx.append(P[r])
    rng=np.random.default_rng(0); ref={frozenset(ids[lab==g]) for g in (1,2)}; h=0
    for _ in range(boot):
        ix=rng.integers(0,G.shape[0],G.shape[0]); l2=clus(dist(G[ix],m[ix])); h+={frozenset(ids[l2==g]) for g in (1,2)}==ref
    print(f"{tag:28s} n={len(ids)} sites {G.shape[0]:5d} split {len(a)}:{len(b)} | withinA {st.median(wa):.0f} withinB {st.median(wb):.0f} pooled {st.median(w):.0f} between {st.median(bt):.0f} ratio {st.median(bt)/st.median(w):.2f} | fixed {len(fx)} (max pos {max(fx)/1e6 if fx else 0:.2f} Mb) | boot {h}/{boot}")
    return G,m,ids,P
S,pos,alts,G=load('gt_matrix.tsv','samples.txt')
for ex in (None,'J73478'):
    t='' if ex is None else ' (-J73478)'
    run('Jindo full'+t,S,pos,alts,G,ex)
    run('Jindo full, no spanning-del'+t,S,pos,alts,G,ex,boot=0,extra=np.array(['*' not in a for a in alts]))
S2,p2,a2,G2=load('trunc4Mb_gt.tsv','trunc4Mb_samples.txt')
S3,p3,a3,G3=load('ros_gt.tsv','ros_samples.txt')
for ex in (None,'J73478'):
    t='' if ex is None else ' (-J73478)'
    run('4-Mb control'+t,S2,p2,a2,G2,ex); run('ROS_Cfam_1.0'+t,S3,p3,a3,G3,ex)
Gf,mf,ids,Pf=run('private-variant check',S,pos,alts,G,None,boot=0)
full=~mf.any(1); pv=[int(np.sum(full&(Gf[:,i]!='0')&(np.sum(Gf!='0',axis=1)==1))) for i in range(len(ids))]
print("private variants (sites called in all 20):", dict(zip(ids,pv)), "| median", st.median(pv))
