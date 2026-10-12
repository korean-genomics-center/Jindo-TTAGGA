import itertools, openpyxl
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
wb=openpyxl.load_workbook('Jindo_Supplementary_Data_v3.xlsx',read_only=True)
lin={str(r[0]):str(r[1]) for r in list(wb['S-Data 7 chrY lineages'].iter_rows(values_only=True))[2:] if r[0]}
S=[l.strip() for l in open('ros_samples.txt')]; n=len(S)
rows=[]
for line in open('ros_gt.tsv'):
    f=line.rstrip('\n').split('\t')
    g=[x.replace('|','/').split('/')[0] for x in f[4:4+n]]; g=['.' if x in ('.','') else x for x in g]
    if g.count('.')<=1: rows.append(g)
D=[[0]*n for _ in range(n)]
for g in rows:
    for i,j in itertools.combinations(range(n),2):
        if g[i]!='.' and g[j]!='.' and g[i]!=g[j]: D[i][j]+=1; D[j][i]+=1
for excl in (None,'J73478'):
    k=[i for i,s in enumerate(S) if s!=excl]
    sub=[[D[i][j] for j in k] for i in k]
    lab=fcluster(linkage(squareform(sub,checks=False),'average'),2,'maxclust')
    ids=[S[i] for i in k]
    same=all((lab[a]==lab[b])==(lin[ids[a]]==lin[ids[b]]) for a,b in itertools.combinations(range(len(ids)),2))
    print(f"ROS sites {len(rows)} | n={len(ids)} | split {sum(lab==1)}:{sum(lab==2)} | identical to Jindo-reference lineages: {same}")
