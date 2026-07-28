#!/usr/bin/env python3
"""Figure 5 v5: chrY paternal lineage structure.
v5: panel b now compares THREE references — Jindo chrY full (21.26 Mb),
    Jindo chrY truncated to the proximal 4 Mb, and ROS_Cfam_1.0 chrY (3.94 Mb).
    The truncation control shows that the failure of ROS_Cfam_1.0 is not a
    length effect: a 4-Mb window of the Jindo chrY resolves the lineages better
    than the full chromosome, whereas the same-sized ROS chrY does not resolve
    them at all.
"""
import sys, gzip, itertools, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
sys.path.insert(0,"${JINDO_ROOT}/Results/Manuscript_Figures/scripts")
import _config
from _config import COLORS

G_DARK="#2C3E50"; G_MID="#5D6D7E"; G_LIGHT="#AEB6BF"
ACCENT="#D55E00"; BLUE=COLORS['jindo_hap2']
C_L1=BLUE; C_L2=ACCENT
C_FULL=BLUE; C_TRUNC="#56B4E9"; C_ROS=G_MID

H='${JINDO_ROOT}/Analysis/Population_20Jindo/chrY_Analysis/Haplogroup'
OUT='${JINDO_ROOT}/Results/Manuscript_Figures/output/fig5'
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size":7,"axes.linewidth":0.6,
    "xtick.major.width":0.6,"ytick.major.width":0.6})

def load(vcf):
    s=[]; g=[]
    with gzip.open(vcf,'rt') as f:
        for line in f:
            if line.startswith('##'): continue
            p=line.rstrip('\n').split('\t')
            if line.startswith('#CHROM'): s=p[9:]; continue
            row=[]
            for x in p[9:]:
                a=x.split(':')[0].replace('|','/').split('/')[0]
                row.append(-1 if a in ('.','') else int(a))
            g.append(row)
    return s, np.array(g)

def dmat(M,n):
    D=np.zeros((n,n))
    for i,j in itertools.combinations(range(n),2):
        a=M[:,i]; b=M[:,j]; ok=(a>=0)&(b>=0)
        D[i,j]=D[j,i]=np.sum(a[ok]!=b[ok])
    return D

def sep(D,cl,n):
    i1=[i for i in range(n) if cl[i]==1]; i2=[i for i in range(n) if cl[i]==2]
    w=np.array([D[i,j] for g in (i1,i2) for i,j in itertools.combinations(g,2)])
    b=np.array([D[i,j] for i in i1 for j in i2])
    return w,b

# 세 데이터셋
DATA=[
 ("full",  f'{H}/chrY_clean.vcf.gz',      "Jindo1-G-TTAGGA\nchrY, full\n(21.26 Mb)",  C_FULL),
 ("trunc", f'{H}/chrY_trunc4Mb.vcf.gz',   "Jindo1-G-TTAGGA\nchrY, proximal 4 Mb\n(truncation control)", C_TRUNC),
 ("ros",   f'{H}/ros_chrY_clean.vcf.gz',  "ROS_Cfam_1.0\nchrY\n(3.94 Mb)",           C_ROS),
]
res={}
for tag,vcf,label,col in DATA:
    s,M=load(vcf); n=len(s)
    D=dmat(M,n); Z=linkage(squareform(D),'average'); cl=fcluster(Z,2,criterion='maxclust')
    w,b=sep(D,cl,n)
    res[tag]=dict(s=s,M=M,D=D,Z=Z,cl=cl,w=w,b=b,
                  ratio=np.median(b)/np.median(w), label=label, col=col, nsite=M.shape[0])
    print(f"  {tag:6s}: {M.shape[0]:,} sites  within {np.median(w):.0f}  between {np.median(b):.0f}  ratio {res[tag]['ratio']:.1f}x")

FIGW=getattr(_config,"FIG_W_DOUBLE",7.2)
fig=plt.figure(figsize=(FIGW,4.1))
gs=GridSpec(1,2,width_ratios=[1.15,1.25],wspace=0.40,
            left=0.09,right=0.97,top=0.90,bottom=0.26)

# ── (a) dendrogram: full chrY ──
axA=fig.add_subplot(gs[0])
r=res['full']; s=r['s']; Z=r['Z']; cl=r['cl']; n=len(s)
lc={1:C_L1,2:C_L2}
link_cols={}
for i,m in enumerate(Z[:,:2].astype(int)):
    c1,c2=m
    col1=lc[cl[c1]] if c1<n else link_cols[c1]
    col2=lc[cl[c2]] if c2<n else link_cols[c2]
    link_cols[i+n]= col1 if col1==col2 else G_LIGHT
dendrogram(Z, labels=s, orientation='left', ax=axA,
           link_color_func=lambda k: link_cols[k], leaf_font_size=5.4)
for lbl in axA.get_ymajorticklabels():
    lbl.set_color(lc[cl[s.index(lbl.get_text())]])
axA.set_xlabel("chrY SNP differences", fontsize=7)
axA.spines[['top','right','left']].set_visible(False)
axA.tick_params(axis='x', labelsize=6)

# ── (b) 3-way boxplot ──
axB=fig.add_subplot(gs[1])
pos=[0,0.85,  2.3,3.15,  4.6,5.45]
order=['full','trunc','ros']
data=[]; cols=[]
for t in order:
    data += [res[t]['w'], res[t]['b']]
    cols += [res[t]['col'], res[t]['col']]
bp=axB.boxplot(data, positions=pos, widths=0.56, patch_artist=True,
               medianprops=dict(color='white',lw=1.1),
               flierprops=dict(marker='o',ms=1.4,mfc=G_LIGHT,mec='none'))
for k,(box,c) in enumerate(zip(bp['boxes'],cols)):
    if k%2==0:
        box.set_facecolor('white'); box.set_edgecolor(c); box.set_linewidth(1.3)
    else:
        box.set_facecolor(c); box.set_edgecolor(c); box.set_linewidth(1.0)
for e in ['whiskers','caps']:
    for ln in bp[e]: ln.set_color(G_MID); ln.set_linewidth(0.6)

axB.set_xticks([0.42, 2.72, 5.02])
axB.set_xticklabels([res[t]['label'] for t in order], fontsize=6.0)
axB.set_ylabel("Pairwise chrY SNP differences", fontsize=7)
axB.spines[['top','right']].set_visible(False)
ymax=max(res[t]['b'].max() for t in order)*1.08
axB.set_ylim(0, ymax)

axB.legend(handles=[Patch(facecolor='white',edgecolor=G_DARK,lw=1.2,label='within lineage'),
                    Patch(facecolor=G_DARK,edgecolor=G_DARK,label='between lineages')],
           loc='upper center', bbox_to_anchor=(0.5,-0.22),
           ncol=2, fontsize=6.0, frameon=False, columnspacing=1.4, handlelength=1.4)

fig.text(0.012, 0.965, "a", fontsize=11, fontweight="bold", va="top", ha="left")
fig.text(0.487, 0.965, "b", fontsize=11, fontweight="bold", va="top", ha="left")

fig.savefig(f"{OUT}/fig5_chrY_haplogroup_v5.png", dpi=300, bbox_inches='tight')
fig.savefig(f"{OUT}/fig5_chrY_haplogroup_v5.pdf", bbox_inches='tight')
print("\nsaved:", f"{OUT}/fig5_chrY_haplogroup_v5.png")
print("\n=== legend/본문용 수치 ===")
for t in order:
    r=res[t]
    print(f"  {t:6s}: sites {r['nsite']:,}  within {np.median(r['w']):.0f}  between {np.median(r['b']):.0f}  ratio {r['ratio']:.1f}x")
