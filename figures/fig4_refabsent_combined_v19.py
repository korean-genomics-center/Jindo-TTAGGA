#!/usr/bin/env python3
"""NEW Fig4 v1 - 4 panels:
  a = UpSet (3,278 absent from all 4 refs)
  b = composition + domain (3,163 = 2,203 family + 960 novel)
  c = g2334 representative annotation-absent gene
  d = callable accessibility (TTAGGA vs ROS) [fig5a_callable_v10]
Layout: a top(full) / b c (middle L|R) / d bottom(full)."""
import sys, os
sys.path.insert(0,"${JINDO_ROOT}/Results/Manuscript_Figures/scripts")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from importlib import import_module
import _config
from _config import COLORS

DATADIR="${JINDO_ROOT}/Results/Gene_Annotation/RefAbsent_Functional"
OUTDIR ="${JINDO_ROOT}/Results/Manuscript_Figures/output/fig4"
os.makedirs(OUTDIR,exist_ok=True)
G_DARK="#2C3E50";G_MID="#5D6D7E";G_LIGHT="#AEB6BF";G_PALE="#D5DBDB"
ACCENT="#D55E00";BLUE=COLORS['jindo_hap2']

# panel a from existing SNP module
fa = import_module("fig5a_callable_v10")

# ---- data for b/c/d ----
upset=[]
with open(f"{DATADIR}/figdata_upset.tsv") as f:
    for l in f: c,n=l.rstrip("\n").split("\t"); upset.append((c,int(n)))
upset.sort(key=lambda x:-x[1])
exons=[]
with open(f"{DATADIR}/figdata_g2334_exons.tsv") as f:
    for l in f: s,e=l.split(); exons.append((int(s),int(e)))
dx=[];dy=[]
with open(f"{DATADIR}/figdata_g2334_depth.tsv") as f:
    for l in f: p,d=l.split(); dx.append(int(p)); dy.append(int(d))
refs=["ros","gsd","boxer","zoey"]
ref_labels={"ros":"ROS_Cfam_1.0","gsd":"UU_Cfam_GSD","boxer":"Dog10K_Boxer","zoey":"UMICH_Zoey"}
set_totals={"ros":3387,"gsd":3514,"boxer":3799,"zoey":3669}
ALL4="boxer,gsd,ros,zoey"

def panel_b(fig, sub):
    g=sub.subgridspec(2,2,height_ratios=[2.6,1.4],width_ratios=[1.0,4.6],hspace=0.05,wspace=0.04)
    axb=fig.add_subplot(g[0,1]); n=len(upset); cs=[d[1] for d in upset]
    axb.bar(range(n),cs,color=[ACCENT if d[0]==ALL4 else G_MID for d in upset],width=0.68,edgecolor="white",linewidth=0.4)
    for x,(c,v) in enumerate(upset):
        axb.text(x,v+max(cs)*0.015,str(v),ha="center",va="bottom",fontsize=5,color=ACCENT if c==ALL4 else G_DARK,fontweight="bold" if c==ALL4 else "normal")
    axb.set_ylabel("Genes absent from\nreference set",fontsize=7)
    axb.set_xlim(-0.6,n-0.4);axb.set_ylim(0,max(cs)*1.15);axb.set_xticks([])
    for s in ["top","right"]:axb.spines[s].set_visible(False)
    # (panel label drawn globally below)
    axd=fig.add_subplot(g[1,1],sharex=axb)
    for x,(c,_) in enumerate(upset):
        mem=c.split(",")
        for yi,ref in enumerate(refs):
            y=len(refs)-1-yi;on=ref in mem
            axd.scatter(x,y,s=16,color=(ACCENT if(on and c==ALL4) else (G_DARK if on else G_PALE)),zorder=3)
        mm=[r for r in refs if r in mem]
        if len(mm)>1:
            ys=[len(refs)-1-refs.index(m) for m in mm]
            axd.plot([x,x],[min(ys),max(ys)],color=(ACCENT if c==ALL4 else G_DARK),lw=0.8,zorder=2)
    axd.set_yticks([]);axd.set_xticks([]);axd.set_xlim(-0.6,n-0.4);axd.set_ylim(-0.6,len(refs)-0.4)
    for s in ["top","right","bottom","left"]:axd.spines[s].set_visible(False)
    axl=fig.add_subplot(g[1,0],sharey=axd)
    for yi,ref in enumerate(refs):
        y=len(refs)-1-yi
        axl.text(0.98,y,ref_labels[ref],ha="right",va="center",fontsize=5.5,color=G_DARK)
        axl.text(1.0,y,"  "+str(set_totals[ref]),ha="left",va="center",fontsize=5,color=G_MID)
    axl.set_xlim(0,1.35);axl.set_ylim(-0.6,len(refs)-0.4);axl.axis("off")

def panel_c(fig, sub):
    g=sub.subgridspec(2,1,height_ratios=[1.0,1.3],hspace=0.7)
    axt=fig.add_subplot(g[0])
    fam,nov=2203,960
    axt.barh(0,fam,color=G_MID,height=0.5);axt.barh(0,nov,left=fam,color=BLUE,height=0.5)
    axt.text(fam/2,0,f"Family member\n{fam}",ha="center",va="center",color="white",fontsize=6,fontweight="bold")
    axt.text(fam+nov/2,0,f"Complete\nnovel\n{nov}",ha="center",va="center",color="white",fontsize=5.5,fontweight="bold")
    axt.set_xlim(0,3163);axt.set_ylim(-0.5,0.5);axt.set_yticks([])
    axt.set_xlabel("Reference-absent genes (n=3,163)",fontsize=6.5)
    for s in ["top","right","left"]:axt.spines[s].set_visible(False)
    axt.tick_params(labelsize=5.5)
    # (panel label drawn globally below)
    axb=fig.add_subplot(g[1])
    cats=[("No annotated domain",1880,G_LIGHT),("Other functional",222,G_MID),("Immune (Ig/TCR)",99,ACCENT),("Other (ZF/GPCR/kinase)",18,G_DARK)]
    for y,(lab,val,col) in enumerate(cats):
        axb.barh(y,val,color=col,height=0.62)
        axb.text(val+25,y,str(val),va="center",ha="left",fontsize=6,color=ACCENT if col==ACCENT else G_DARK,fontweight="bold" if col==ACCENT else "normal")
    axb.set_yticks(range(len(cats)));axb.set_yticklabels([c[0] for c in cats],fontsize=6)
    axb.invert_yaxis();axb.set_xlim(0,2100)
    axb.set_xlabel("Family members by domain (n=2,203)",fontsize=6.5)
    for s in ["top","right"]:axb.spines[s].set_visible(False)
    axb.tick_params(labelsize=5.5)

def panel_d(fig, sub):
    g=sub.subgridspec(2,1,height_ratios=[2.5,1],hspace=0.12)
    gstart,gend=8418075,8422810
    axt=fig.add_subplot(g[0])
    axt.fill_between(dx,dy,color=G_MID,linewidth=0)
    axt.set_xlim(gstart-100,gend+100);axt.set_ylim(0,max(dy)*1.1)
    axt.set_ylabel("RNA-seq depth",fontsize=6.5);axt.set_xticks([])
    for s in ["top","right","bottom"]:axt.spines[s].set_visible(False)
    # (panel label drawn globally below)
    axg=fig.add_subplot(g[1],sharex=axt)
    yy=0.5;axg.plot([gstart,gend],[yy,yy],color=G_DARK,lw=0.8,zorder=1)
    for s,e in exons: axg.add_patch(plt.Rectangle((s,yy-0.24),e-s,0.48,color=G_DARK,zorder=2))
    axg.set_xlim(gstart-100,gend+100);axg.set_ylim(0,1);axg.set_yticks([])
    axg.set_xlabel("chr12:8,418,075-8,422,810   (g2334, TFIIH p52, 13 exons; hap2 pair g2255)",fontsize=6.5)
    for s in ["top","right","left"]:axg.spines[s].set_visible(False)
    axg.plot([gstart+150,gstart+1150],[0.13,0.13],color=G_DARK,lw=1.3)
    axg.text(gstart+650,0.22,"1 kb",ha="center",va="bottom",fontsize=5.5,color=G_DARK)

# assemble (v19): panel order follows first-citation order in the text.
#   a = UpSet        (top-left)
#   b = composition  (top-right)
#   c = g2334        (middle, full width)
#   d = callable/SNP (bottom, full width)
fig=plt.figure(figsize=(_config.FIG_W_DOUBLE,7.4))
outer=GridSpec(3,1,height_ratios=[1.35,1.0,1.35],hspace=0.55,left=0.08,right=0.97,top=0.93,bottom=0.07)
# top band: a (UpSet, left wide) | b (composition, right)
top=outer[0].subgridspec(1,2,width_ratios=[1.7,1.0],wspace=0.38)
panel_b(fig, top[0])                       # UpSet       -> a
b_slot=top[1].subgridspec(2,1,height_ratios=[0.45,2.0],hspace=0.0)
panel_c(fig, b_slot[1])                    # composition -> b
panel_d(fig, outer[1])                     # g2334       -> c
_a_axes = fa.make_panel_a(fig, outer[2])   # callable/SNP-> d
# ---- unified panel labels (figure coords) ----
LBL_X_LEFT = 0.012
LBL_X_B    = 0.60
_lab=dict(fontsize=11,fontweight="bold",va="top",ha="left")
fig.text(LBL_X_LEFT, 0.965, "a", **_lab)
fig.text(LBL_X_B,    0.965, "b", **_lab)
fig.text(LBL_X_LEFT, 0.620, "c", **_lab)
fig.text(LBL_X_LEFT, 0.365, "d", **_lab)
plt.savefig(f"{OUTDIR}/fig4_refabsent_combined_v19.png",dpi=300,facecolor="white",bbox_inches="tight")
plt.savefig(f"{OUTDIR}/fig4_refabsent_combined_v19.pdf",facecolor="white",bbox_inches="tight")
print("saved fig4_refabsent_combined_v19")
