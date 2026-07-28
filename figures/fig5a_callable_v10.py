#!/usr/bin/env python3
"""Fig4 panel a (callable accessibility) as a module: make_panel_a(fig, gs_subplot).
Drop-in replacement for fig5a_TTAGGA_vs_ROS_v7.make_panel_a.
3 sub-panels inside the given gs_subplot:
  (1) total callable sequence  : TTAGGA > ROS
  (2) total PASS SNPs          : TTAGGA < ROS (same-breed bias)
  (3) per-chromosome callable gain (chr1-38, X, Y; uniform color, no labels)
"""
import sys, csv
import numpy as np
sys.path.insert(0,"${JINDO_ROOT}/Results/Manuscript_Figures/scripts")
import _config
from _config import COLORS
from matplotlib.gridspec import GridSpecFromSubplotSpec

G_MID="#5D6D7E"; G_LIGHT="#AEB6BF"; ACCENT="#D55E00"; BLUE=COLORS['jindo_hap2']
C_T2T=BLUE; C_ROS=G_MID; C_GAIN=BLUE; C_NEG=G_LIGHT

TSV="${JINDO_ROOT}/Analysis/Population_20Jindo/Reanalysis_Combined/Callable/callable_perchr.tsv"
SNP_T2T=12600556; SNP_ROS=13136401

def _load():
    chrom,t2t,ros,diff=[],[],[],[]
    with open(TSV) as fh:
        for row in csv.DictReader(fh,delimiter="\t"):
            chrom.append(row["chr"]); t2t.append(int(row["t2t_callable_bp"]))
            ros.append(int(row["ros_callable_bp"])); diff.append(int(row["diff"]))
    return chrom,t2t,ros,diff

def make_panel_a(fig, gs_subplot):
    chrom,t2t,ros,diff=_load()
    diffM=np.array(diff)/1e6
    idxY=chrom.index("chrY")
    order=[i for i,c in enumerate(chrom) if c!="chrY"]+[idxY]
    TOT_T2T=sum(t2t)/1e6; TOT_ROS=sum(ros)/1e6

    g=GridSpecFromSubplotSpec(3,1,subplot_spec=gs_subplot,
        height_ratios=[0.42,0.42,1.35],hspace=1.05)

    def hbar(ax,vals,xlabel,fmt):
        ax.barh([1,0],vals,color=[C_T2T,C_ROS],height=0.62)
        ax.set_yticks([1,0]); ax.set_yticklabels(["Jindo1-G-TTAGGA","ROS_Cfam_1.0"],fontsize=6.4)
        for y,v in zip([1,0],vals):
            ax.text(v+max(vals)*0.012,y,fmt(v),va="center",ha="left",fontsize=6.2)
        ax.set_xlim(0,max(vals)*1.20); ax.set_xlabel(xlabel,fontsize=7)
        for s in ["top","right"]: ax.spines[s].set_visible(False)

    ax1=fig.add_subplot(g[0]); hbar(ax1,[TOT_T2T,TOT_ROS],
        "Total callable sequence (Mb)",lambda v:f"{v:,.0f} Mb")
    ax2=fig.add_subplot(g[1]); hbar(ax2,[SNP_T2T,SNP_ROS],
        "Total PASS SNPs",lambda v:f"{v:,.0f}")
    ax2.ticklabel_format(style="plain",axis="x")

    ax3=fig.add_subplot(g[2])
    labels=[chrom[i].replace("chr","") for i in order]
    d=np.array([diffM[i] for i in order])
    colors=[C_NEG if diffM[i]<0 else C_GAIN for i in order]
    x=np.arange(len(order))
    ax3.bar(x,d,color=colors,width=0.72)
    ax3.axhline(0,color="black",lw=0.6)
    ax3.set_xticks(x); ax3.set_xticklabels(labels,fontsize=5.2)
    ax3.set_ylabel("Additional callable\nsequence in TTAGGA (Mb)",fontsize=7)
    ax3.set_xlabel("Chromosome",fontsize=7)
    ax3.set_ylim(min(d)-0.5,max(d)+0.8)
    for s in ["top","right"]: ax3.spines[s].set_visible(False)
    return (ax1,ax2,ax3)

if __name__=="__main__":
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    fig=plt.figure(figsize=(getattr(_config,"FIG_W_DOUBLE",7.2),5.0))
    gs=GridSpec(1,1)
    make_panel_a(fig,gs[0])
    plt.savefig("${JINDO_ROOT}/Results/Manuscript_Figures/output/fig4/fig4a_callable_v10_standalone.png",
                dpi=300,bbox_inches="tight")
    print("standalone saved")
