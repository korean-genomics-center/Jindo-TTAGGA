#!/usr/bin/env python
"""
Figure 1B v14: Chromosome landscape with gene-density shading.

Changes from v13:
  - Gene-density colorbars pushed further right (toward the 90-120 Mb end of
    the axis) to sit in the far right blank corner beside chr31/chr33, better
    filling the lower-right empty space. Placement only; no data change.

Carried over from v12/v13:
  - chr1 centromere visibility (min display width 0.8 Mb, centromere zorder
    above telomere).
  - chr27/chr32 reverse-complement preview (REVERSE_CHRS); set REVERSE_CHRS =
    set() to revert to ROS_Cfam_1.0 orientation.
"""
import sys
from pathlib import Path

ROOT = Path("${JINDO_ROOT}")
SCRIPTS_DIR = ROOT / "Results/Manuscript_Figures/scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

from _config import COLORS, FIG_W_DOUBLE

OUT_DIR = ROOT / "Results/Manuscript_Figures/output/fig1"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FAI_HAP1 = ROOT / "Results/Repeat_Annotation/RepeatMasker/Hap1/FinalASM.J495799_Child.Trio_Specific.hap1.fasta.masked.fai"
FAI_HAP2 = ROOT / "Results/Repeat_Annotation/RepeatMasker/Hap2/FinalASM.J495799_Child.Trio_Specific.hap2.fasta.masked.fai"
TELO_HAP1 = ROOT / "Results/Manuscript_Figures/data/telomere/final/FinalASM.hap1_telomeric_repeat_windows.tsv"
TELO_HAP2 = ROOT / "Results/Manuscript_Figures/data/telomere/final/FinalASM.hap2_telomeric_repeat_windows.tsv"
CENT_HAP1 = ROOT / "Results/Manuscript_Figures/data/centromere/Hap1.centromere.bed"
CENT_HAP2 = ROOT / "Results/Manuscript_Figures/data/centromere/Hap2.centromere.bed"
GD_HAP1 = ROOT / "Results/Manuscript_Figures/data/gene_density/hap1_genedensity_1Mb.tsv"
GD_HAP2 = ROOT / "Results/Manuscript_Figures/data/gene_density/hap2_genedensity_1Mb.tsv"

UNION_CHRS = [f"chr{i}" for i in range(1, 39)] + ["chrX", "chrY"]
TELO_THRESHOLD = 30

REVERSE_CHRS = {"chr27", "chr32"}
CENT_MIN_MB = 0.8

BACKBONE_EDGE = "#444444"
TELOMERE_COLOR = "#722F37"
CENTROMERE_COLOR = "#17BECF"

CMAP_H1 = LinearSegmentedColormap.from_list("h1", ["#FBEAEA", "#C0392B"])
CMAP_H2 = LinearSegmentedColormap.from_list("h2", ["#EAF0F7", "#2C5F8A"])
GD_VMAX = 60

_cache = {}

def _parse_fai(path):
    out = {}
    for line in open(path):
        p = line.strip().split("\t")
        if len(p) >= 2:
            out[p[0]] = int(p[1])
    return out

def _load_gd(path):
    gd = {}
    for line in open(path):
        c, w, n = line.strip().split("\t")
        gd.setdefault(c, {})[int(w)] = int(n)
    return gd

def _load():
    if "loaded" in _cache:
        return _cache
    _cache.update({
        "chrlen_h1": _parse_fai(FAI_HAP1), "chrlen_h2": _parse_fai(FAI_HAP2),
        "telo_h1": pd.read_csv(TELO_HAP1, sep="\t"), "telo_h2": pd.read_csv(TELO_HAP2, sep="\t"),
        "cent_h1": pd.read_csv(CENT_HAP1, sep="\t", comment="#",
                  names=["chrom","start","end","monomer_bp","array_bp","copy_number"]),
        "cent_h2": pd.read_csv(CENT_HAP2, sep="\t", comment="#",
                  names=["chrom","start","end","monomer_bp","array_bp","copy_number"]),
        "gd_h1": _load_gd(GD_HAP1), "gd_h2": _load_gd(GD_HAP2),
        "loaded": True,
    })
    for d in ("telo_h1","telo_h2"):
        df=_cache[d]; df["total"]=df["forward_repeat_number"]+df["reverse_repeat_number"]
    return _cache

def _detect_telomere(telo_df, chrom, chrlen):
    sub = telo_df[telo_df["id"] == chrom]
    if len(sub)==0: return None,None
    strong = sub[sub["total"]>=TELO_THRESHOLD]
    if len(strong)==0: return None,None
    pos = strong["window"].values
    left = pos[pos<1_000_000]; right = pos[pos>chrlen-1_000_000]
    return (left.min()/1e6 if len(left)>0 else None,
            right.max()/1e6 if len(right)>0 else None)

def draw_chr_row(ax, y, chrom, hap, data, height=0.6):
    chrlen = data[f"chrlen_h{hap}"].get(chrom, 0)
    if chrlen==0: return
    chrlen_mb = chrlen/1e6
    rev = chrom in REVERSE_CHRS
    cmap = CMAP_H1 if hap==1 else CMAP_H2
    gd = data[f"gd_h{hap}"].get(chrom, {})
    norm = Normalize(vmin=0, vmax=GD_VMAX)
    nwin_total = int(np.ceil(chrlen_mb))

    for w in range(nwin_total):
        x0 = (chrlen_mb - 1 - w) if rev else w
        x1 = min(x0 + 1, chrlen_mb)
        x0 = max(x0, 0.0)
        val = gd.get(w, 0)
        ax.add_patch(Rectangle((x0, y-height/2), x1-x0, height,
                     facecolor=cmap(norm(val)), edgecolor="none", linewidth=0))
    ax.add_patch(Rectangle((0, y-height/2), chrlen_mb, height,
                 facecolor="none", edgecolor=BACKBONE_EDGE, linewidth=0.4))

    cc = data[f"cent_h{hap}"][data[f"cent_h{hap}"]["chrom"]==chrom]
    if len(cc)>0:
        cs=cc.iloc[0]["start"]/1e6; ce=cc.iloc[0]["end"]/1e6
        if rev:
            cs, ce = chrlen_mb - ce, chrlen_mb - cs
        cw = max(ce-cs, CENT_MIN_MB)
        cx = min(cs, chrlen_mb - cw)
        cx = max(cx, 0.0)
        ax.add_patch(Rectangle((cx, y-height/2*1.1), cw, height*1.1,
                     facecolor=CENTROMERE_COLOR, edgecolor="none", zorder=11))
    lt, rt = _detect_telomere(data[f"telo_h{hap}"], chrom, chrlen)
    TW=0.5
    for t in (lt, rt):
        if t is not None:
            tpos = (chrlen_mb - t) if rev else t
            ax.add_patch(Rectangle((tpos-TW/2, y-height/2*1.1), TW, height*1.1,
                         facecolor=TELOMERE_COLOR, edgecolor="none", zorder=10))

def make_panel_b(fig, gs_subplot):
    data=_load()
    max_mb=max(max(data["chrlen_h1"].values()),max(data["chrlen_h2"].values()))/1e6
    BAR_HEIGHT=0.6; GAP_WITHIN=0.15; GAP_BETWEEN=0.7
    ax=fig.add_subplot(gs_subplot)
    y_positions=[]; cur=0
    for chrom in UNION_CHRS:
        y1=cur; y2=cur+BAR_HEIGHT+GAP_WITHIN
        y_positions.append((chrom,y1,y2)); cur=y2+BAR_HEIGHT+GAP_BETWEEN
    total=cur
    for chrom,y1,y2 in y_positions:
        if chrom in data["chrlen_h1"]: draw_chr_row(ax,y1,chrom,1,data,BAR_HEIGHT)
        if chrom in data["chrlen_h2"]: draw_chr_row(ax,y2,chrom,2,data,BAR_HEIGHT)
        ax.text(-1.5,(y1+y2)/2,chrom,fontsize=6,ha="right",va="center",fontweight="bold")
    ax.set_xlim(-3,max_mb+5); ax.set_ylim(-1,total+1); ax.invert_yaxis()
    ax.set_xlabel("Position (Mb)",fontsize=7)
    ax.tick_params(axis="x",length=2,width=0.4,labelsize=6)
    ax.tick_params(axis="y",length=0); ax.set_yticks([])
    for sp in ("top","right","left"): ax.spines[sp].set_visible(False)

    leg=[Rectangle((0,0),1,1,facecolor=CMAP_H1(0.7),edgecolor=BACKBONE_EDGE,label="Hap1"),
         Rectangle((0,0),1,1,facecolor=CMAP_H2(0.7),edgecolor=BACKBONE_EDGE,label="Hap2"),
         Rectangle((0,0),1,1,facecolor=CENTROMERE_COLOR,label="Centromere"),
         Rectangle((0,0),1,1,facecolor=TELOMERE_COLOR,label="Telomere (TTAGGG)n")]
    ax.legend(handles=leg,loc="upper right",frameon=False,fontsize=6,ncol=4,
              bbox_to_anchor=(1.0,1.02))
    sm1=ScalarMappable(norm=Normalize(0,GD_VMAX),cmap=CMAP_H1); sm1.set_array([])
    sm2=ScalarMappable(norm=Normalize(0,GD_VMAX),cmap=CMAP_H2); sm2.set_array([])

    y_chr31 = [(yy1+yy2)/2 for c,yy1,yy2 in y_positions if c=="chr31"][0]
    y_chr33 = [(yy1+yy2)/2 for c,yy1,yy2 in y_positions if c=="chr33"][0]
    cb_w = max_mb * 0.26
    cb_h = 0.45
    cb_x = max_mb * 0.72
    cax1=ax.inset_axes([cb_x, y_chr31, cb_w, cb_h], transform=ax.transData)
    cax2=ax.inset_axes([cb_x, y_chr33, cb_w, cb_h], transform=ax.transData)
    cb1=fig.colorbar(sm1,cax=cax1,orientation="horizontal")
    cb2=fig.colorbar(sm2,cax=cax2,orientation="horizontal")
    cb1.set_label("Hap1 gene density (genes/Mb)",fontsize=5.5)
    cb2.set_label("Hap2 gene density (genes/Mb)",fontsize=5.5)
    for cb in (cb1,cb2):
        cb.ax.tick_params(labelsize=5,length=1.5)
        cb.ax.xaxis.set_label_position("top")
    ticks=[0,10,20,30,40,50,60]
    cb1.set_ticks(ticks); cb2.set_ticks(ticks)
    return ax

def main():
    fig=plt.figure(figsize=(FIG_W_DOUBLE,8.0))
    gs=gridspec.GridSpec(1,1,figure=fig,left=0.05,right=0.97,top=0.97,bottom=0.04)
    make_panel_b(fig,gs[0,0])
    plt.savefig(OUT_DIR/"fig1b_landscape_v14.pdf",bbox_inches="tight")
    plt.savefig(OUT_DIR/"fig1b_landscape_v14.png",bbox_inches="tight",dpi=600)
    print("Saved v14 (colorbars pushed to far-right corner)")

if __name__=="__main__":
    main()
