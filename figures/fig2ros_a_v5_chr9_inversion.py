#!/usr/bin/env python
"""
Fig 2 Panel A v5: chr9 inversion (Jindo TTAGGA vs ROS_Cfam_1.0).

v5 changes from v4:
  - FIX: v4 read data/synteny/hap1_vs_ROScfam.syri.out as Maternal, but that file
    links to the label-inverted SyRI run "Jindo_Mat", which actually carries chrY
    (= TTAGGA-P). v4 therefore drew the paternal inversion on the maternal track.
  - Input is now the SyRI run on the deposited sequences (Final_SyRI_v1), whose
    haplotype identity was verified by chromosome content:
        TTAGGA-M (hap1, carries chrX): main INV chr9:11,050,752-19,079,396 (8.03 Mb)
        TTAGGA-P (hap2, carries chrY): main INV chr9:15,287,239-23,012,340 (7.73 Mb)
  - Inversion coordinates are no longer hard-coded; the largest chr9 INV record is
    read from each SyRI file.
Data: SyRI records (chr9 only); ref = Jindo, qry = ROS_Cfam_1.0.
"""
import sys
from pathlib import Path

ROOT = Path("${JINDO_ROOT}")
SCRIPTS_DIR = ROOT / "Results/Manuscript_Figures/scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon, Rectangle
import numpy as np
import pandas as pd

from _config import FIG_W_DOUBLE

SYN_DIR = ROOT / "Results/Manuscript_Figures/data/synteny"
OUT_DIR = ROOT / "Results/Manuscript_Figures/output/fig2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SYRI1 = Path("${JINDO_ROOT}/Results/Genome_Assembly/Assembly_Evaluation/Alignment_Evaluation/Final_SyRI_v1/TTAGGA-M.R7orient.ROS_Cfam_1.0/syri_v3.ref_TTAGGA-M.qry_ROS_Cfam_1.0syri.out")   # TTAGGA-M, carries chrX (verified)
SYRI2 = Path("${JINDO_ROOT}/Results/Genome_Assembly/Assembly_Evaluation/Alignment_Evaluation/Final_SyRI_v1/TTAGGA-P.R7orient.ROS_Cfam_1.0/syri_v3.ref_TTAGGA-P.qry_ROS_Cfam_1.0syri.out")   # TTAGGA-P, carries chrY (verified)

HAP1_CHR9_LEN = 62_373_562
HAP2_CHR9_LEN = 66_475_277
ROS_CHR9_LEN  = 62_002_293

def read_main_inv(path, chrom="chr9"):
    best=None
    with open(path) as fh:
        for line in fh:
            p=line.rstrip("\n").split("\t")
            if len(p)<11 or p[0]!=chrom or p[10]!="INV": continue
            rs,re_,qs,qe=int(p[1]),int(p[2]),int(p[6]),int(p[7])
            if best is None or re_-rs>best["ref"][1]-best["ref"][0]:
                best={"ref":(rs,re_),"ros":(min(qs,qe),max(qs,qe)),"size_mb":round((re_-rs)/1e6,2)}
    return best

HAP1_INV = read_main_inv(SYRI1)
HAP2_INV = read_main_inv(SYRI2)

Y_HAP1, Y_ROS, Y_HAP2 = 2.5, 1.5, 0.5
BAR_H = 0.12
COLOR_SYN = "#CCCCCC"
COLOR_INV = "#D62728"
COLOR_HAP1_BAR = "#D62728"
COLOR_HAP2_BAR = "#1F77B4"
COLOR_ROS_BAR  = "#666666"
DISPLAY_MAX_MB = 27
display_max_bp = DISPLAY_MAX_MB * 1e6


def parse_chr9_alignments(path, ref_chr="chr9", types=("SYNAL", "INVAL")):
    rows = []
    with open(path) as fh:
        for line in fh:
            p = line.strip().split("\t")
            if len(p) < 12:
                continue
            if p[0] != ref_chr or p[5] != ref_chr:
                continue
            if p[10] not in types:
                continue
            try:
                rows.append({"type": p[10],
                             "ref_start": int(p[1]), "ref_end": int(p[2]),
                             "qry_start": int(p[6]), "qry_end": int(p[7])})
            except ValueError:
                continue
    return pd.DataFrame(rows)


def to_x(pos_bp, total_bp=display_max_bp):
    return pos_bp / total_bp


def _draw_chrom_bar(ax, y, length_bp, color, label):
    display_len = min(length_bp, display_max_bp)
    ax.add_patch(Rectangle((0, y - BAR_H/2), to_x(display_len), BAR_H,
                 facecolor=color, edgecolor="black", linewidth=0.5, zorder=5))
    if length_bp > display_max_bp:
        ax.text(1.005, y, "...", ha="left", va="center", fontsize=10, fontweight="bold", zorder=6)
    ax.text(-0.01, y, label, ha="right", va="center", fontsize=8, fontweight="bold", zorder=6)


def _ribbon(ax, x1s, x1e, y1, x2s, x2e, y2, color, alpha, z, inv=False):
    if inv:
        verts = np.array([[x1s, y1], [x1e, y1], [x2s, y2], [x2e, y2]])
    else:
        verts = np.array([[x1s, y1], [x1e, y1], [x2e, y2], [x2s, y2]])
    ax.add_patch(Polygon(verts, closed=True, facecolor=color, alpha=alpha, edgecolor="none", zorder=z))


def make_panel_a(fig, gs_subplot):
    hap1_aln = parse_chr9_alignments(SYRI1)
    hap2_aln = parse_chr9_alignments(SYRI2)

    ax = fig.add_subplot(gs_subplot)
    _draw_chrom_bar(ax, Y_HAP1, HAP1_CHR9_LEN, COLOR_HAP1_BAR, "TTAGGA-M chr9")
    _draw_chrom_bar(ax, Y_ROS,  ROS_CHR9_LEN,  COLOR_ROS_BAR,  "ROS_Cfam_1.0 chr9")
    _draw_chrom_bar(ax, Y_HAP2, HAP2_CHR9_LEN, COLOR_HAP2_BAR, "TTAGGA-P chr9")

    for _, r in hap1_aln[hap1_aln["type"] == "SYNAL"].iterrows():
        if r["ref_start"] > display_max_bp or r["qry_start"] > display_max_bp:
            continue
        _ribbon(ax, to_x(r["ref_start"]), to_x(r["ref_end"]), Y_HAP1 - BAR_H/2,
                to_x(r["qry_start"]), to_x(r["qry_end"]), Y_ROS + BAR_H/2,
                COLOR_SYN, 0.35, 1)
    for _, r in hap2_aln[hap2_aln["type"] == "SYNAL"].iterrows():
        if r["ref_start"] > display_max_bp or r["qry_start"] > display_max_bp:
            continue
        _ribbon(ax, to_x(r["qry_start"]), to_x(r["qry_end"]), Y_ROS - BAR_H/2,
                to_x(r["ref_start"]), to_x(r["ref_end"]), Y_HAP2 + BAR_H/2,
                COLOR_SYN, 0.35, 1)

    h1r, h1o = HAP1_INV["ref"], HAP1_INV["ros"]
    _ribbon(ax, to_x(h1r[0]), to_x(h1r[1]), Y_HAP1 - BAR_H/2,
            to_x(h1o[0]), to_x(h1o[1]), Y_ROS + BAR_H/2, COLOR_INV, 0.75, 4, inv=True)
    h2r, h2o = HAP2_INV["ref"], HAP2_INV["ros"]
    _ribbon(ax, to_x(h2o[0]), to_x(h2o[1]), Y_ROS - BAR_H/2,
            to_x(h2r[0]), to_x(h2r[1]), Y_HAP2 + BAR_H/2, COLOR_INV, 0.75, 4, inv=True)

    ax.text(to_x((h1r[0]+h1r[1])/2), (Y_HAP1+Y_ROS)/2, f"{HAP1_INV['size_mb']:.2f} Mb\ninversion",
            ha="center", va="center", fontsize=8, fontweight="bold", color=COLOR_INV,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=COLOR_INV, linewidth=0.7), zorder=10)
    ax.text(to_x((h2r[0]+h2r[1])/2), (Y_ROS+Y_HAP2)/2, f"{HAP2_INV['size_mb']:.2f} Mb\ninversion",
            ha="center", va="center", fontsize=8, fontweight="bold", color=COLOR_INV,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=COLOR_INV, linewidth=0.7), zorder=10)

    xtick_mb = list(range(0, DISPLAY_MAX_MB + 1, 5))
    ax.set_xticks([to_x(m*1e6) for m in xtick_mb])
    ax.set_xticklabels([f"{m}" for m in xtick_mb], fontsize=7)
    ax.set_xlabel("Position (Mb)", fontsize=9, fontweight="bold")
    ax.set_xlim(-0.12, 1.05)
    ax.set_ylim(0, 3.2)
    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)

    legend_elements = [
        mpatches.Patch(facecolor=COLOR_SYN, alpha=0.5, label="Syntenic alignment"),
        mpatches.Patch(facecolor=COLOR_INV, alpha=0.75, label="Inversion (>1 Mb)"),
    ]
    ax.legend(handles=legend_elements, loc="lower center", bbox_to_anchor=(0.5, 1.02),
              ncol=2, fontsize=7.5, frameon=False)
    return ax


def main():
    fig = plt.figure(figsize=(FIG_W_DOUBLE, 4.0))
    gs = gridspec.GridSpec(1, 1, figure=fig)
    make_panel_a(fig, gs[0])
    plt.savefig(OUT_DIR / "fig2ros_a_v5.pdf", dpi=600, bbox_inches="tight")
    plt.savefig(OUT_DIR / "fig2ros_a_v5.png", dpi=300, bbox_inches="tight")
    print("Saved v5", HAP1_INV, HAP2_INV)


if __name__ == "__main__":
    main()
