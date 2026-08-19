#!/usr/bin/env python
"""
Fig 2 Panel C v6: Hap1 vs Hap2 whole-genome dot-plot (intra-Jindo), autosomes only.

v2 changes from v1 (Jeju co-author feedback + restructure):
  - chrX and chrY removed from both axes (they confused readers; the X/Y
    comparison is not the point of an intra-haplotype synteny check).
    -> keep set is chr1..chr38 only.
  - Refactored into make_panel_c(fig, gs_subplot) for the combined figure.
  - Paths moved to backup base.
  - In-panel title removed (NC style; goes to caption).
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
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np

from _config import FIG_W_DOUBLE

DATA = ROOT / "Results/Manuscript_Figures/arxiv_figures/data"
PAF      = DATA / "hap2_vs_hap1.10kb.paf"
HAP1_FAI = DATA / "hap1.fasta.fai"
HAP2_FAI = DATA / "hap2.fasta.fai"
OUT_DIR  = ROOT / "Results/Manuscript_Figures/output/fig2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

COLOR_FWD = "#0072B2"
COLOR_REV = "#CC79A7"

# all chromosomes incl. chrX/chrY (Jeju feedback was to drop the '(chrX)/(chrY)'
# text from the AXIS LABELS, not to remove the chromosomes themselves)
KEEP = set([f"chr{i}" for i in range(1, 39)] + ["chrX", "chrY"])


def chrom_sort_key(c):
    s = str(c).replace("chr", "")
    if s == "X": return (1, 0)
    if s == "Y": return (2, 0)
    try:    return (0, int(s))
    except: return (3, s)


def load_fai_canonical(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            p = line.split("\t")
            if p[0] in KEEP:
                d[p[0]] = int(p[1])
    return d


def make_offsets(chr_list, fai):
    offsets = {}; cum = 0
    for c in chr_list:
        offsets[c] = cum; cum += fai[c]
    return offsets, cum


def make_panel_c(fig, gs_subplot):
    hap1_fai = load_fai_canonical(HAP1_FAI)
    hap2_fai = load_fai_canonical(HAP2_FAI)
    hap1_chr = sorted(hap1_fai.keys(), key=chrom_sort_key)
    hap2_chr = sorted(hap2_fai.keys(), key=chrom_sort_key)
    x_off, x_total = make_offsets(hap1_chr, hap1_fai)
    y_off, y_total = make_offsets(hap2_chr, hap2_fai)

    paf_cols = ["q_name","q_len","q_start","q_end","strand","r_name","r_len","r_start","r_end"]
    df = pd.read_csv(PAF, sep="\t", header=None, usecols=range(9), names=paf_cols)
    df = df[df["q_name"].isin(hap2_fai) & df["r_name"].isin(hap1_fai)].copy()
    df["x_start"] = (df["r_start"] + df["r_name"].map(x_off)) / 1e6
    df["x_end"]   = (df["r_end"]   + df["r_name"].map(x_off)) / 1e6
    df["y_start"] = (df["q_start"] + df["q_name"].map(y_off)) / 1e6
    df["y_end"]   = (df["q_end"]   + df["q_name"].map(y_off)) / 1e6
    n_fwd = (df["strand"]=="+").sum(); n_rev = (df["strand"]=="-").sum()

    ax = fig.add_subplot(gs_subplot)
    for strand, colour in [("+", COLOR_FWD), ("-", COLOR_REV)]:
        sub = df[df["strand"]==strand]
        xs = np.column_stack([sub["x_start"], sub["x_end"]]).T
        if strand=="+":
            ys = np.column_stack([sub["y_start"], sub["y_end"]]).T
        else:
            ys = np.column_stack([sub["y_end"], sub["y_start"]]).T
        ax.plot(xs, ys, color=colour, linewidth=0.5, alpha=0.85, solid_capstyle="round")

    for c in hap1_chr[1:]:
        ax.axvline(x_off[c]/1e6, color="#BBBBBB", linewidth=0.25, zorder=0)
    for c in hap2_chr[1:]:
        ax.axhline(y_off[c]/1e6, color="#BBBBBB", linewidth=0.25, zorder=0)

    def midpoints(chr_list, fai, offsets):
        return [(c, (offsets[c]+fai[c]/2)/1e6) for c in chr_list]
    x_mids = midpoints(hap1_chr, hap1_fai, x_off)
    y_mids = midpoints(hap2_chr, hap2_fai, y_off)
    # x and y ticks: same rule for both. Label every chromosome with a small
    # font; rotate x labels vertical so neighbours don't collide horizontally,
    # and keep y labels horizontal. Tick marks on all boundaries.
    ax.set_xticks([m for _,m in x_mids])
    ax.set_xticklabels([c.replace("chr","") for c,_ in x_mids], fontsize=3.5, rotation=90, va="top")
    ax.set_yticks([m for _,m in y_mids])
    ax.set_yticklabels([c.replace("chr","") for c,_ in y_mids], fontsize=3.5)
    ax.tick_params(axis="x", pad=2.5, length=2)
    ax.tick_params(axis="y", pad=2.5, length=2)
    ax.set_xlim(0, x_total/1e6); ax.set_ylim(0, y_total/1e6)
    ax.set_xlabel("Jindo1-G-TTAGGA-M (Maternal)", fontsize=8, fontweight="bold", labelpad=4)
    ax.set_ylabel("Jindo1-G-TTAGGA-P (Paternal)", fontsize=8, fontweight="bold", labelpad=4)
    ax.tick_params(axis="both", length=2.5, width=0.5)

    legend_elements = [
        Line2D([0],[0], color=COLOR_FWD, lw=2.0, label=f"Forward (+)  {n_fwd:,}"),
        Line2D([0],[0], color=COLOR_REV, lw=2.0, label=f"Reverse (\u2212)  {n_rev:,}"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(0.01,0.99),
              fontsize=6, frameon=True, edgecolor="#888888", facecolor="white",
              framealpha=0.92, handlelength=2.0)
    return ax


def main():
    fig = plt.figure(figsize=(FIG_W_DOUBLE*0.5, FIG_W_DOUBLE*0.5))
    gs = gridspec.GridSpec(1,1,figure=fig, left=0.10, right=0.97, top=0.95, bottom=0.10)
    make_panel_c(fig, gs[0])
    plt.savefig(OUT_DIR/"fig2ros_c_v6.pdf", bbox_inches="tight")
    plt.savefig(OUT_DIR/"fig2ros_c_v6.png", dpi=600, bbox_inches="tight")
    print("Saved v2 (autosomes only)")


if __name__=="__main__":
    main()
