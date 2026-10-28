#!/usr/bin/env python
"""
Figure 3C v7 (no title): chrY synteny dotplot — Jindo Hap2 vs ROS_Cfam_1.0
                                                  with correct strand
                                                  + newly-resolved-region emphasis
                                                  + external legend.

Changes from v5 (fig3c_v5_chrY_synteny_strand_correct.py, 2026-05-09):
    - Refactored into make_panel_c(fig, gs_subplot, ...) function for use
      in combined multi-panel figure (Fig 3 final).
    - Standalone behaviour preserved: running this script directly produces
      fig3c_v9_chrY_synteny.pdf with same layout/data as v5.
    - All v5 features kept: corrected strand from .coords (forward 31% /
      reverse 69%), newly resolved (~14 Mb) shading, centromere strip,
      external legend, Wong 2011 colours.

Usage:
    Standalone: python fig3c_v9_chrY_synteny_refactored.py
    Imported:   from fig3c_v9_chrY_synteny_refactored import make_panel_c
                make_panel_c(fig, gs_subplot)
"""

import sys
from pathlib import Path

SCRIPTS_DIR = Path("${JINDO_ROOT}/Results/Manuscript_Figures/scripts")
sys.path.insert(0, str(SCRIPTS_DIR))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from _config import COLORS, FIG_W_DOUBLE

# ============================================================
# Constants
# ============================================================
ROOT = Path("${JINDO_ROOT}")
BED  = ROOT / "Results/Gene_Annotation/Synteny/MUMmer_hap2/chrY_synteny_with_strand_v2.bed"
OUT_DIR = ROOT / "Results/Manuscript_Figures/output/fig3"

ROS_CHRY_LEN   = 3_937_623
JINDO_CHRY_LEN = 21_255_890

CENTROMERE_START = 7_270_000
CENTROMERE_END   = 10_830_000

NEW_REGION_START = 7_270_000
NEW_REGION_END   = JINDO_CHRY_LEN

# Wong 2011 colors via _config.py
COLOR_FWD  = COLORS.get("telomere_fwd", "#0072B2")
COLOR_REV  = COLORS.get("telomere_rev", "#CC79A7")
COLOR_CENT = "#AACCEE"
COLOR_NEW  = "#FFF7C2"

# ============================================================
# Data loading (cached)
# ============================================================
_data_cache = {}

def _load_data():
    if "loaded" in _data_cache:
        return _data_cache

    print("[Panel C] Reading chrY synteny BED with strand ...")
    df = pd.read_csv(BED, sep="\t", header=None,
                     names=["qry_chr", "qry_start", "qry_end",
                            "ref_chr", "ref_start", "ref_end",
                            "strand", "pct_id"])

    df["aln_len_q"] = df["qry_end"] - df["qry_start"]
    df_large = df[df["aln_len_q"] >= 1000].copy()

    n_fwd = (df_large["strand"] == "+").sum()
    n_rev = (df_large["strand"] == "-").sum()

    _data_cache["df_large"] = df_large
    _data_cache["df_fwd"]   = df_large[df_large["strand"] == "+"]
    _data_cache["df_rev"]   = df_large[df_large["strand"] == "-"]
    _data_cache["n_fwd"]    = n_fwd
    _data_cache["n_rev"]    = n_rev
    _data_cache["loaded"]   = True

    print(f"[Panel C] Loaded {len(df_large):,} alignment blocks (>=1 kb)")
    print(f"  Forward (+): {n_fwd:,} ({n_fwd/len(df_large)*100:.0f}%)")
    print(f"  Reverse (-): {n_rev:,} ({n_rev/len(df_large)*100:.0f}%)")
    return _data_cache


# ============================================================
# Main panel function
# ============================================================
def make_panel_c(fig, gs_subplot, show_title=True, show_legend=True):
    """
    Render Panel C (chrY synteny dotplot) into the given gridspec subplot.

    Parameters
    ----------
    fig : matplotlib Figure
    gs_subplot : SubplotSpec
    show_title : bool
    show_legend : bool
        Whether to draw the legend on the right side. For combined figures
        with limited space, may be turned off and replaced with a shared
        legend elsewhere.

    Returns
    -------
    Axes
    """
    data = _load_data()
    df_fwd = data["df_fwd"]
    df_rev = data["df_rev"]
    df_large = data["df_large"]
    n_fwd = data["n_fwd"]
    n_rev = data["n_rev"]

    ax = fig.add_subplot(gs_subplot)

    # Background shading
    ax.axvspan(NEW_REGION_START/1e6, NEW_REGION_END/1e6,
               color=COLOR_NEW, alpha=0.55, zorder=0)
    ax.axvspan(CENTROMERE_START/1e6, CENTROMERE_END/1e6,
               color=COLOR_CENT, alpha=0.55, zorder=0.5)

    # Forward alignments
    for _, row in df_fwd.iterrows():
        ax.plot([row["qry_start"]/1e6, row["qry_end"]/1e6],
                [row["ref_start"]/1e6, row["ref_end"]/1e6],
                color=COLOR_FWD, linewidth=0.8, alpha=0.7,
                solid_capstyle="round", zorder=2)

    # Reverse alignments (qry forward, ref reverse direction)
    for _, row in df_rev.iterrows():
        ax.plot([row["qry_start"]/1e6, row["qry_end"]/1e6],
                [row["ref_end"]/1e6, row["ref_start"]/1e6],
                color=COLOR_REV, linewidth=0.8, alpha=0.7,
                solid_capstyle="round", zorder=2.5)

    # Annotations
    new_mid_x = (NEW_REGION_START + NEW_REGION_END) / 2 / 1e6
    ax.annotate(
        "Newly resolved sequence\n(~14 Mb: centromere +\ndistal arm; absent from\nROS_Cfam_1.0 chrY)",
        xy=(new_mid_x, ROS_CHRY_LEN / 1e6 * 0.45),
        xytext=(new_mid_x + 1, ROS_CHRY_LEN / 1e6 * 0.78),
        fontsize=6, ha="center", va="center", style="italic",
        color="#7B6F00",
        arrowprops=dict(arrowstyle="-", linewidth=0.4, color="#7B6F00",
                        connectionstyle="arc3,rad=0"),
        zorder=3,
    )

    cent_mid_x = (CENTROMERE_START + CENTROMERE_END) / 2 / 1e6
    ax.text(cent_mid_x, ROS_CHRY_LEN / 1e6 * 0.96,
            "predicted\ncentromere",
            ha="center", va="top",
            fontsize=5.5, color="#3366AA", style="italic",
            zorder=3)

    # Axes
    ax.set_xlim(0, JINDO_CHRY_LEN/1e6)
    ax.set_ylim(0, ROS_CHRY_LEN/1e6)
    ax.set_xlabel(
        f"Jindo Hap2 chrY position (Mb)  [total {JINDO_CHRY_LEN/1e6:.2f} Mb]",
        fontsize=7,
    )
    ax.set_ylabel(
        f"ROS_Cfam_1.0 chrY position (Mb)\n[total {ROS_CHRY_LEN/1e6:.2f} Mb]",
        fontsize=7,
    )
    ax.grid(linestyle=":", alpha=0.3, linewidth=0.4)
    ax.tick_params(axis="both", which="major", length=2.5, width=0.5, labelsize=6)

    # Legend (outside, right column)
    if show_legend:
        legend_handles = [
            Line2D([0], [0], color=COLOR_FWD, linewidth=1.5,
                   label=f"Forward alignment\n(n={n_fwd:,}, {n_fwd/len(df_large)*100:.0f}%)"),
            Line2D([0], [0], color=COLOR_REV, linewidth=1.5,
                   label=f"Reverse alignment\n(n={n_rev:,}, {n_rev/len(df_large)*100:.0f}%)"),
            Patch(facecolor=COLOR_CENT, edgecolor="none", alpha=0.55,
                  label="Predicted centromere"),
            Patch(facecolor=COLOR_NEW,  edgecolor="none", alpha=0.55,
                  label="Newly resolved (~14 Mb)"),
        ]
        ax.legend(
            handles=legend_handles,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            fontsize=6, ncol=1, handlelength=1.5, handleheight=0.9,
            borderaxespad=0,
            labelspacing=1.0,
        )

    return ax


# ============================================================
# Standalone main: produces fig3c_v9_chrY_synteny.pdf
# ============================================================
def main():
    fig_w = FIG_W_DOUBLE * 0.65
    fig_h = 3.2
    fig = plt.figure(figsize=(fig_w, fig_h))
    outer_gs = gridspec.GridSpec(1, 1,
                                  left=0.10, right=0.78,
                                  top=0.92, bottom=0.13,
                                  figure=fig)
    make_panel_c(fig, outer_gs[0, 0], show_title=True, show_legend=True)

    out_pdf = OUT_DIR / "fig3c_v9_chrY_synteny.pdf"
    out_png = OUT_DIR / "fig3c_v9_chrY_synteny.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, bbox_inches="tight", dpi=600)
    print(f"\nSaved: {out_pdf}")
    print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()
