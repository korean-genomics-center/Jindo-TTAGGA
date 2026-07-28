#!/usr/bin/env python
"""
Figure 3A v6 (MERGED b+c): chrY (Hap2; 21.3 Mb) multi-track map WITH the
CpG 5mC methylation track added as a final data track sharing the same chrY
x-axis (co-author request: merge structure panel b and methylation panel c).

Changes from v5 (fig3a_v5_chrY_map_refactored.py):
    - Added a CpG 5mC methylation track (purple) below the five repeat tracks
      and above the position axis, on the same 0-21.3 Mb chrY coordinate.
    - inner_gs 7 rows -> 8 rows (gene + 5 repeats + methylation + position).
    - Methylation data + loader merged in from fig3_methyl_v1.py.

Methylation provenance (corrected): single assembly individual (J495799 /
Baeksan) PacBio HiFi 5mC via pb-CpG-tools (aligned_bam_to_cpg_scores),
Hap2 alignment, 364,725 CpG sites. NOT the 20-individual WGS cohort.
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
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec

from _config import COLORS, FIG_W_DOUBLE

ROOT = Path("${JINDO_ROOT}")
OUT_DIR = ROOT / "Results/Manuscript_Figures/output/fig3"
GENES_GFF = OUT_DIR / "chrY_genes.gff3"
REPEATS_TSV = OUT_DIR / "chrY_repeats.tsv"
METH_TSV = ROOT / "Analysis/Population_20Jindo/Methylation/Analysis/chrY_methylation_hap2.tsv"

CHRY_LENGTH = 21_255_890
CHRY_LENGTH_MB = CHRY_LENGTH / 1e6

CENTROMERE_START = 7_270_000
CENTROMERE_END = 10_830_000

WINDOW_SIZE = 100_000
N_BINS = int(np.ceil(CHRY_LENGTH / WINDOW_SIZE))
METH_WINDOW = 100_000

COLORS_REPEAT = {
    "LINE/L1": COLORS.get("repeat_LINE",   "#0072B2"),
    "Simple":  COLORS.get("repeat_simple", "#F0E442"),
    "LTR":     COLORS.get("repeat_LTR",    "#E69F00"),
    "SINE":    COLORS.get("repeat_SINE",   "#56B4E9"),
    "DNA TE":  COLORS.get("repeat_DNA_TE", "#CC79A7"),
}
GENE_COLOR = "#5D2E8C"
METH_LINE  = "#7B3294"
METH_FILL  = "#C2A5CF"

CATEGORIES = ["LINE/L1", "Simple", "LTR", "SINE", "DNA TE"]

_data_cache = {}

def _load_methylation():
    pos = []
    meth = []
    with open(METH_TSV) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            try:
                pos.append(int(p[1]))
                meth.append(float(p[3]))
            except ValueError:
                continue
    pos = np.array(pos); meth = np.array(meth)
    nwin = int(pos.max() // METH_WINDOW) + 1
    wmean = np.full(nwin, np.nan)
    wx = (np.arange(nwin) * METH_WINDOW + METH_WINDOW / 2) / 1e6
    for w in range(nwin):
        m = (pos >= w * METH_WINDOW) & (pos < (w + 1) * METH_WINDOW)
        if m.sum() > 0:
            wmean[w] = meth[m].mean()
    return wx, wmean

def _load_data():
    if "loaded" in _data_cache:
        return _data_cache

    print("[Panel A] Reading genes ...")
    genes_df = pd.read_csv(GENES_GFF, sep="\t", comment="#",
                           names=["chrom", "src", "type", "start", "end", "score",
                                  "strand", "frame", "attr"])
    genes_df = genes_df[genes_df["type"] == "gene"].reset_index(drop=True)

    print("[Panel A] Reading repeats ...")
    rep_df = pd.read_csv(REPEATS_TSV, sep=r"\s+", header=None,
                         names=["sw", "div", "del", "ins", "chrom", "start", "end",
                                "left", "strand", "rep_name", "rep_class"]
                              + [f"x{i}" for i in range(10)],
                         engine="python", on_bad_lines="skip")
    rep_df = rep_df[["chrom", "start", "end", "rep_name", "rep_class"]].dropna()
    rep_df["start"] = pd.to_numeric(rep_df["start"], errors="coerce")
    rep_df["end"]   = pd.to_numeric(rep_df["end"], errors="coerce")
    rep_df = rep_df.dropna(subset=["start", "end"])
    rep_df["start"] = rep_df["start"].astype(int)
    rep_df["end"]   = rep_df["end"].astype(int)

    def categorize(rcls):
        rcls = str(rcls)
        if "LINE/L1" in rcls or rcls == "LINE/L1": return "LINE/L1"
        if "Simple_repeat" in rcls: return "Simple"
        if "LTR" in rcls: return "LTR"
        if "SINE" in rcls: return "SINE"
        if "DNA" in rcls: return "DNA TE"
        return "Other"

    rep_df["category"] = rep_df["rep_class"].apply(categorize)

    def density_track(df, n_bins, window):
        bins = np.zeros(n_bins)
        for _, row in df.iterrows():
            s_bin = row["start"] // window
            e_bin = min(row["end"] // window, n_bins - 1)
            for b in range(s_bin, e_bin + 1):
                overlap_s = max(row["start"], b * window)
                overlap_e = min(row["end"], (b + 1) * window)
                bins[b] += max(0, overlap_e - overlap_s)
        return bins / window

    densities = {cat: density_track(rep_df[rep_df["category"] == cat], N_BINS, WINDOW_SIZE)
                 for cat in CATEGORIES}
    bin_centers_mb = (np.arange(N_BINS) + 0.5) * WINDOW_SIZE / 1e6

    print("[Panel A] Reading methylation (single-individual HiFi 5mC) ...")
    meth_wx, meth_wmean = _load_methylation()

    _data_cache["genes_df"] = genes_df
    _data_cache["densities"] = densities
    _data_cache["bin_centers_mb"] = bin_centers_mb
    _data_cache["meth_wx"] = meth_wx
    _data_cache["meth_wmean"] = meth_wmean
    _data_cache["loaded"] = True
    print(f"[Panel A] Loaded {len(genes_df)} genes, "
          f"{sum(len(v) for v in densities.values())} bins, "
          f"{np.sum(~np.isnan(meth_wmean))} methylation windows.")
    return _data_cache

def make_panel_a(fig, gs_subplot, show_title=True):
    data = _load_data()
    genes_df = data["genes_df"]
    densities = data["densities"]
    bin_centers_mb = data["bin_centers_mb"]
    meth_wx = data["meth_wx"]
    meth_wmean = data["meth_wmean"]

    inner_gs = gridspec.GridSpecFromSubplotSpec(
        8, 1,
        subplot_spec=gs_subplot,
        hspace=0.18,
        height_ratios=[0.9, 1.3, 1.3, 1.3, 1.3, 1.3, 1.6, 0.5],
    )

    axes_list = []

    ax_genes = fig.add_subplot(inner_gs[0])
    # Individual genes are only a few kb wide; on a 21.26-Mb axis this is below
    # one rendered pixel, so isolated genes (the 11 proximal to 11 Mb) vanished
    # in v6 while the dense 11-15 Mb block remained visible. Enforce a minimum
    # drawn width so every gene is represented.
    MIN_GENE_MB = 0.045
    for _, gene in genes_df.iterrows():
        w = max((gene["end"] - gene["start"]) / 1e6, MIN_GENE_MB)
        ax_genes.add_patch(Rectangle((gene["start"]/1e6, 0.3), w, 0.4,
                                     facecolor=GENE_COLOR, edgecolor="none",
                                     zorder=3))
    ax_genes.set_xlim(0, CHRY_LENGTH_MB)
    ax_genes.set_ylim(0, 1)
    ax_genes.set_yticks([0.5])
    ax_genes.set_yticklabels([f"Genes\n(n={len(genes_df)})"], fontsize=6)
    ax_genes.set_xticks([])
    ax_genes.spines["bottom"].set_visible(False)
    ax_genes.spines["left"].set_visible(False)
    ax_genes.tick_params(axis="y", which="both", length=0)
    axes_list.append(ax_genes)

    track_axes = []
    for i, cat in enumerate(CATEGORIES, start=1):
        ax = fig.add_subplot(inner_gs[i])
        track_axes.append(ax)

        density_pct = densities[cat] * 100
        ax.fill_between(bin_centers_mb, density_pct, color=COLORS_REPEAT[cat],
                        alpha=0.75, linewidth=0)
        ax.set_xlim(0, CHRY_LENGTH_MB)
        ymax = max(density_pct) * 1.1 + 0.5
        ax.set_ylim(0, ymax)

        half_pct = max(density_pct) * 0.5
        full_pct = max(density_pct)
        if full_pct >= 10:
            ax.set_yticks([round(half_pct), round(full_pct)])
            ax.set_yticklabels([f"{round(half_pct)}%", f"{round(full_pct)}%"], fontsize=5)
        else:
            ax.set_yticks([round(half_pct, 1), round(full_pct, 1)])
            ax.set_yticklabels([f"{half_pct:.1f}%", f"{full_pct:.1f}%"], fontsize=5)

        ax.set_ylabel(cat, fontsize=6, rotation=0, ha="right", va="center", labelpad=8)
        ax.set_xticks([])
        ax.grid(axis="y", linestyle=":", alpha=0.3, linewidth=0.4)
        ax.tick_params(axis="y", which="both", length=1.5, width=0.4)
    axes_list.extend(track_axes)

    ax_meth = fig.add_subplot(inner_gs[6])
    valid = ~np.isnan(meth_wmean)
    ax_meth.fill_between(meth_wx[valid], 0, meth_wmean[valid],
                         color=METH_FILL, alpha=0.5, zorder=1)
    ax_meth.plot(meth_wx[valid], meth_wmean[valid],
                 color=METH_LINE, linewidth=1.0, zorder=2)
    ax_meth.axhline(50, color="#888888", linewidth=0.4, linestyle=":", zorder=0)
    ax_meth.set_xlim(0, CHRY_LENGTH_MB)
    ax_meth.set_ylim(0, 100)
    ax_meth.set_yticks([0, 50, 100])
    ax_meth.set_yticklabels(["0", "50", "100"], fontsize=5)
    ax_meth.set_ylabel("CpG 5mC\n(%)", fontsize=6, rotation=0, ha="right", va="center", labelpad=8)
    ax_meth.set_xticks([])
    ax_meth.grid(axis="y", linestyle=":", alpha=0.3, linewidth=0.4)
    ax_meth.tick_params(axis="y", which="both", length=1.5, width=0.4)
    ax_meth.spines["top"].set_visible(False)
    ax_meth.spines["right"].set_visible(False)
    axes_list.append(ax_meth)

    for ax in [ax_genes] + track_axes + [ax_meth]:
        ax.axvspan(CENTROMERE_START/1e6, CENTROMERE_END/1e6,
                   color="#AACCEE", alpha=0.25, zorder=0)

    ax_pos = fig.add_subplot(inner_gs[7])
    ax_pos.set_xlim(0, CHRY_LENGTH_MB)
    ax_pos.set_yticks([])
    ax_pos.set_xlabel("Position on chrY (Mb)", fontsize=7)
    ax_pos.tick_params(axis="x", which="major", length=2.5, width=0.5, labelsize=6)
    ax_pos.spines["left"].set_visible(False)
    ax_pos.axvspan(CENTROMERE_START/1e6, CENTROMERE_END/1e6,
                   color="#AACCEE", alpha=0.4, zorder=0)
    ax_pos.text((CENTROMERE_START + CENTROMERE_END) / 2 / 1e6, 0.5,
                "predicted centromere", ha="center", va="center",
                fontsize=6, color="#3366AA", style="italic")
    axes_list.append(ax_pos)

    return axes_list

def main():
    fig_w = FIG_W_DOUBLE
    fig_h = 5.2
    fig = plt.figure(figsize=(fig_w, fig_h))
    outer_gs = gridspec.GridSpec(1, 1,
                                  left=0.10, right=0.98,
                                  top=0.95, bottom=0.09,
                                  figure=fig)
    make_panel_a(fig, outer_gs[0, 0], show_title=True)

    out_pdf = OUT_DIR / "fig3a_v6_chrY_map_methyl.pdf"
    out_png = OUT_DIR / "fig3a_v6_chrY_map_methyl.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, bbox_inches="tight", dpi=600)
    print(f"\nSaved: {out_pdf}")
    print(f"Saved: {out_png}")

if __name__ == "__main__":
    main()
