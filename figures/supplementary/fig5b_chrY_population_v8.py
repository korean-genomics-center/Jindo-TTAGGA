#!/usr/bin/env python
"""
Figure 5B v8: chrY POPULATION variation across 20 male Jindo - per-individual heatmap.

Redesign from v7 (density tracks) to make the *population* nature explicit and to
visually distinguish this panel from Fig3c (chrY synteny dot-plot).

Design:
    - Two stacked heatmap blocks: ROS_Cfam_1.0 (top) and Jindo1-G-TTAGGA Hap2 (bottom).
    - Within each block, 20 rows = 20 male Jindo individuals.
    - x-axis = chrY position; cell color = per-individual SNP count in each 100-kb bin.
    - ROS block: signal only in the first ~4 Mb (rest is unassembled -> empty).
    - TTAGGA block: signal across the full 21.26 Mb -> "population variation visible
      only with the complete reference."

Data: Results/Manuscript_Figures/data/chry_population/{t2t,ros}_chrY_per_sample.tsv
      (columns: sample<TAB>position), produced from the GATK PASS chrY SNPs.
All paths on backup.
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
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap

from _config import COLORS, FIG_W_DOUBLE

DATA_DIR = ROOT / "Results/Manuscript_Figures/data/chry_population"
OUT_DIR = ROOT / "Results/Manuscript_Figures/output/fig5"
OUT_DIR.mkdir(parents=True, exist_ok=True)

T2T_TSV = DATA_DIR / "t2t_chrY_per_sample.tsv"
ROS_TSV = DATA_DIR / "ros_chrY_per_sample.tsv"

CHRY_LENGTH = 21_255_906
CHRY_LENGTH_MB = CHRY_LENGTH / 1e6
CENTROMERE_START_MB = 7.27
CENTROMERE_END_MB = 10.83
WINDOW_SIZE = 100_000
N_BINS = int(np.ceil(CHRY_LENGTH / WINDOW_SIZE))
ROS_COVERAGE_END_MB = 4.0

# Fixed sample order (matches VCF header)
SAMPLES = ["J12952315","J13044625","J288326","J288508","J352388","J433555",
           "J495439","J496128","J498404","J543672","J61435","J61783","J62428",
           "J62957","J63253","J63542","J73201","J73478","J76290","J90474"]

COLOR_T2T = "#0072B2"  # TTAGGA chrY heatmap = blue (single track, not Hap1/Hap2 contrast)
COLOR_ROS = COLORS.get("ros",        "#5D6D7E")  # ROS = gray


def _build_matrix(tsv_path):
    """Return 20 x N_BINS matrix of per-sample SNP counts per 100-kb bin."""
    idx = {s: i for i, s in enumerate(SAMPLES)}
    mat = np.zeros((len(SAMPLES), N_BINS))
    with open(tsv_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            s, pos = parts[0], parts[1]
            if s not in idx:
                continue
            b = min(int(pos) // WINDOW_SIZE, N_BINS - 1)
            mat[idx[s], b] += 1
    return mat


def _cmap(base_hex):
    return LinearSegmentedColormap.from_list("c", ["white", base_hex])


def make_panel_b(fig, gs_subplot):
    t2t_mat = _build_matrix(T2T_TSV)
    ros_mat = _build_matrix(ROS_TSV)
    print(f"[Panel B] T2T total SNP records: {int(t2t_mat.sum()):,}")
    print(f"[Panel B] ROS total SNP records: {int(ros_mat.sum()):,}")

    # shared color scale (95th percentile of nonzero to avoid outlier washout)
    vmax = np.percentile(np.concatenate([t2t_mat[t2t_mat > 0], ros_mat[ros_mat > 0]]), 95)

    inner_gs = gridspec.GridSpecFromSubplotSpec(
        3, 1, subplot_spec=gs_subplot, hspace=0.18, height_ratios=[1.0, 1.0, 0.18])

    extent = [0, CHRY_LENGTH_MB, 0, len(SAMPLES)]

    # ---- ROS block ----
    ax_ros = fig.add_subplot(inner_gs[0])
    ax_ros.imshow(ros_mat, aspect="auto", cmap=_cmap(COLOR_ROS),
                  vmin=0, vmax=vmax, extent=extent, origin="lower", interpolation="nearest")
    ax_ros.axvspan(ROS_COVERAGE_END_MB, CHRY_LENGTH_MB, facecolor="lightgray",
                   alpha=0.30, hatch="///", edgecolor="gray", linewidth=0.3, zorder=3)
    ax_ros.set_xlim(0, CHRY_LENGTH_MB)
    ax_ros.set_xticks([])
    ax_ros.set_ylabel("ROS_Cfam_1.0\n20 individuals", fontsize=6.5)
    ax_ros.set_yticks([])
    for sp in ["top", "right", "bottom"]:
        ax_ros.spines[sp].set_visible(False)

    # ---- TTAGGA block ----
    ax_t2t = fig.add_subplot(inner_gs[1])
    ax_t2t.imshow(t2t_mat, aspect="auto", cmap=_cmap(COLOR_T2T),
                  vmin=0, vmax=vmax, extent=extent, origin="lower", interpolation="nearest")
    ax_t2t.axvspan(CENTROMERE_START_MB, CENTROMERE_END_MB, color="#3366AA",
                   alpha=0.12, zorder=3)
    ax_t2t.set_xlim(0, CHRY_LENGTH_MB)
    ax_t2t.set_xticks([])
    ax_t2t.set_ylabel("Jindo1-G-TTAGGA Hap2\n20 individuals", fontsize=6.5)
    ax_t2t.set_yticks([])
    for sp in ["top", "right", "bottom"]:
        ax_t2t.spines[sp].set_visible(False)

    # ---- position axis ----
    ax_pos = fig.add_subplot(inner_gs[2])
    ax_pos.set_xlim(0, CHRY_LENGTH_MB)
    ax_pos.set_yticks([])
    ax_pos.set_xlabel("Position on chrY (Mb)", fontsize=7)
    ax_pos.tick_params(axis="x", length=2.5, width=0.5, labelsize=6)
    for sp in ["left", "top", "right"]:
        ax_pos.spines[sp].set_visible(False)
    ax_pos.axvspan(CENTROMERE_START_MB, CENTROMERE_END_MB, color="#AACCEE", alpha=0.4)
    ax_pos.text((CENTROMERE_START_MB + CENTROMERE_END_MB) / 2, 0.5, "centromere",
                ha="center", va="center", fontsize=5, color="#3366AA", style="italic")

    return [ax_ros, ax_t2t, ax_pos]


def main():
    fig_w = FIG_W_DOUBLE * 0.55
    fig_h = 3.2
    fig = plt.figure(figsize=(fig_w, fig_h))
    outer_gs = gridspec.GridSpec(1, 1, left=0.20, right=0.97, top=0.95, bottom=0.13, figure=fig)
    make_panel_b(fig, outer_gs[0, 0])
    out_pdf = OUT_DIR / "fig5b_chrY_population_v8.pdf"
    out_png = OUT_DIR / "fig5b_chrY_population_v8.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, bbox_inches="tight", dpi=600)
    print(f"\nSaved: {out_pdf}")
    print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()
