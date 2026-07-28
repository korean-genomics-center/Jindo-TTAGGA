#!/usr/bin/env python
"""
Fig 2 Panel B v6: Per-chromosome length comparison vs ROS.

v6 changes from v5:
  - Refactored into make_panel_b(fig, gs_subplot) for the combined figure.
  - All paths moved to backup base (${JINDO_ROOT}).
  - Order/colors unchanged: Maternal (red) / Paternal (blue) / ROS (gray),
    already matching the Jeju co-author feedback.
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
import numpy as np

from _config import FIG_W_DOUBLE

HAP1_FAI = ROOT / "Results/Manuscript_Figures/data/assembly_fai/JindoHap1.fai"
HAP2_FAI = ROOT / "Results/Manuscript_Figures/data/assembly_fai/JindoHap2.fai"
ROS_FAI  = ROOT / "Results/Manuscript_Figures/data/reference_fai/ROS_Cfam1.0.fai"
OUT_DIR  = ROOT / "Results/Manuscript_Figures/output/fig2"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def load_fai(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            p = line.strip().split("\t")
            d[p[0]] = int(p[1])
    return d


def _load():
    hap1 = load_fai(HAP1_FAI)
    hap2 = load_fai(HAP2_FAI)
    ros_raw = load_fai(ROS_FAI)
    canonical = [f"chr{i}" for i in range(1, 39)] + ["chrX", "chrY"]
    hap1 = {c: hap1[c] for c in canonical if c in hap1}
    hap2 = {c: hap2[c] for c in canonical if c in hap2}
    ros = {}
    for i in range(38):
        cm_id = f"CM025{100+i:03d}.1"
        if cm_id in ros_raw:
            ros[f"chr{i+1}"] = ros_raw[cm_id]
    if "CM025138.1" in ros_raw:
        ros["chrX"] = ros_raw["CM025138.1"]
    if "CM025139.1" in ros_raw:
        ros["chrY"] = ros_raw["CM025139.1"]
    return hap1, hap2, ros


def make_panel_b(fig, gs_subplot):
    hap1, hap2, ros = _load()
    chr_order = [f"chr{i}" for i in range(1, 39)] + ["chrX", "chrY"]
    ros_lens  = [ros.get(c, 0)  / 1e6 for c in chr_order]
    hap1_lens = [hap1.get(c, 0) / 1e6 for c in chr_order]
    hap2_lens = [hap2.get(c, 0) / 1e6 for c in chr_order]

    ax = fig.add_subplot(gs_subplot)
    x = np.arange(len(chr_order))
    width = 0.22
    ax.bar(x - width, hap1_lens, width, label="Jindo1-G-TTAGGA-M (Maternal)",
           color="#D62728", edgecolor="black", linewidth=0.3)
    ax.bar(x, hap2_lens, width, label="Jindo1-G-TTAGGA-P (Paternal)",
           color="#1F77B4", edgecolor="black", linewidth=0.3)
    ax.bar(x + width, ros_lens, width, label="ROS_Cfam_1.0",
           color="#888888", edgecolor="black", linewidth=0.3)
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("chr", "") for c in chr_order], rotation=0, fontsize=6)
    ax.set_xlabel("Chromosome", fontsize=8, fontweight="bold")
    ax.set_ylabel("Length (Mb)", fontsize=8, fontweight="bold")
    ax.legend(fontsize=7, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.12), frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", linestyle=":", alpha=0.3, linewidth=0.3)
    ax.set_xlim(-0.5, len(chr_order) - 0.5)
    ymax = max(max(ros_lens), max(hap1_lens), max(hap2_lens))
    ax.set_ylim(0, ymax * 1.10)
    return ax


def main():
    fig = plt.figure(figsize=(FIG_W_DOUBLE, 4.0))
    gs = gridspec.GridSpec(1, 1, figure=fig)
    make_panel_b(fig, gs[0])
    plt.savefig(OUT_DIR / "fig2ros_b_v6.pdf", dpi=600, bbox_inches="tight")
    plt.savefig(OUT_DIR / "fig2ros_b_v6.png", dpi=300, bbox_inches="tight")
    print("Saved v6")


if __name__ == "__main__":
    main()
