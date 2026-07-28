#!/usr/bin/env python
"""
Figure 1 combined v17 — 2 panels (a / b).

Changes from v16:
  - Panel c (completeness metrics: BUSCO + Merqury QV + GCI) REMOVED from Fig 1
    and moved to Supplementary Figure 7 (Sfig7_completeness_v1.py).
  - Fig 1 now carries a single message: "who, and how complete the assembly is
    structurally":
       a  sample provenance (Korea map + Baeksan photo + trio pedigree)
       b  chromosome landscape (Hap1/Hap2 stacked, gene-density shaded,
          centromere + telomere)  -- the HERO panel
  - Figure height reduced 12.5 -> 10.6 in and height_ratios rebalanced so that
    panel b keeps the same physical height as in v16 (no trailing whitespace).
  - Panel letter for b moved up accordingly.
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

from _config import FIG_W_DOUBLE

from fig1a_v13_final     import make_panel_a
from fig1b_landscape_v15 import make_panel_b as make_panel_landscape

OUT_DIR = ROOT / "Results/Manuscript_Figures/output/fig1"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    fig_w = FIG_W_DOUBLE

    # v16 geometry, for reference:
    #   fig_h = 12.5, top=0.985, bottom=0.04, hspace=0.22
    #   height_ratios = [1.4, 7.0, 1.6]  (a / b / c)
    fig_h = 10.6
    fig = plt.figure(figsize=(fig_w, fig_h))

    # 2 rows: a (sample, compact) / b (landscape, HERO)
    outer = gridspec.GridSpec(
        2, 1, figure=fig,
        left=0.06, right=0.97, top=0.985, bottom=0.045,
        hspace=0.20,
        height_ratios=[1.4, 7.0],
    )

    make_panel_a(fig, outer[0, 0], show_seq_strategy=False, show_panel_letter=False)
    make_panel_landscape(fig, outer[1, 0])

    # panel letters (figure coordinates; b raised because c is gone)
    fig.text(0.012, 0.978, "a", fontsize=12, fontweight="bold", va="top")
    fig.text(0.012, 0.815, "b", fontsize=12, fontweight="bold", va="top")

    out_pdf = OUT_DIR / "fig1_combined_v17.pdf"
    out_png = OUT_DIR / "fig1_combined_v17.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, bbox_inches="tight", dpi=600)
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()
