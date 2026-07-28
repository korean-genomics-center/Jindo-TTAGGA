#!/usr/bin/env python
"""
Figure 3 combined v10 — the complete canine chrY. THREE panels.
v10 change from v9: panel c uses fig3_tspy_v7 (block-band TSPY with partial
copies marked by arrowheads) instead of fig3_tspy_v4 (squashed arrows).

Panels:
  a  chrY length across published male canid assemblies + Jindo (fig3d_v16;
     chrY-absent assemblies excluded, unified breed labels, two-colour legend)
  b  chrY multi-track map (genes + 5 repeat classes) + CpG 5mC methylation on
     the shared 0-21.3 Mb x-axis (fig3a_v6; structure panel b + methylation
     panel c merged)
  c  TSPY ampliconic coding array, dog vs human CHM13v2.0, block-band with
     partial copies arrow-marked (fig3_tspy_v7)
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
from fig3d_v18_canid_landscape   import make_panel_d
from fig3a_v7_chrY_map_methyl    import make_panel_a
from fig3_tspy_v10                import make_panel as make_panel_tspy
OUT_DIR = ROOT / "Results/Manuscript_Figures/output/fig3"
OUT_DIR.mkdir(parents=True, exist_ok=True)
def main():
    fig = plt.figure(figsize=(FIG_W_DOUBLE, 12.5))
    outer = gridspec.GridSpec(
        3, 1, figure=fig,
        left=0.10, right=0.95, top=0.965, bottom=0.045,
        hspace=0.30,
        height_ratios=[2.5, 4.6, 2.0],
    )
    make_panel_d(fig, outer[0, 0], show_title=False, show_legend=True)   # a
    make_panel_a(fig, outer[1, 0], show_title=False)                     # b
    make_panel_tspy(fig, outer[2, 0], show_xlabel=True)                  # c
    fig.text(0.025, 0.960, "a", fontsize=14, fontweight="bold", va="bottom")
    fig.text(0.025, 0.660, "b", fontsize=14, fontweight="bold", va="bottom")
    fig.text(0.025, 0.235, "c", fontsize=14, fontweight="bold", va="bottom")
    out_pdf = OUT_DIR / "fig3_combined_v14.pdf"
    out_png = OUT_DIR / "fig3_combined_v14.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, bbox_inches="tight", dpi=300)
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")
if __name__ == "__main__":
    main()
