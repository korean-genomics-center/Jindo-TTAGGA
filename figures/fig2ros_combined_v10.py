#!/usr/bin/env python
"""
Figure 2 combined v9 — ROS/reference comparison. 3 full-width rows (a/b/c).
v9 change from v8: panel a now uses fig2ros_a_v4 (chr9 inversion coordinate
fix — HAP1/HAP2 inversion coordinates were swapped in v3/v8; corrected
against actual SyRI INVAL records). Panels b/c unchanged.
Figure 2 = a chr9 inversion / b per-chromosome length / c continuity
benchmark vs ROS_Cfam_1.0.
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
from fig2ros_a_v4_chr9_inversion import make_panel_a
from fig2ros_b_v6_length         import make_panel_b
from fig2ros_d_v8_benchmark      import make_panel_c as make_panel_benchmark
OUT_DIR = ROOT / "Results/Manuscript_Figures/output/fig2"
OUT_DIR.mkdir(parents=True, exist_ok=True)
def main():
    fig = plt.figure(figsize=(FIG_W_DOUBLE*1.05, 8.5))
    outer = gridspec.GridSpec(
        3, 1, figure=fig,
        left=0.07, right=0.97, top=0.955, bottom=0.07,
        hspace=0.42,
        height_ratios=[2.4, 2.4, 2.2],
    )
    make_panel_a(fig, outer[0, 0])              # a: chr9 inversion (v4, fixed)
    make_panel_b(fig, outer[1, 0])              # b: per-chromosome length
    make_panel_benchmark(fig, outer[2, 0], show_footnote=False)  # c: continuity
    fig.text(0.012, 0.945, "a", fontsize=13, fontweight="bold", va="top")
    fig.text(0.012, 0.620, "b", fontsize=13, fontweight="bold", va="top")
    fig.text(0.012, 0.300, "c", fontsize=13, fontweight="bold", va="top")
    out_pdf = OUT_DIR / "fig2ros_combined_v10.pdf"
    out_png = OUT_DIR / "fig2ros_combined_v10.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, bbox_inches="tight", dpi=300)
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")
if __name__ == "__main__":
    main()
