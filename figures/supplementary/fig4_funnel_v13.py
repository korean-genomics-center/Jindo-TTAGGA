#!/usr/bin/env python
"""
Fig4 funnel v10 - reference-absent gene pipeline (clean vertical spine).
Supersedes v9: removed the left/right OrthoFinder family/novel boxes (2,317/961);
added an explicit artefact-filtering box (-115 artefacts) on the spine.
Flow:
  36,006 -> missing>=1 ref 4,046 -> absent from all 4 = 3,278
  -> OrthoFinder classification -> [artefact filtering -115] -> 3,163 clean
  -> family 2,203 + orthogroup-unassigned (label renamed from complete-novel, v13)
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
from matplotlib.patches import FancyBboxPatch

try:
    from _config import FIG_W_DOUBLE
except Exception:
    FIG_W_DOUBLE = 7.2

OUT_DIR = ROOT / "Results/Manuscript_Figures/output/fig4"
OUT_DIR.mkdir(parents=True, exist_ok=True)

C_SERIAL="#9CC3E0"; C_SERIAL_E="#0072B2"
C_CORE="#7FBFB3";  C_CORE_E="#1B7A6B"
C_FILTER="#F0E2C8"; C_FILTER_E="#B8923E"
C_FINAL="#1B7A6B"

# Vertical spine (top -> bottom)
steps = [
    ("BRAKER3 gene predictions",          "36,006", "serial"),
    ("Missing from \u22651 dog reference", "4,046",  "serial"),
    ("Absent from all 4 dog references",  "3,278",  "core"),
    ("OrthoFinder classification",        "",       "corehdr"),
    ("Artefact filtering (\u2212115)",    "",       "filter"),
    ("After artefact removal",            "3,163",  "serial"),
    ("Present on re-search\n(50 targets; \u221245)", "",  "filter"),
    ("Not evaluable\n(<50 aa or low complexity; \u2212747)", "", "filter"),
    ("Reference-absent (refined)",        "2,371",  "final"),
]

def make_panel_funnel(fig, gs_subplot):
    ax = fig.add_subplot(gs_subplot)
    n = len(steps)
    box_w, box_h = 0.54, 0.60
    x0 = 0.42
    gap = 1.02

    def box(xc, yc, w, h, label, num, fc, ec, tc, bold=False, num_fs=8.0, lab_fs=5.8):
        ax.add_patch(FancyBboxPatch((xc-w/2, yc-h/2), w, h,
                     boxstyle="round,pad=0.02,rounding_size=0.06",
                     facecolor=fc, edgecolor=ec, linewidth=0.9))
        if num:
            ax.text(xc, yc+0.14, label, ha="center", va="center", fontsize=lab_fs,
                    color=tc, fontweight=("bold" if bold else "normal"))
            ax.text(xc, yc-0.16, num, ha="center", va="center", fontsize=num_fs,
                    color=tc, fontweight="bold")
        else:
            ax.text(xc, yc, label, ha="center", va="center", fontsize=lab_fs,
                    color=tc, fontweight=("bold" if bold else "normal"))

    yc_of = {}
    for i,(label,num,grp) in enumerate(steps):
        yc = (n-1-i)*gap
        yc_of[i] = yc
        if grp=="serial":    fc,ec,tc = C_SERIAL,  C_SERIAL_E, "#0a2740"
        elif grp=="core":    fc,ec,tc = C_CORE,    C_CORE_E,   "#0a2d28"
        elif grp=="corehdr": fc,ec,tc = "#D8ECE7", C_CORE_E,   "#0a2d28"
        elif grp=="filter":  fc,ec,tc = C_FILTER,  C_FILTER_E, "#5a4416"
        else:                fc,ec,tc = C_FINAL,   C_CORE_E,   "white"
        # filter box a bit shorter
        h = 0.50 if grp in ("corehdr","filter") else box_h
        box(x0, yc, box_w, h, label, num, fc, ec, tc, bold=(grp=="final"))
        if i < n-1:
            ax.annotate("", xy=(x0, yc-h/2-(gap-h)+0.02),
                        xytext=(x0, yc-h/2-0.01),
                        arrowprops=dict(arrowstyle="-|>", color="#888888", lw=1.0))

    # Final split: clean 3,163 -> family 2,203 + complete-novel 960
    y_clean = yc_of[len(steps)-1]
    y_leaf = y_clean - gap
    lfx = x0 - 0.34
    lnx = x0 + 0.34
    lw, lh = 0.46, 0.56
    ax.annotate("", xy=(lfx, y_leaf+lh/2), xytext=(x0-0.10, y_clean-box_h/2),
                arrowprops=dict(arrowstyle="-|>", color="#888888", lw=1.0))
    ax.annotate("", xy=(lnx, y_leaf+lh/2), xytext=(x0+0.10, y_clean-box_h/2),
                arrowprops=dict(arrowstyle="-|>", color="#888888", lw=1.0))
    ax.add_patch(FancyBboxPatch((lfx-lw/2, y_leaf-lh/2), lw, lh,
                 boxstyle="round,pad=0.02,rounding_size=0.06",
                 facecolor=C_CORE, edgecolor=C_CORE_E, linewidth=0.8))
    ax.text(lfx, y_leaf+0.13, "Family members", ha="center", va="center", fontsize=5.0, color="#0a2d28")
    ax.text(lfx, y_leaf-0.14, "1,708", ha="center", va="center", fontsize=7.2, color="#0a2d28", fontweight="bold")
    ax.add_patch(FancyBboxPatch((lnx-lw/2, y_leaf-lh/2), lw, lh,
                 boxstyle="round,pad=0.02,rounding_size=0.06",
                 facecolor=C_FINAL, edgecolor=C_CORE_E, linewidth=0.8))
    ax.text(lnx, y_leaf+0.13, "Orthogroup-\nunassigned", ha="center", va="center", fontsize=4.4, color="white")
    ax.text(lnx, y_leaf-0.14, "663", ha="center", va="center", fontsize=7.2, color="white", fontweight="bold")

    ax.set_xlim(lfx-lw/2-0.12, lnx+lw/2+0.12)
    ax.set_ylim(y_leaf-lh/2-0.4, (n-1)*gap+0.7)
    ax.axis("off")
    return ax

def main():
    fig = plt.figure(figsize=(FIG_W_DOUBLE*0.52, 7.0))
    gs = gridspec.GridSpec(1,1,figure=fig)
    make_panel_funnel(fig, gs[0])
    plt.savefig(OUT_DIR/"fig4_funnel_v13.pdf", bbox_inches="tight")
    plt.savefig(OUT_DIR/"fig4_funnel_v13.png", bbox_inches="tight", dpi=600)
    print("Saved v10 (clean vertical funnel)")

if __name__=="__main__":
    main()
