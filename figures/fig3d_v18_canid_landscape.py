#!/usr/bin/env python
"""
Figure 3D v17: chrY landscape across male canid assemblies + Jindo Hap2.

Changes from v16:
    - ALL chrY values re-measured directly from assembly FASTA (2026-07-21).
      Definition applied uniformly: gap = run of >= 10 consecutive N.
    - ROS_Cfam_1.0 : unplaced 2.37 -> 2.79 Mb, gaps 7 -> 15, N 687,168 -> 456,949
    - Yella v2     : gaps 14 -> 16, N 1,400 -> 1,596
    - Bernese OD   : gaps 18 -> 19, N 1,800 -> 1,898
    - Cairn CA611, Basenji, Jindo Hap2: measured values match v16; unchanged.
    - Dead "is_ros" key removed (the bar-colour branch was dropped in v16).
    - Legend states the chrY-assignment rule and the gap definition.

Measured values (chr-scale / unplaced / N / gaps, gap = N-run >= 10 bp):
    Jindo Hap2  21,255,906 /         0 /       0 /  0
    ROS_Cfam1.0  3,937,623 / 2,790,983 / 456,949 / 15
    CA611        3,540,084 / 2,139,221 /   1,300 / 13
    Yella v2             0 / 5,366,740 /   1,596 / 16
    OD           3,314,329 / 1,960,291 /   1,898 / 19
    Basenji      3,626,350 /         0 /     600 /  6

Component sums:
    ROS unplaced   1,569,522 + 1,221,461 = 2,790,983
    ROS N          146 + 189,047 + 267,756 = 456,949   gaps 5+6+4 = 15
    CA611 unplaced 1,521,396 + 617,825 = 2,139,221
    CA611 N        1,000 + 0 + 300 = 1,300             gaps 10+0+3 = 13
    OD unplaced    1,366,912 + 593,379 = 1,960,291
    OD N           1,000 + 598 + 300 = 1,898           gaps 10+6+3 = 19
    Yella unplaced 3,348,343 + 1,426,613 + 591,784 = 5,366,740
    Yella N        1,096 + 200 + 300 = 1,596           gaps 11+2+3 = 16

Why v11-v16 were wrong:
    v11 ("NCBI corrected") changed the ROS row to match manuscript V15 text
    rather than re-measuring. Scaffold_34 (NW_024010443.1) was recorded as
    1.15 Mb instead of its actual 1,569,522 bp, which shrank unplaced to
    2.37 Mb while simultaneously inflating N to 687,168. Yella and OD gap/N
    were carried from v9 without measurement. v9/v10 ROS N (456,877) was in
    fact close to the true value (456,949); v11 overwrote a correct number
    with an incorrect one.

chrY assignment rule (differs by assembly, stated in the legend):
    ROS_Cfam_1.0 : chromosome-scale chrY + the two NCBI chrY-unlocalized
                   scaffolds (NW_024010443.1, NW_024010444.1).
    Basenji      : chromosome-scale chrY only (no chrY-assigned scaffolds).
    CA611 / OD   : chromosome-scale chrY + the two largest unplaced scaffolds
                   absent from chrX (JARDRD/JARDRE ...041, ...042).
    Yella v2     : no NCBI chrY assignment; the three largest unplaced contigs
                   absent from chrX (DAOUOP010000040/041/042). Marked "*".

Gap-size caveat:
    The ONT assemblies (CA611, OD, Yella v2) and the Basenji assembly use
    fixed 100-bp N placeholders, so their "N total" is notation, not missing
    sequence. ROS_Cfam_1.0 gaps retain estimated physical sizes (mean ~30 kb),
    so its 457 kb of N represents genuinely unresolved sequence.

Usage:
    Standalone: python fig3d_v18_canid_landscape.py
    Imported:   from fig3d_v18_canid_landscape import make_panel_d
                make_panel_d(fig, gs_subplot)
"""

import sys
from pathlib import Path

SCRIPTS_DIR = Path("${JINDO_ROOT}/Results/Manuscript_Figures/scripts")
sys.path.insert(0, str(SCRIPTS_DIR))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

from _config import COLORS, FIG_W_DOUBLE

ROOT = Path("${JINDO_ROOT}")
OUT_DIR = ROOT / "Results/Manuscript_Figures/output/fig3"
OUT_DIR.mkdir(parents=True, exist_ok=True)

EXCLUDE_ABSENT = True

COLOR_JINDO   = COLORS.get("jindo_hap2", "#E69F00")
COLOR_OTHER   = "#5D6D7E"
COLOR_ABSENT  = "#D7DBDD"
HATCH_UNPL    = "////"

# ============================================================
# Data: measured directly from assembly FASTA (2026-07-21)
# gap = run of >= 10 consecutive N; see header for per-assembly rule.
# DO NOT edit these numbers to match manuscript text. If the manuscript and
# this table disagree, re-measure from FASTA and fix the manuscript.
# ============================================================
DATA_RAW = [
    {"name": "Jindo Hap2 (this study)",     "breed": "Jindo",                "asm": "Hap2, this study", "tech2": "HiFi + ONT-UL, Trio",   "main_mb": 21.26, "unplaced_mb": 0.00, "gaps": 0,    "total_n": 0,        "year": 2026, "tech": "HiFi + ONT-UL + Trio",  "absent": False, "highlight": True},
    {"name": "ROS_Cfam_1.0 (Labrador)",     "breed": "Labrador",             "asm": "ROS_Cfam_1.0",     "tech2": "PacBio Sequel + HiC",   "main_mb": 3.94,  "unplaced_mb": 2.79, "gaps": 15,   "total_n": 456949,   "year": 2020, "tech": "PacBio Sequel + HiC",   "absent": False},
    {"name": "Cairn Terrier (CA611)",       "breed": "Cairn Terrier",        "asm": "CA611",            "tech2": "ONT PromethION",        "main_mb": 3.54,  "unplaced_mb": 2.14, "gaps": 13,   "total_n": 1300,     "year": 2023, "tech": "ONT PromethION",        "absent": False},
    {"name": "Yella v2 (Labrador)",         "breed": "Labrador",             "asm": "Yella v2",         "tech2": "ONT GridION",           "main_mb": 0.00,  "unplaced_mb": 5.37, "gaps": 16,   "total_n": 1596,     "year": 2023, "tech": "ONT GridION",           "absent": False, "tentative": True},
    {"name": "Bernese Mtn Dog (OD)",        "breed": "Bernese Mountain Dog", "asm": "OD",               "tech2": "ONT PromethION",        "main_mb": 3.31,  "unplaced_mb": 1.96, "gaps": 19,   "total_n": 1898,     "year": 2023, "tech": "ONT PromethION",        "absent": False},
    {"name": "Basenji (breed-1.1, Wags)",   "breed": "Basenji",              "asm": "breed-1.1",        "tech2": "PacBio Sequel",         "main_mb": 3.63,  "unplaced_mb": 0.00, "gaps": 6,    "total_n": 600,      "year": 2019, "tech": "PacBio Sequel",         "absent": False},
    {"name": "Yella alt (Labrador)",        "breed": "Labrador",             "asm": "Yella alt",        "tech2": "ONT GridION",           "main_mb": 0.00,  "unplaced_mb": 0.00, "gaps": None, "total_n": None,     "year": 2020, "tech": "ONT GridION",           "absent": True},
    {"name": "Yella principal (Labrador)",  "breed": "Labrador",             "asm": "Yella principal",  "tech2": "ONT GridION",           "main_mb": 0.00,  "unplaced_mb": 0.00, "gaps": None, "total_n": None,     "year": 2020, "tech": "ONT GridION",           "absent": True},
    {"name": "Whippet (TBG_BS_Lino)",       "breed": "Whippet",              "asm": "TBG_BS_Lino",      "tech2": "PacBio HiFi Revio",     "main_mb": 0.00,  "unplaced_mb": 0.00, "gaps": None, "total_n": None,     "year": 2024, "tech": "PacBio HiFi Revio",     "absent": True},
]


def _sort_key(row):
    if row["absent"]:
        return (1, row["name"])
    total = row["main_mb"] + row["unplaced_mb"]
    return (0, -total)


_DATA_SRC = [r for r in DATA_RAW if not (EXCLUDE_ABSENT and r["absent"])]
DATA = sorted(_DATA_SRC, key=_sort_key)


def _fmt_gap_annotation(gaps, total_n):
    if gaps is None or total_n is None:
        return ""
    if gaps == 0:
        return "0 gaps"
    if total_n >= 10000:
        return f"{gaps} gaps, {total_n/1000:.0f} kb total N"
    elif total_n >= 1000:
        return f"{gaps} gaps, {total_n/1000:.1f} kb total N"
    else:
        return f"{gaps} gaps, {total_n} bp total N"


def _row_label(row):
    line1 = f"{row['breed']} ({row['asm']})"
    line2 = f"({row['year']}, {row['tech2']})"
    return f"{line1}\n{line2}"


def make_panel_d(fig, gs_subplot=None, ax=None, show_title=True, show_legend=True):
    if ax is None:
        if gs_subplot is None:
            ax = fig.add_subplot(111)
        else:
            ax = fig.add_subplot(gs_subplot)

    n = len(DATA)
    y_pos = np.arange(n)
    any_absent = False

    for i, row in enumerate(DATA):
        if row["absent"]:
            any_absent = True
            ax.barh(i, 0.4, color=COLOR_ABSENT, edgecolor="#909497",
                    linewidth=0.5, height=0.6)
            ax.text(0.55, i, "chrY absent in assembly",
                    va="center", ha="left",
                    fontsize=6, style="italic", color="#7B7D7D")
        else:
            bar_color = COLOR_JINDO if row.get("highlight") else COLOR_OTHER

            if row["main_mb"] > 0:
                ax.barh(i, row["main_mb"], color=bar_color,
                        edgecolor="black", linewidth=0.4, height=0.6)

            if row["unplaced_mb"] > 0:
                ax.barh(i, row["unplaced_mb"], left=row["main_mb"],
                        color=bar_color, edgecolor="black",
                        linewidth=0.4, height=0.6,
                        hatch=HATCH_UNPL, alpha=0.55)

            total = row["main_mb"] + row["unplaced_mb"]
            gap_txt = _fmt_gap_annotation(row["gaps"], row["total_n"])
            annot = f"{total:.2f} Mb  ({gap_txt})" if gap_txt else f"{total:.2f} Mb"
            if row.get("tentative"):
                annot += "*"
            ax.text(total + 0.4, i, annot,
                    va="center", ha="left", fontsize=6, color="#2C3E50")

    y_labels = [_row_label(row) for row in DATA]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("chrY length (Mb)", fontsize=7)
    ax.set_xlim(0, 24)
    ax.tick_params(axis='x', labelsize=6)
    ax.tick_params(axis='y', labelsize=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    legend_elements = [
        Patch(facecolor=COLOR_JINDO, edgecolor='black', linewidth=0.4, label='Jindo Hap2 (this study)'),
        Patch(facecolor=COLOR_OTHER, edgecolor='black', linewidth=0.4, label='Published male dog assemblies'),
        Patch(facecolor='white', edgecolor='black', linewidth=0.4, hatch=HATCH_UNPL, label='Unplaced chrY-assigned contigs'),
    ]
    if any_absent:
        legend_elements.append(
            Patch(facecolor=COLOR_ABSENT, edgecolor='#909497', linewidth=0.5, label='chrY absent in assembly')
        )
    ax.legend(handles=legend_elements, loc='lower right', fontsize=5.5,
              frameon=True, framealpha=0.95, edgecolor='#BDC3C7')

    if any(r.get("tentative") for r in DATA):
        ax.text(0.0, -0.14,
                "* Yella v2 carries no NCBI chrY assignment; the three largest "
                "unplaced contigs absent from chrX were tentatively assigned.",
                transform=ax.transAxes, fontsize=5, style="italic",
                color="#5D6D7E", va="top")
    return ax


if __name__ == "__main__":
    fig_w = FIG_W_DOUBLE * 0.75
    fig_h = 3.2
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    make_panel_d(fig, ax=ax)
    plt.tight_layout()
    out_pdf = OUT_DIR / "fig3d_v18.pdf"
    out_png = OUT_DIR / "fig3d_v18.png"
    plt.savefig(out_pdf, dpi=600, bbox_inches='tight')
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")
    print()
    for row in DATA:
        total = row["main_mb"] + row["unplaced_mb"]
        t = " (tentative chrY assignment)" if row.get("tentative") else ""
        print(f"  {row['name']:34} : main={row['main_mb']:5.2f} + unplaced={row['unplaced_mb']:5.2f} "
              f"= {total:5.2f} Mb, {row['gaps']} gaps, {row['total_n']} N bp{t}")
