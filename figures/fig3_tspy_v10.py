#!/usr/bin/env python
"""
Fig 3 (TSPY ampliconic array panel) v7 — dog vs human comparison.

Changes from v6:
  - Partial (truncated) TSPY copies are now marked with a small downward
    arrowhead above the copy, because at ~0.8 kb the block is too thin for a
    colour difference alone to read. The lighter shade is kept as a secondary
    cue. In the dog array, 2 of 46 copies are ~821 bp (vs ~3.2 kb) and sit at
    array-segment starts. Human copies are all ~2.8 kb (no partials).
  - Caption should state: arrowheads mark partial (truncated, ~0.8 kb) TSPY
    copies.

Units (carried from v6): blocks are TSPY CODING regions (dog mean ~3.1 kb,
human ~2.8 kb), same definition for both species. The total-length difference
(493 vs 896 kb) comes from repeat-unit SPACING (dog ~10.7 kb, human ~20.8 kb).

  Dog (Jindo)      : 46 TSPY coding copies (2 partial), 493 kb array span
  Human (CHM13v2.0): 44 TSPY coding copies in the array, 896 kb array span
  Orthologous TSPY proteins (47% aa identity, e=1.8e-67).

  NOTE (v8): v6/v7 used human_tspy_coding.tsv (45 copies), which includes
  TSPY2 at chrY:10,224,454-10,227,249 (- strand). TSPY2 lies ~385 kb distal
  to the array in IR3 and is NOT part of the tandem array (Rhie et al. 2023,
  T2T-Y). Including it inflated the human span to 1,284 kb and the apparent
  periodicity to 28.5 kb/copy. v8 uses human_tspy_array_only.tsv (44 copies,
  the array proper), giving 895.7 kb and 20.8 kb/copy, consistent with the
  ~20.2 kb composite repeat unit reported for T2T-Y. The dog array is
  measured the same way (array only), so the two species are comparable.
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
from matplotlib.patches import Rectangle
import numpy as np

from _config import FIG_W_DOUBLE

DOG_TSV   = ROOT / "Results/Manuscript_Figures/data/tspy/tspy_array.tsv"
HUMAN_TSV = ROOT / "Results/Manuscript_Figures/data/tspy/human/human_tspy_array_only.tsv"
OUT_DIR = ROOT / "Results/Manuscript_Figures/output/fig3"
OUT_DIR.mkdir(parents=True, exist_ok=True)

DOG_COLOR           = "#C44E52"
DOG_COLOR_PARTIAL   = "#E8A0A2"
HUMAN_COLOR         = "#4C72B0"
HUMAN_COLOR_PARTIAL = "#A9C0DD"
ARROW_COLOR         = "#333333"

PARTIAL_BP = 1500

BAND_H   = 0.34
TRACK_BG = "#EDEDED"


def _load_dog():
    rows = []
    for line in open(DOG_TSV):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            rows.append((int(p[1]), int(p[2]), p[3]))
    return rows


def _load_human():
    rows = []
    for line in open(HUMAN_TSV):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3:
            rows.append((int(p[0]), int(p[1]), p[2]))
    return rows


def _draw_array(ax, rows, y, color_full, color_partial, label):
    start0 = min(r[0] for r in rows)
    # Span is measured across full-length copies only, so that dog (44 full +
    # 2 partial) and human (44 full) are compared on the same basis.
    _full = [r for r in rows if (r[1] - r[0]) >= PARTIAL_BP]
    span_kb = (max(r[1] for r in _full) - min(r[0] for r in _full)) / 1000.0

    ax.add_patch(Rectangle((0, y - BAND_H/2), span_kb, BAND_H,
                 facecolor=TRACK_BG, edgecolor="none", zorder=1))

    n_partial = 0
    partial_x = []
    for s, e, strand in rows:
        x0 = (s - start0) / 1000.0
        w  = (e - s) / 1000.0
        is_partial = (e - s) < PARTIAL_BP
        if is_partial:
            n_partial += 1
            partial_x.append(x0 + w/2.0)
        if w <= 0:
            w = 0.3
        ax.add_patch(Rectangle((x0, y - BAND_H/2), w, BAND_H,
                     facecolor=(color_partial if is_partial else color_full),
                     edgecolor="black", linewidth=0.25, zorder=3))

    for px in partial_x:
        ax.plot(px, y + BAND_H/2 + 0.18, marker="v", markersize=4,
                color=ARROW_COLOR, markeredgecolor="none", zorder=5,
                clip_on=False)

    n = len(rows)
    ax.text(-30, y, label, ha="right", va="center", fontsize=8, fontweight="bold")
    n_full = n - n_partial
    if n_partial > 0:
        annot = f"{n_full} full-length + {n_partial} partial \u00b7 {span_kb:.0f} kb"
    else:
        annot = f"{n_full} full-length \u00b7 {span_kb:.0f} kb"
    ax.text(span_kb + 20, y, annot, ha="left", va="center", fontsize=7)
    return span_kb


def make_panel(fig, gs_subplot, show_xlabel=True):
    dog = _load_dog()
    human = _load_human()

    ax = fig.add_subplot(gs_subplot)
    s_h = _draw_array(ax, human, 1.0, HUMAN_COLOR, HUMAN_COLOR_PARTIAL, "Human\n(CHM13v2.0)")
    s_d = _draw_array(ax, dog,   0.0, DOG_COLOR,   DOG_COLOR_PARTIAL,   "Dog\n(Jindo)")

    ax.set_xlim(-260, max(s_h, s_d) + 230)
    ax.set_ylim(-0.6, 1.7)
    ax.set_yticks([])
    if show_xlabel:
        ax.set_xlabel("Array position (kb, each array zeroed at its start)",
                      fontsize=8, fontweight="bold")
    ax.tick_params(labelsize=7)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    return ax


def main():
    fig = plt.figure(figsize=(FIG_W_DOUBLE, 2.4))
    gs = gridspec.GridSpec(1, 1, figure=fig, left=0.12, right=0.93, top=0.92, bottom=0.22)
    make_panel(fig, gs[0])
    plt.savefig(OUT_DIR / "fig3_tspy_v10.pdf", bbox_inches="tight")
    plt.savefig(OUT_DIR / "fig3_tspy_v10.png", dpi=300, bbox_inches="tight")
    print("Saved tspy v10 (partial copies marked with arrowheads)")


if __name__ == "__main__":
    main()
