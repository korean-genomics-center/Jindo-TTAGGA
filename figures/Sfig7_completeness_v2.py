#!/usr/bin/env python
"""
Supplementary Figure 7 — assembly completeness and consensus accuracy.

Content moved out of main Figure 1 panel c (v16) without changing any value:
    BUSCO completeness (carnivora_odb12, n = 13,727) for TTAGGA-M, TTAGGA-P
    and three published canine references; Merqury k-mer consensus quality
    (QV, k = 21); GCI score (HiFi + ONT UL).

Implementation note
-------------------
This script deliberately re-uses make_panel_b() from fig1b_v4_busco_refactored,
i.e. exactly the same function that drew panel c of Figure 1 v16. Numbers are
therefore guaranteed to be identical to the previously circulated main figure;
nothing is re-entered by hand.

Do NOT use the older scripts/suppfig1_busco_qv.py: its values do not match
Figure 1 v16 and it omits the GCI sub-panel entirely.

NC style: no figure-level title, lower-case panel letters, no interpretive
annotation.
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

from fig1b_v4_busco_refactored import make_panel_b as make_panel_completeness

OUT_DIR = ROOT / "Results/Manuscript_Figures/output/supp"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    fig_w = FIG_W_DOUBLE
    fig_h = 3.4
    fig = plt.figure(figsize=(fig_w, fig_h))

    outer = gridspec.GridSpec(
        1, 1, figure=fig,
        left=0.06, right=0.97, top=0.90, bottom=0.14,
    )

    make_panel_completeness(fig, outer[0, 0])

    out_pdf = OUT_DIR / "S_Fig7_completeness_v2.pdf"
    out_png = OUT_DIR / "S_Fig7_completeness_v2.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, bbox_inches="tight", dpi=600)
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()
