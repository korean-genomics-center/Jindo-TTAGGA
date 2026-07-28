#!/usr/bin/env python3
"""
Figure 1A v13: Korea map + Baeksan photo + Trio pedigree.

Layout (3 sub-panels horizontally):
    LEFT-LEFT  : Korean peninsula (cartopy Natural Earth)
    LEFT-RIGHT : Baeksan photo + caption
    RIGHT      : Trio pedigree (Sire square -- Dam circle -> Baeksan square)

Changes from v9:
  - Trio sub-panel coordinates compacted (12x9.5 -> 10x8) for tighter scaling
  - All Trio fonts reduced (matching combined figure aspect ratio):
      * "Trio family" title: 11 -> 9.5
      * "Sire/Dam (paternal/maternal)": 9.5 -> 7.5
      * "Baekho/Beodeul": 9 -> 7
      * "Baeksan (offspring, F1)": 9.5 -> 7.5
      * "assembled individual": 8 -> 6.5
      * Symbol legend: 8 -> 6.5
  - Symbol legend repositioned to bottom (y=-0.4) with more clearance
  - width_ratios [0.9, 1.0, 1.5] -> [0.85, 1.0, 1.8] (more room for trio)
  - Vertical spacing between Sire/Dam labels and Baekho/Beodeul tightened
  - Reduces overlap with both standalone and combined figure usage

Usage: same as v9.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import cartopy.crs as ccrs
import cartopy.feature as cfeature

SCRIPTS_DIR = Path("${JINDO_ROOT}/Results/Manuscript_Figures/scripts")
sys.path.insert(0, str(SCRIPTS_DIR))


# ============================================================
# Constants
# ============================================================
ROOT = Path("${JINDO_ROOT}")
DATA_DIR = ROOT / "Results/Manuscript_Figures/data"
OUT_DIR  = ROOT / "Results/Manuscript_Figures/output/fig1/panels"

BAEKSAN_PHOTO = DATA_DIR / "photos" / "Baeksan.png"

KOREA_EXTENT = [124.5, 131.0, 33.0, 39.5]
JINDO_LAT, JINDO_LON = 34.48, 126.27

COLOR_LAND      = "#E8E8E8"
COLOR_GRAY_EDGE = "#7F7F7F"


# ============================================================
# Drawing primitives
# ============================================================
def _draw_korea_map(ax_map):
    ax_map.set_extent(KOREA_EXTENT, crs=ccrs.PlateCarree())
    ax_map.add_feature(cfeature.LAND,
                       facecolor=COLOR_LAND, edgecolor=COLOR_GRAY_EDGE, linewidth=0.6)
    ax_map.add_feature(cfeature.OCEAN, facecolor='white')
    ax_map.add_feature(cfeature.COASTLINE, edgecolor=COLOR_GRAY_EDGE, linewidth=0.6)
    ax_map.plot(JINDO_LON, JINDO_LAT, marker='*', color='black',
                markersize=11, markeredgewidth=0.5,
                transform=ccrs.PlateCarree())


def _draw_baeksan_photo(ax_photo):
    ax_photo.set_xlim(0, 1)
    ax_photo.set_ylim(0, 1)
    ax_photo.set_aspect('equal')
    ax_photo.axis('off')

    if BAEKSAN_PHOTO.exists():
        img = mpimg.imread(str(BAEKSAN_PHOTO))
        import numpy as np
        img = np.fliplr(img)
        ax_photo.imshow(img, extent=(0.05, 0.95, 0.30, 0.95),
                        aspect='auto', zorder=2)
    else:
        ax_photo.add_patch(mpatches.Rectangle(
            (0.05, 0.30), 0.90, 0.65,
            linewidth=1.5, edgecolor='black', facecolor='white'))
        ax_photo.text(0.5, 0.625, 'Baeksan\nphoto\n(missing)',
                      ha='center', va='center', fontsize=8.5,
                      color=COLOR_GRAY_EDGE, style='italic')

    # Caption (smaller fonts to match v10 trio)
    ax_photo.text(0.5, 0.20, 'Korean Jindo dog',
                  ha='center', va='top', fontsize=8.5, fontweight='bold')
    ax_photo.text(0.5, 0.10, 'Korea natural monument no. 53',
                  ha='center', va='top', fontsize=6.5, style='italic',
                  color=COLOR_GRAY_EDGE)


def _draw_trio_pedigree(ax_right):
    """Trio pedigree v10: COMPACT coords + REDUCED fonts."""
    # Compact coordinate system: 10x8 (was 12x9.5)
    ax_right.set_xlim(-0.3, 10.3)
    ax_right.set_ylim(-1.5, 7.5)
    # aspect 'equal' 제거: combined에서 가로로 눌려 라벨이 겹치던 문제 해결
    ax_right.axis('off')

    # (Trio family title removed per review)

    # Symbol positions (re-centered for 10-unit width)
    SIRE_X, SIRE_Y = 3.0, 5.0
    DAM_X,  DAM_Y  = 7.0, 5.0
    PROB_X, PROB_Y = 5.0, 2.2
    SYMBOL_SIZE = 0.45        # slightly smaller (was 0.55)

    # ==================== Sire ====================
    ax_right.add_patch(mpatches.Rectangle(
        (SIRE_X - SYMBOL_SIZE, SIRE_Y - SYMBOL_SIZE),
        2 * SYMBOL_SIZE, 2 * SYMBOL_SIZE,
        linewidth=1.2, edgecolor='black', facecolor='white'))
    # "Sire (paternal)" — closer to symbol
    ax_right.text(SIRE_X, SIRE_Y + SYMBOL_SIZE + 1.05,
                  'Sire (paternal)',
                  ha='center', va='bottom', fontsize=7.5, fontweight='bold')
    ax_right.text(SIRE_X, SIRE_Y + SYMBOL_SIZE + 0.30,
                  'Baekho',
                  ha='center', va='bottom', fontsize=7, fontstyle='italic')

    # ==================== Dam ====================
    ax_right.add_patch(mpatches.Circle(
        (DAM_X, DAM_Y), SYMBOL_SIZE,
        linewidth=1.2, edgecolor='black', facecolor='white'))
    ax_right.text(DAM_X, DAM_Y + SYMBOL_SIZE + 1.05,
                  'Dam (maternal)',
                  ha='center', va='bottom', fontsize=7.5, fontweight='bold')
    ax_right.text(DAM_X, DAM_Y + SYMBOL_SIZE + 0.30,
                  'Beodeul',
                  ha='center', va='bottom', fontsize=7, fontstyle='italic')

    # ==================== Mating + drop lines ====================
    ax_right.plot([SIRE_X + SYMBOL_SIZE, DAM_X - SYMBOL_SIZE],
                  [SIRE_Y, DAM_Y], color='black', linewidth=1.0)
    MID_X = (SIRE_X + DAM_X) / 2.0
    PROB_BOX = SYMBOL_SIZE * 1.1
    DROP_Y = PROB_Y + PROB_BOX + 0.25
    ax_right.plot([MID_X, MID_X], [SIRE_Y, DROP_Y],
                  color='black', linewidth=1.0)
    ax_right.plot([MID_X, PROB_X], [DROP_Y, DROP_Y],
                  color='black', linewidth=1.0)
    ax_right.plot([PROB_X, PROB_X], [DROP_Y, PROB_Y + PROB_BOX],
                  color='black', linewidth=1.0)

    # ==================== Proband ====================
    # Proband: white square with a red star = the individual we targeted/assembled
    ax_right.add_patch(mpatches.Rectangle(
        (PROB_X - PROB_BOX, PROB_Y - PROB_BOX),
        2 * PROB_BOX, 2 * PROB_BOX,
        linewidth=1.6, edgecolor='black', facecolor='white'))
    ax_right.plot(PROB_X, PROB_Y, marker='*', markersize=11,
                  color='#D62728', markeredgecolor='#D62728', zorder=5)

    # "Baeksan (offspring, F1)" — closer to proband symbol
    ax_right.text(PROB_X, PROB_Y - PROB_BOX - 0.35,
                  'Baeksan (offspring, F1)',
                  ha='center', va='top', fontsize=7.5, fontweight='bold')
    ax_right.text(PROB_X, PROB_Y - PROB_BOX - 1.00,
                  'assembled individual',
                  ha='center', va='top', fontsize=6.5, style='italic',
                  color=COLOR_GRAY_EDGE)

    # ==================== Symbol legend (bottom, well below proband) ====================
    LEG_Y = -0.95           # well below proband label
    LEG_BOX_HALF = 0.10

    # male (square)
    ax_right.add_patch(mpatches.Rectangle(
        (0.3, LEG_Y - LEG_BOX_HALF), 2 * LEG_BOX_HALF, 2 * LEG_BOX_HALF,
        linewidth=0.8, edgecolor='black', facecolor='white'))
    ax_right.text(0.65, LEG_Y, 'male', fontsize=6.5, va='center')

    # female (circle)
    ax_right.add_patch(mpatches.Circle(
        (2.6, LEG_Y), LEG_BOX_HALF,
        linewidth=0.8, edgecolor='black', facecolor='white'))
    ax_right.text(2.9, LEG_Y, 'female', fontsize=6.5, va='center')

    # assembled (white square + red star)
    ax_right.add_patch(mpatches.Rectangle(
        (5.1, LEG_Y - LEG_BOX_HALF), 2 * LEG_BOX_HALF, 2 * LEG_BOX_HALF,
        linewidth=0.8, edgecolor='black', facecolor='white'))
    ax_right.plot(5.1 + LEG_BOX_HALF, LEG_Y, marker='*', markersize=6,
                  color='#D62728', markeredgecolor='#D62728', zorder=5)
    ax_right.text(5.55, LEG_Y, 'assembled (this study)', fontsize=6.5, va='center')


# ============================================================
# Main panel function
# ============================================================
def make_panel_a(fig, gs_subplot, show_seq_strategy=True, show_panel_letter=True):
    """
    Render Panel A v10 with compact Trio pedigree.
    """
    inner_gs = gridspec.GridSpecFromSubplotSpec(
        1, 3, subplot_spec=gs_subplot,
        width_ratios=[1.0, 1.0, 1.2],   # near-equal; trio slightly wider for legend
        wspace=0.10,
    )

    ax_map = fig.add_subplot(inner_gs[0, 0], projection=ccrs.PlateCarree())
    _draw_korea_map(ax_map)

    ax_photo = fig.add_subplot(inner_gs[0, 1])
    _draw_baeksan_photo(ax_photo)

    ax_right = fig.add_subplot(inner_gs[0, 2])
    _draw_trio_pedigree(ax_right)

    if show_panel_letter:
        ax_map.text(-0.18, 1.08, 'a',
                    fontsize=14, fontweight='bold',
                    transform=ax_map.transAxes,
                    va='top', ha='left')

    if show_seq_strategy:
        fig.text(0.5, 0.04,
                 'Sequencing strategy: '
                 'PacBio HiFi (Revio) ~150x and ONT ultra-long (PromethION P2 Solo) '
                 '~103x for the F1 individual, plus RNA-seq;  '
                 'Illumina (NovaSeq 6000) ~41x for each parent and ~41x for Baeksan '
                 '(used for trio-binning k-mer phasing).',
                 ha='center', va='bottom', fontsize=8,
                 color='#404040', wrap=True)

    return ax_map, ax_photo, ax_right


# ============================================================
# Standalone main
# ============================================================
def main():
    fig = plt.figure(figsize=(11, 4.5))
    outer_gs = gridspec.GridSpec(1, 1, figure=fig,
                                  left=0.04, right=0.98,
                                  top=0.93, bottom=0.13)
    make_panel_a(fig, outer_gs[0, 0], show_seq_strategy=True, show_panel_letter=True)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_pdf = OUT_DIR / "fig1a_v11_final.pdf"
    out_png = OUT_DIR / "fig1a_v11_final.png"
    plt.savefig(out_pdf, bbox_inches='tight', dpi=300)
    plt.savefig(out_png, bbox_inches='tight', dpi=200)
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")


if __name__ == '__main__':
    main()
