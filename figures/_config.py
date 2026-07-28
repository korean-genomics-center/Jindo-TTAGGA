"""
Common configuration for all manuscript figures.
Nature Communications style guide compliant + color-blind safe.

Color palette: based on Wong (2011) Nature Methods 8:441
URL: https://www.nature.com/articles/nmeth.1618
Tested for deuteranopia, protanopia, and tritanopia.

NO red+green pairings. NO pure green for any categorical encoding.
"""
import matplotlib as mpl
from pathlib import Path

# ============================================================
# Project paths
# ============================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR     = PROJECT_ROOT / "data"
OUTPUT_DIR   = PROJECT_ROOT / "output"

# ============================================================
# Nature Communications style — fonts, sizes
# ============================================================
mpl.rcParams.update({
    'font.family':       'Arial',
    'font.size':         7,
    'axes.labelsize':    7,
    'axes.titlesize':    8,
    'xtick.labelsize':   6,
    'ytick.labelsize':   6,
    'legend.fontsize':   6,
    'legend.title_fontsize': 6,

    'pdf.fonttype':      42,
    'ps.fonttype':       42,
    'svg.fonttype':      'none',

    'axes.linewidth':    0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.minor.width': 0.4,
    'ytick.minor.width': 0.4,
    'xtick.major.size':  2.5,
    'ytick.major.size':  2.5,

    'figure.facecolor':  'white',
    'axes.facecolor':    'white',
    'savefig.facecolor': 'white',
    'savefig.dpi':       600,
    'savefig.bbox':      'tight',
    'savefig.pad_inches': 0.05,
})

# ============================================================
# Figure widths
# ============================================================
FIG_W_SINGLE = 89  / 25.4
FIG_W_HALF   = 120 / 25.4
FIG_W_DOUBLE = 183 / 25.4

# ============================================================
# Color-blind safe palette (Wong 2011 + extensions)
#
# Categorical colors used across figures:
#   #000000  black
#   #E69F00  orange
#   #56B4E9  sky blue
#   #009E73  bluish green     <-- USE WITH CAUTION (avoid red pair)
#   #F0E442  yellow
#   #0072B2  blue
#   #D55E00  vermillion       <-- "red" but distinguishable from green
#   #CC79A7  reddish purple
#
# We deliberately AVOID green-vs-red pairings throughout.
# ============================================================
COLORS = {
    # ----- Sample colors (Jindo emphasis = warm; references = cool gray) -----
    'jindo_hap1': '#D62728',   # red (Hap1, maternal, has chrX) — Taegeuk red
    'jindo_hap2': '#1F77B4',   # blue (Hap2, paternal, has chrY) — Taegeuk blue
    'ros':        '#5D6D7E',   # slate gray (ROS_Cfam1.0 reference)
    'canfam6':    '#5D6D7E',   # slate gray
    'gsd':        '#85929E',   # silver gray
    'canfam3':    '#B2BABB',   # light gray

    # ----- BUSCO categories (4 categories, all distinguishable for CVD) -----
    # Strategy: blue gradient for "good" (Complete), orange for problem (Fragmented),
    # vermillion for serious problem (Missing).
    # NO green anywhere.
    'busco_C_S':  '#0072B2',   # complete single-copy   - dark blue
    'busco_C_D':  '#56B4E9',   # complete duplicated    - sky blue
    'busco_F':    '#F0E442',   # fragmented             - yellow
    'busco_M':    '#D55E00',   # missing                - vermillion

    # ----- Telomere (forward / reverse) -----
    'telomere_fwd': '#0072B2',   # blue
    'telomere_rev': '#CC79A7',   # reddish purple (NOT pure red)
    'gap':          '#000000',

    # ----- SV categories (Fig 4D synteny: inv / trans / dup) -----
    'sv_syntenic':     '#B2BABB',   # gray
    'sv_inversion':    '#E69F00',   # orange
    'sv_translocation':'#56B4E9',   # sky blue
    'sv_duplication':  '#CC79A7',   # reddish purple

    # ----- Repeat element categories (Fig 2A) — needs many distinguishable colors -----
    'repeat_LINE':    '#0072B2',
    'repeat_SINE':    '#56B4E9',
    'repeat_LTR':     '#E69F00',
    'repeat_DNA_TE':  '#CC79A7',
    'repeat_simple':  '#F0E442',
    'repeat_satellite':'#D55E00',
    'repeat_other':   '#B2BABB',
}

# ============================================================
# Standard chromosome ordering helper
# ============================================================
def chrom_sort_key(c):
    """Sort key for canine chromosomes: chr1..chr38 < chrX < chrY < others."""
    s = str(c).replace("chr", "")
    if s == "X":
        return (1, 0)
    if s == "Y":
        return (2, 0)
    try:
        return (0, int(s))
    except ValueError:
        return (3, s)
