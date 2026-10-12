"""
Fig 5 Panel a v1: Canid BUSCO single-copy ortholog phylogeny.

Refactored from Render_phylogeny_paper_grade3.py (Sup S8 base) to importable
make_panel_a form for combined Fig 5.

Narrative: Jindo Hap1/Hap2 cluster as the most wolf-proximal domestic lineage
across 33 canid taxa (12,277 BUSCO single-copy orthologs; MACSE+trimAl+RAxML-NG+
ASTRAL species tree).

Key finding (this study): Jindo diverges immediately after the wolf clade,
basal to all other domestic dog breeds including the Dingo. Both maternal
(Hap1, chrX) and paternal (Hap2, chrY) haplotypes share this position.

Refactor changes from Sup S8 v3:
  - Function form: make_panel_a(fig, ax, show_title=False, show_legend=True)
  - NC compact layout (no large title; main figure font sizes)
  - Wolf-proximal annotation arrow on Jindo position
  - Reduced label font size for double-column width
"""
import os
import sys
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

SCRIPTS_DIR = Path("${JINDO_ROOT}/Results/Manuscript_Figures/scripts")
sys.path.insert(0, str(SCRIPTS_DIR))
from _config import FIG_W_DOUBLE


# ============================================================
# Paths
# ============================================================
ROOT = Path("${JINDO_ROOT}")
TREE_PATH = ROOT / "Results/Phylogeny_BUSCO/From_whale/BUSCO/Merged_aln_parallel/trimmed/rename_trimmed/raxml_out/astral.species0.set_outgroup.tre"
OUT_DIR = ROOT / "Results/Manuscript_Figures/output/fig5"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# Taxon map (same as Sup S8 v3)
# ============================================================
TAXON_MAP = {
    "GCA034620435hap1": ("Woyote (wolf × coyote)",                     "outgroup"),
    "GCA905319855hap2": ("Greenland wolf",                              "wolf"),
    "GCA007922845hap1": ("Gray wolf (UniMelb)",                         "wolf"),
    "GCA048164845hap1": ("Mexican gray wolf (Hap2)",                    "wolf"),
    "GCA048164855hap1": ("Mexican gray wolf (Hap1)",                    "wolf"),
    "GCF048164855hap1": ("Mexican gray wolf (RefSeq)",                  "wolf"),
    "GCAJindoDoghap1":  ("Jindo dog Hap1 (maternal, chrX)",             "jindo"),
    "GCAJindoDoghap2":  ("Jindo dog Hap2 (paternal, chrY)",             "jindo"),
    "GCA003254725hap2": ("Australian Dingo (Sandy)",                    "dingo"),
    "GCF003254725hap2": ("Australian Dingo (RefSeq)",                   "dingo"),
    "GCA012295265hap2": ("Alpine Dingo (Cooindah)",                     "dingo"),
    "GCA000002285hap4": ("Boxer (Tasha, CanFam3.1)",                    "western"),
    "GCF000002285hap5": ("Boxer (Tasha, RefSeq)",                       "western"),
    "GCA014441545hap1": ("Labrador (ROS_Cfam_1.0)",                     "western"),
    "GCF014441545hap1": ("Labrador (ROS_Cfam, RefSeq)",                 "western"),
    "GCA031165255hap1": ("Labrador (Yella v2)",                         "western"),
    "GCA011100685hap1": ("German Shepherd (Mischka)",                   "western"),
    "GCF011100685hap1": ("German Shepherd (Mischka, RefSeq)",           "western"),
    "GCA008641055hap3": ("German Shepherd (Nala)",                      "western"),
    "GCA004886185hap2": ("Basenji (breed-1.1)",                         "western"),
    "GCA013276365hap2": ("Basenji (CanFam_Bas)",                        "western"),
    "GCA005444595hap1": ("Great Dane (Zoey)",                           "western"),
    "GCF005444595hap1": ("Great Dane (Zoey, RefSeq)",                   "western"),
    "GCA044048985hap1": ("Beagle (near-T2T 2025)",                      "western"),
    "GCA031010635hap1": ("Bernese Mountain (OD)",                       "western"),
    "GCA031010765hap1": ("Bernese Mountain (BD)",                       "western"),
    "GCA031010295hap1": ("Cairn Terrier (CA611)",                       "western"),
    "GCA031771975hap1": ("Irish Wolfhound (N220234)",                   "western"),
    "GCA040939265hap1": ("Irish Wolfhound (N210636)",                   "western"),
    "GCA043643935hap1": ("Whippet",                                     "western"),
    "GCA004027395hap1": ("Dog (BS72)",                                  "western"),
    "GCA012044875hap1": ("Dog (JHU/APL 1)",                             "western"),
    "GCA012045015hap1": ("Dog (JHU/APL 2)",                             "western"),
}

GROUP_COLORS = {
    "outgroup": "#7F0000",
    "wolf":     "#D55E00",
    "jindo":    "#0072B2",   # Wong blue (Fig 2 consistency)
    "dingo":    "#009E73",   # Wong green
    "western":  "#666666",
}


# ============================================================
# Newick parser
# ============================================================
class Node:
    __slots__ = ("name", "branch_length", "support", "children", "x", "y")
    def __init__(self, name="", branch_length=0.0, support=None):
        self.name = name
        self.branch_length = branch_length
        self.support = support
        self.children = []
        self.x = 0.0
        self.y = 0.0


def parse_newick(s):
    s = s.strip().rstrip(";").strip()
    pos = [0]
    def parse_subtree():
        node = Node()
        if s[pos[0]] == "(":
            pos[0] += 1
            while True:
                node.children.append(parse_subtree())
                if s[pos[0]] == ",":
                    pos[0] += 1
                elif s[pos[0]] == ")":
                    pos[0] += 1
                    break
        label = ""
        while pos[0] < len(s) and s[pos[0]] not in ":,()":
            label += s[pos[0]]; pos[0] += 1
        if node.children:
            node.support = label if label else None
        else:
            node.name = label
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            bl = ""
            while pos[0] < len(s) and s[pos[0]] not in ",()":
                bl += s[pos[0]]; pos[0] += 1
            try: node.branch_length = float(bl)
            except: node.branch_length = 0.0
        return node
    return parse_subtree()


def get_leaves(node):
    if not node.children: return [node]
    out = []
    for c in node.children: out.extend(get_leaves(c))
    return out


def layout(node, x=0.0, y_counter=[0]):
    node.x = x
    if not node.children:
        node.y = y_counter[0]; y_counter[0] += 1; return
    for c in node.children:
        layout(c, x + c.branch_length, y_counter)
    node.y = sum(c.y for c in node.children) / len(node.children)


def reorder(node):
    for c in node.children: reorder(c)
    if node.children:
        def rank(c):
            leaves = get_leaves(c)
            groups = set(TAXON_MAP.get(l.name, (l.name, "western"))[1] for l in leaves)
            if "outgroup" in groups: return 0
            if "wolf" in groups and "jindo" not in groups: return 1
            if "jindo" in groups: return 2
            if "dingo" in groups: return 3
            return 4
        node.children.sort(key=rank)


# ============================================================
# Main panel function
# ============================================================
def make_panel_a(fig, gs_subplot=None, ax=None, show_title=True, show_legend=True):
    """
    Render Panel a (33-canid BUSCO phylogeny) into given gridspec subplot.
    Compact layout for NC double-column main figure.
    """
    if ax is None:
        if gs_subplot is None:
            ax = fig.add_subplot(111)
        else:
            ax = fig.add_subplot(gs_subplot)

    # Parse tree
    with open(TREE_PATH) as f:
        nwk = f.read().strip()
    nwk_clean = re.sub(r"\[[^\]]*\]", "", nwk)
    root = parse_newick(nwk_clean)
    reorder(root)
    # Reset y_counter for repeatable layout
    layout.__defaults__ = (0.0, [0])
    layout(root)

    leaves = get_leaves(root)
    n_leaves = len(leaves)

    # Tree edges
    def plot_node(node):
        if node.children:
            ys = [c.y for c in node.children]
            ax.plot([node.x, node.x], [min(ys), max(ys)],
                    color="#222", lw=0.6, zorder=2)
            for c in node.children:
                ax.plot([node.x, c.x], [c.y, c.y],
                        color="#222", lw=0.6, zorder=2)
                plot_node(c)
            # Support values (rare display)
            if node.support and node.support not in (None, ""):
                try:
                    sup = float(node.support)
                    if 0.5 <= sup < 0.95:  # show only ambiguous nodes
                        ax.text(node.x, node.y, f"{sup:.2f}",
                                fontsize=4.5, color="#aaa",
                                ha="left", va="center", zorder=3)
                except ValueError:
                    pass

    plot_node(root)

    # Leaf labels
    leaf_xmax = max(l.x for l in leaves)
    for leaf in leaves:
        friendly, group = TAXON_MAP.get(leaf.name, (leaf.name, "western"))
        col = GROUP_COLORS[group]
        weight = "bold" if group == "jindo" else "normal"
        fs = 5.5 if group == "jindo" else 5.0
        ax.text(leaf.x + leaf_xmax * 0.008, leaf.y,
                "  " + friendly,
                ha="left", va="center",
                fontsize=fs, color=col, fontweight=weight, zorder=4,
                family="Arial")

    # ============================================================
    # Highlight Jindo position with annotation
    # ============================================================
    jindo_leaves = [l for l in leaves if TAXON_MAP.get(l.name, (None, None))[1] == "jindo"]
    if jindo_leaves:
        jindo_ys = [l.y for l in jindo_leaves]
        jindo_y_center = sum(jindo_ys) / len(jindo_ys)
        jindo_x = max(l.x for l in jindo_leaves)

        # Highlight box around Jindo labels
        ax.add_patch(Rectangle(
            (jindo_x + leaf_xmax * 0.005, min(jindo_ys) - 0.4),
            leaf_xmax * 0.55, max(jindo_ys) - min(jindo_ys) + 0.8,
            facecolor="#E0F0FF", edgecolor="#0072B2",
            linewidth=0.6, alpha=0.4, zorder=1
        ))

        # [NC-style: annotation removed] # Annotation arrow pointing to Jindo
        # [NC-style: annotation removed] ax.annotate(
        # [NC-style: annotation removed] "Most wolf-proximal\ndomestic lineage\n(basal to Dingo)",
        # [NC-style: annotation removed] xy=(jindo_x + leaf_xmax * 0.6, jindo_y_center),
        # [NC-style: annotation removed] xytext=(leaf_xmax * 1.55, jindo_y_center),
        # [NC-style: annotation removed] ha="left", va="center",
        # [NC-style: annotation removed] fontsize=6, color="#0072B2",
        # [NC-style: annotation removed] fontweight="bold",
        # [NC-style: annotation removed] arrowprops=dict(arrowstyle="->",
        # [NC-style: annotation removed] color="#0072B2", lw=0.8),
        # [NC-style: annotation removed] bbox=dict(boxstyle="round,pad=0.4",
        # [NC-style: annotation removed] facecolor="#F0F8FF",
        # [NC-style: annotation removed] edgecolor="#0072B2", linewidth=0.6),
        # [NC-style: annotation removed] )

    # Axes
    ax.set_yticks([])
    ax.set_ylim(-1, n_leaves)
    ax.invert_yaxis()
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", labelsize=5.5)
    ax.set_xlabel("ASTRAL coalescent units (terminal branches = 1; not substitution rate)",
                  fontsize=6.5)
    ax.set_xlim(-leaf_xmax * 0.02, leaf_xmax * 2.3)

    if show_title:
        ax.set_title(
            "Jindo clusters as the most wolf-proximal domestic lineage\n"
            "(33 canid taxa, 12,277 BUSCO orthologs)",
            fontsize=8, pad=8, loc="left", fontweight="bold"
        )

    if show_legend:
        legend_elems = [
            Rectangle((0, 0), 1, 1, fc=GROUP_COLORS["jindo"],
                      label="Jindo (this study)"),
            Rectangle((0, 0), 1, 1, fc=GROUP_COLORS["wolf"],
                      label="Wolf"),
            Rectangle((0, 0), 1, 1, fc=GROUP_COLORS["dingo"],
                      label="Dingo"),
            Rectangle((0, 0), 1, 1, fc=GROUP_COLORS["western"],
                      label="Other domestic dog breeds"),
            Rectangle((0, 0), 1, 1, fc=GROUP_COLORS["outgroup"],
                      label="Outgroup (Woyote)"),
        ]
        ax.legend(handles=legend_elems,
                  loc="lower right",
                  bbox_to_anchor=(1.0, 0.0),
                  fontsize=5.5, frameon=True, framealpha=0.95,
                  edgecolor="#888",
                  ncol=1,
                  handlelength=1.2,
                  handleheight=1.0)

    return ax


# ============================================================
# Standalone main
# ============================================================
def main():
    fig_w = FIG_W_DOUBLE
    fig_h = fig_w * 0.85  # taller for 33 taxa stack
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    make_panel_a(fig, ax=ax, show_title=False)

    plt.tight_layout()

    out_pdf = OUT_DIR / "Sfig5_phylogeny_v3.pdf"
    out_png = OUT_DIR / "Sfig5_phylogeny_v3.png"
    plt.savefig(out_pdf, dpi=600, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")

    print()
    print("=== Panel a v1 Summary ===")
    print(f"  33 canid taxa, BUSCO single-copy ortholog ASTRAL phylogeny")
    print(f"  Jindo Hap1/Hap2 highlighted as wolf-proximal domestic lineage")
    print(f"  Basal to Dingo, all other domestic dog breeds")


if __name__ == "__main__":
    main()
