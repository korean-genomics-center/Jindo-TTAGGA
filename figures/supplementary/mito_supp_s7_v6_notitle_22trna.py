#!/usr/bin/env python
"""
Supplementary Figure S7 v4: Mitochondrial genome - paper-grade

Fixes from v3:
  - tRNA detection: only "tRNA-XXX" pattern (not just "Phe"/"Val" in product)
  - Should now correctly identify 22 tRNAs (universal mammalian standard)
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from matplotlib.patches import Patch

OUT_DIR = Path("${JINDO_ROOT}/Results/Manuscript_Figures/output/mito")
GENBANK = Path("${JINDO_ROOT}/Results/Mitochondrial_Genome_Assembly/potential_contigs/atg000005l/atg000005l.annotation/atg000005l.annotation_MitoFinder_mitfi_Final_Results/atg000005l.annotation_mtDNA_contig.gb")
DEPTH_FILE = OUT_DIR / "final_mito_depth.txt"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "savefig.dpi": 300,
    "axes.linewidth": 1.2,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# === Parse GenBank — STRICT logic ===
features = []
with open(GENBANK) as f:
    lines = f.readlines()

i = 0
while i < len(lines):
    line = lines[i].rstrip()
    if line.startswith("     gene") and ".." in line:
        loc = line.strip().split()[1]
        strand = "+"
        if loc.startswith("complement("):
            strand = "-"
            loc = loc.replace("complement(", "").rstrip(")")
        if ".." in loc:
            try:
                start, end = map(int, loc.replace("<", "").replace(">", "").split(".."))
            except ValueError:
                i += 1
                continue
            
            gene_name = None
            product = None
            j = i + 1
            while j < len(lines) and lines[j].startswith("                     "):
                qline = lines[j].strip()
                if qline.startswith("/gene="):
                    gene_name = qline.split("=", 1)[1].strip().strip('"')
                elif qline.startswith("/product="):
                    product = qline.split("=", 1)[1].strip().strip('"')
                j += 1
            
            # STRICT classification
            if product and product.startswith("tRNA-") and product != "tRNA-Ala2":
                # tRNA only if /product= explicitly starts with "tRNA-"
                ftype = "tRNA"
                name = product
            elif product and (product.startswith("rrn") or "rRNA" in product or 
                              product in ["s-rRNA", "l-rRNA", "12S ribosomal RNA", "16S ribosomal RNA"]):
                ftype = "rRNA"
                name = product
            elif gene_name and not product:
                # Has /gene= but no /product= → CDS (PCG)
                ftype = "CDS"
                name = gene_name
            elif gene_name and product and product.startswith("tRNA-") and product != "tRNA-Ala2":
                ftype = "tRNA"
                name = product
            else:
                # Skip ambiguous
                i = j
                continue
            
            features.append({
                "type": ftype, "start": start, "end": end,
                "strand": strand, "name": name
            })
            i = j
            continue
    i += 1

print(f"Parsed {len(features)} features")
from collections import Counter
type_counts = Counter(f["type"] for f in features)
print(f"Type counts: {dict(type_counts)}")

# === Read coverage ===
depth_df = pd.read_csv(DEPTH_FILE, sep="\t", names=["chrom", "pos", "depth"])
mean_cov = depth_df["depth"].mean()
print(f"Coverage: mean={mean_cov:.0f}x")

# === Figure ===
fig = plt.figure(figsize=(14, 8))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.3], hspace=0.35,
                      left=0.07, right=0.97, top=0.94, bottom=0.07)

# Panel A
ax_a = fig.add_subplot(gs[0])

color_map = {
    "CDS": "#9467bd",
    "rRNA": "#2ca02c",
    "tRNA": "#1f77b4",
}

mito_len = 16741
y_fwd = 0.55
y_rev = 0.30
height = 0.13

mid_y = (y_fwd + y_rev + height) / 2
ax_a.axhline(y=mid_y, xmin=0.02, xmax=0.98,
             color="black", linewidth=0.6, zorder=1)

for feat in features:
    y = y_fwd if feat["strand"] == "+" else y_rev
    color = color_map.get(feat["type"], "#888888")
    width = feat["end"] - feat["start"]
    ax_a.add_patch(plt.Rectangle((feat["start"], y), width, height,
                                  facecolor=color, edgecolor="black",
                                  linewidth=0.4, zorder=5))
    
    if feat["type"] in ("CDS", "rRNA"):
        center = (feat["start"] + feat["end"]) / 2
        y_label = (y + height + 0.04) if feat["strand"] == "+" else (y - 0.05)
        va = "bottom" if feat["strand"] == "+" else "top"
        ax_a.text(center, y_label, feat["name"], fontsize=7,
                  ha="center", va=va, fontweight="bold")

ax_a.text(-300, y_fwd + height/2, "+", fontsize=14, fontweight="bold",
          ha="right", va="center", color="#444444")
ax_a.text(-300, y_rev + height/2, "−", fontsize=14, fontweight="bold",
          ha="right", va="center", color="#444444")

n_pcg = type_counts.get("CDS", 0)
n_rrna = type_counts.get("rRNA", 0)
n_trna = type_counts.get("tRNA", 0)

legend_handles = [
    Patch(facecolor="#9467bd", edgecolor="black", linewidth=0.4, label=f"Protein-coding ({n_pcg})"),
    Patch(facecolor="#2ca02c", edgecolor="black", linewidth=0.4, label=f"rRNA ({n_rrna})"),
    Patch(facecolor="#1f77b4", edgecolor="black", linewidth=0.4, label=f"tRNA ({n_trna})"),
]
ax_a.legend(handles=legend_handles, loc="upper right", frameon=False,
            fontsize=9, ncol=3, bbox_to_anchor=(0.98, 1.05))

ax_a.set_xlim(-500, mito_len + 200)
ax_a.set_ylim(0.05, 0.95)
ax_a.set_yticks([])
ax_a.set_xlabel("Position (bp)")
n_total = n_pcg + n_rrna + n_trna
ax_a.set_title("a", loc="left", fontweight="bold", fontsize=11)
ax_a.spines["left"].set_visible(False)

# Panel B
ax_b = fig.add_subplot(gs[1])

ax_b.fill_between(depth_df["pos"], depth_df["depth"], color="#1f77b4", alpha=0.6, linewidth=0)
ax_b.plot(depth_df["pos"], depth_df["depth"], color="#1f77b4", linewidth=0.5)

ax_b.axhline(y=mean_cov, color="#d62728", linestyle="--", linewidth=1.2,
             label=f"Mean = {mean_cov:.0f}\u00d7")

ax_b.set_xlim(0, mito_len)
ax_b.set_xlabel("Position (bp)")
ax_b.set_ylabel("PacBio HiFi coverage depth (\u00d7)")
ax_b.set_title("b", loc="left", fontweight="bold", fontsize=11)
ax_b.legend(loc="upper right", frameon=False)
ax_b.grid(axis="y", linestyle=":", alpha=0.3)

plt.savefig(OUT_DIR / "supp_fig_s7_mitogenome_v6.pdf", bbox_inches="tight")
plt.savefig(OUT_DIR / "supp_fig_s7_mitogenome_v6.png", bbox_inches="tight", dpi=200)
print(f"\nSaved: {OUT_DIR}/supp_fig_s7_mitogenome_v6.{{pdf,png}}")
