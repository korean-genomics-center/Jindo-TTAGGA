#!/usr/bin/env python3
"""Supplementary Figure 8: recorded paternal lines of the 20 cohort males, grouped by chrY lineage.
Source: official pedigree certificates issued by Jindo County (18 of 20 males).
Ancestor names and owners are withheld; cohort males are shown by sample ID.
chrY lineage and paternal-line codes as in Supplementary Data 7 (sample, paternal_lineage, paternal_line)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

A_COL, B_COL = "#1F77B4", "#D45B00"   # chrY lineage colours as in Fig. 5
S = 0.42
G = {0:5.0, 1:3.9, 2:2.8, 3:1.7}
YS = 0.25                              # singleton row

fig, ax = plt.subplots(figsize=(7.4, 4.6))

def sq(x, y, label=None, lin=None, cert=True):
    fc = {"A":A_COL, "B":B_COL}.get(lin, "white")
    ax.add_patch(Rectangle((x-S/2, y-S/2), S, S, facecolor=fc, edgecolor="black",
                           lw=0.9, ls="-" if cert else (0,(2,1.5)), zorder=3,
                           hatch=None if cert else "////"))
    if label:
        ax.text(x, y-S/2-0.1, label, ha="center", va="top", fontsize=6.6, fontweight="bold", zorder=4)

def link(parent, kids):
    px, py = parent; ymid = py - 0.55
    ax.plot([px, px], [py-S/2, ymid], color="black", lw=0.9, zorder=1)
    xs = [k[0] for k in kids]
    ax.plot([min(xs+[px]), max(xs+[px])], [ymid, ymid], color="black", lw=0.9, zorder=1)
    for kx, ky in kids:
        ax.plot([kx, kx], [ymid, ky+S/2], color="black", lw=0.9, zorder=1)

def tag(x, y, t):
    ax.text(x, y, t, ha="center", va="bottom", fontsize=7, color="#444444", style="italic")

# ---------------- Lineage A (left) ----------------
# PL2: two grandsons of a common paternal ancestor
f=(1.3,G[0]); a=(0.7,G[1]); b=(1.9,G[1]); c1=(0.7,G[2]); c2=(1.9,G[2])
sq(*f); sq(*a); sq(*b); sq(*c1,"J352388","A"); sq(*c2,"J498404","A")
link(f,[a,b]); link(a,[c1]); link(b,[c2]); tag(1.3,G[0]+0.32,"Paternal line 2")
# PL3: son and great-grandson
f=(4.2,G[0]); k1=(3.6,G[1]); m1=(4.8,G[1]); m2=(4.8,G[2]); k2=(4.8,G[3])
sq(*f); sq(*k1,"J13044625","A"); sq(*m1); sq(*m2); sq(*k2,"J63253","A")
link(f,[k1,m1]); link(m1,[m2]); link(m2,[k2]); tag(4.2,G[0]+0.32,"Paternal line 3")
# singletons A
sA=["J12952315","J288508","J496128","J61783","J62428"]
for i,s in enumerate(sA): sq(0.4+i*1.12, YS, s, "A")
sq(0.4+5*1.12, YS, "J76290", "A", cert=False)

# ---------------- Lineage B (right) ----------------
X0=7.6
f=(X0+1.5,G[0]); sire=(X0+0.8,G[1]); s1=(X0+0.2,G[2]); s2=(X0+1.4,G[2]); s3=(X0+2.3,G[1])
sq(*f); sq(*sire); sq(*s1,"J433555","B"); sq(*s2,"J63542","B"); sq(*s3,"J543672","B")
link(f,[sire,s3]); link(sire,[s1,s2]); tag(X0+1.5,G[0]+0.32,"Paternal line 1")
sB=["J288326","J495439","J61435","J62957","J73478","J90474"]
for i,s in enumerate(sB): sq(X0+0.05+i*0.98, YS, s, "B")
sq(X0+0.05+6*0.98, YS, "J73201", "B", cert=False)

# column headers, divider and row label
ax.text(3.0, 6.05, "chrY lineage A (10 males)", ha="center", fontsize=9, fontweight="bold", color=A_COL)
ax.text(X0+2.9, 6.05, "chrY lineage B (10 males)", ha="center", fontsize=9, fontweight="bold", color=B_COL)
ax.plot([7.0,7.0],[-0.3,6.2], color="#BBBBBB", lw=0.8, ls=(0,(4,3)))
ax.plot([-0.3,14.2],[0.95,0.95], color="#DDDDDD", lw=0.6)
ax.text(-0.3, 0.98, "No other cohort male on the recorded paternal line", fontsize=6.5, color="#444444", va="bottom", style="italic")

handles=[Rectangle((0,0),1,1,facecolor=A_COL,edgecolor="black",label="Cohort male, lineage A"),
         Rectangle((0,0),1,1,facecolor=B_COL,edgecolor="black",label="Cohort male, lineage B"),
         Rectangle((0,0),1,1,facecolor="white",edgecolor="black",label="Recorded paternal ancestor (not sequenced)"),
         Rectangle((0,0),1,1,facecolor="#999999",edgecolor="black",hatch="////",ls=(0,(2,1.5)),label="No pedigree certificate")]
fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=6.6, frameon=False, bbox_to_anchor=(0.5,0.0), handlelength=1.4, columnspacing=1.2)
ax.set_xlim(-0.4,14.3); ax.set_ylim(-0.35,6.35); ax.set_aspect("equal"); ax.axis("off")
plt.tight_layout(rect=(0,0.03,1,1))
OUT = "${JINDO_ROOT}/Results/Manuscript_Figures/Supplementary/figures"
plt.savefig(f"{OUT}/Sfig8_pedigree_v4.pdf", bbox_inches="tight"); plt.savefig(f"{OUT}/Sfig8_pedigree_v4.png", dpi=300, bbox_inches="tight")
print("saved")
