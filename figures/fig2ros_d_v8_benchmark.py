"""
Fig 2D v7 (BUSCO removed): Multi-metric assembly benchmark — Jindo Hap1/Hap2
vs three canine references (ROS_Cfam1.0, CanFam6, GSD1.0).

Four panels:
    1. Total length (Gb)
    2. Number of scaffolds   (fasta entries)
    3. Contig N50 (Mb)       (contigs split at gap definition)
    4. Number of gaps        (>=100 bp N stretches per assembly, NC standard)

Changes from v3:
  - Removed BUSCO panel (moved to Fig 1C: BUSCO + QV + GCI panel)
  - 5 metrics -> 4 metrics
  - Renamed: fig1c_v3 -> fig2ros_d_v7

Visual emphasis: Jindo (warm Wong-palette) vs References (uniform gray).

Changes from v2 (fig1c_benchmark.py, 2026-04-29):
  - Gap count source: gap_stats_v2.tsv (NEW, >=100 bp N stretch)
    instead of gap_stats.tsv (old, >=10 N stretch).
    Rationale: aligns with NC/T2T-CHM13 convention (>=100 bp for "true gaps")
    and with Panel E (Fig 1) and Fig 3 Panel D measurement.
  - Refactored into make_panel_c(fig, gs_subplot) function for use in
    combined Fig 1 figure.
  - Output filename versioned to preserve v2 output.

Usage:
    Standalone: python fig1c_v3_benchmark_refactored.py
    Imported:   from fig1c_v3_benchmark_refactored import make_panel_c
                make_panel_c(fig, gs_subplot)

Inputs:
   data/stats/gap_stats_v2.tsv (NEW, >=100 bp N stretch)
   data/busco/{JindoHap1,JindoHap2,ROS_Cfam1.0,CanFam6,GSD1.0}.carnivora_odb12.txt

Output: output/fig1/panels/fig2ros_d_v7.{pdf,png}
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _config import COLORS, FIG_W_DOUBLE, DATA_DIR, OUTPUT_DIR
from _utils import save_fig, parse_busco_short_summary


# ============================================================
# Sample list — preserved from v2
# ============================================================
REF_GRAY = '#999999'

SAMPLES = [
    ('JindoHap1',   'Jindo Hap1\n(maternal)',  'JindoHap1.carnivora_odb12.txt',  '#D62728'),
    ('JindoHap2',   'Jindo Hap2\n(paternal)',  'JindoHap2.carnivora_odb12.txt',  '#1F77B4'),
    ('ROS_Cfam1.0', 'ROS_Cfam_1.0',             'ROS_Cfam1.0.carnivora_odb12.txt', REF_GRAY),
    ('CanFam6',     'Dog10K_Boxer_Tasha',                 'CanFam6.carnivora_odb12.txt',    REF_GRAY),
    ('GSD1.0',      'UU_Cfam_GSD_1.0',                  'GSD1.0.carnivora_odb12.txt',     REF_GRAY),
]


# ============================================================
# Cached data loading
# ============================================================
_data_cache = {}

def _load_data():
    if "loaded" in _data_cache:
        return _data_cache

    print("[Panel C] Loading benchmark data ...")

    # 1. Gap stats (v2: >=100 bp N stretch, NC standard)
    gap_path = DATA_DIR / 'stats' / 'gap_stats_v2.tsv'
    has_gaps = gap_path.exists()
    if has_gaps:
        gap_df = pd.read_csv(gap_path, sep='\t').set_index('sample')
        print(f"  gap_stats_v2 loaded: {list(gap_df.index)}")
    else:
        gap_df = None
        print("  WARN: gap_stats_v2.tsv not found")

    # 2. BUSCO
    busco_data = {}
    for name, label, busco_file, _ in SAMPLES:
        f = DATA_DIR / 'busco' / busco_file
        d = parse_busco_short_summary(f)
        busco_data[name] = d
        print(f"  {label.replace(chr(10), ' ')}: BUSCO C={d['C']:.1f}%")

    _data_cache.update({
        "gap_df": gap_df,
        "busco_data": busco_data,
        "has_gaps": has_gaps,
        "loaded": True,
    })
    return _data_cache


# ============================================================
# Compute metrics from cached data
# ============================================================
def _compute_metrics(data):
    busco_data = data["busco_data"]
    gap_df = data["gap_df"]
    has_gaps = data["has_gaps"]

    keys = [s[0] for s in SAMPLES]

    # Length: from BUSCO short summary (total_length)
    total_lengths_gb = []
    for k in keys:
        # We need total length from somewhere — use BUSCO summary
        # If BUSCO file doesn't have length, fall back to gap_df
        if has_gaps and 'total_length_bp' in gap_df.columns:
            total_lengths_gb.append(gap_df.loc[k, 'total_length_bp'] / 1e9)
        else:
            # Use known values (Jindo, references)
            length_lookup = {
                'JindoHap1':   2441.6e6,
                'JindoHap2':   2340.5e6,
                'ROS_Cfam1.0': 2396.9e6,
                'CanFam6':     2312.8e6,
                'GSD1.0':      2482.0e6,
            }
            total_lengths_gb.append(length_lookup.get(k, 0) / 1e9)

    # Number of scaffolds: from gap_stats (contig_count)
    if has_gaps:
        n_scaffolds = [int(gap_df.loc[k, 'contig_count']) for k in keys]
        n_gaps      = [int(gap_df.loc[k, 'gap_count']) for k in keys]
        n50_mb      = [gap_df.loc[k, 'contig_n50_bp'] / 1e6 for k in keys]
    else:
        n_scaffolds = None
        n_gaps      = None
        n50_mb      = None

    # BUSCO C
    busco_c = [busco_data[k]['C'] for k in keys]

    return {
        "labels":      [s[1] for s in SAMPLES],
        "colors":      [s[3] for s in SAMPLES],
        "total_gb":    total_lengths_gb,
        "n_scaffolds": n_scaffolds,
        "n50_mb":      n50_mb,
        "n_gaps":      n_gaps,
        "busco_c":     busco_c,
    }


# ============================================================
# Drawing
# ============================================================
def _draw_metric(ax, title, values, fmt, log_x, labels, colors, fontsize=7):
    n = len(labels)
    y_positions = np.arange(n)[::-1]

    if values is None:
        ax.text(0.5, 0.5, 'N/A', ha='center', va='center',
                fontsize=8, transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks(y_positions)
        ax.set_yticklabels(labels, fontsize=fontsize)
        ax.set_title(title, fontsize=8, fontweight='bold')
        return

    bars = ax.barh(y_positions, values, height=0.65, color=colors,
                   edgecolor='white', linewidth=0.4)

    # Annotate values inline
    for i, (y, v) in enumerate(zip(y_positions, values)):
        if v == 0:
            label_x = 0
        else:
            label_x = v * 1.05 if not log_x else v * 1.4
        formatted = fmt.format(v)
        # For Jindo highlight: bold + color match
        is_jindo = i < 2
        ax.text(label_x, y, formatted,
                va='center', ha='left',
                fontsize=fontsize - 0.5,
                fontweight='bold' if is_jindo else 'normal',
                color='#000')

    if log_x:
        ax.set_xscale('log')

    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=fontsize)
    ax.set_title(title, fontsize=8, fontweight='bold', pad=4)

    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.tick_params(axis='x', length=2, width=0.5, labelsize=6.5)
    ax.tick_params(axis='y', length=0)


def _draw_panel_c(axes, metrics):
    ax_len, ax_scf, ax_n50, ax_gaps = axes

    # 1. Total length
    _draw_metric(ax_len, 'Total length (Gb)',
                 metrics["total_gb"], '{:.2f}',
                 log_x=False,
                 labels=metrics["labels"], colors=metrics["colors"])
    ax_len.set_xlim(0, max(metrics["total_gb"]) * 1.25)

    # 2. Number of scaffolds (log scale)
    _draw_metric(ax_scf, 'Number of scaffolds',
                 metrics["n_scaffolds"], '{:,}',
                 log_x=True,
                 labels=[''] * len(metrics["labels"]),
                 colors=metrics["colors"])
    if metrics["n_scaffolds"]:
        ax_scf.set_xlim(left=10, right=max(metrics["n_scaffolds"]) * 5)

    # 3. Contig N50 (Mb)
    _draw_metric(ax_n50, 'Contig N50 (Mb)',
                 metrics["n50_mb"], '{:.1f}',
                 log_x=False,
                 labels=[''] * len(metrics["labels"]),
                 colors=metrics["colors"])
    if metrics["n50_mb"]:
        ax_n50.set_xlim(0, max(metrics["n50_mb"]) * 1.30)

    # 4. Number of gaps (log scale, since 0..585)
    if metrics["n_gaps"] is not None:
        # For log scale, replace 0 with small value (display only)
        gaps_for_plot = [max(v, 0.5) for v in metrics["n_gaps"]]
        _draw_metric(ax_gaps, 'Number of gaps',
                     gaps_for_plot, '{:,.0f}',
                     log_x=True,
                     labels=[''] * len(metrics["labels"]),
                     colors=metrics["colors"])
        # Override formatted text for the 0-bar (would show "1" otherwise)
        for txt, real_v in zip(ax_gaps.texts, metrics["n_gaps"]):
            if real_v == 0:
                txt.set_text('0')
        ax_gaps.set_xlim(left=0.5, right=max(metrics["n_gaps"]) * 5)
    else:
        _draw_metric(ax_gaps, 'Number of gaps',
                     None, '{:,}', log_x=False,
                     labels=[''] * len(metrics["labels"]),
                     colors=metrics["colors"])

# ============================================================
# Main panel function
# ============================================================
def make_panel_c(fig, gs_subplot, show_footnote=True):
    """
    Render Panel C (5-metric benchmark) into the given gridspec subplot.

    Layout: 4 sub-axes side-by-side.

    Returns
    -------
    (ax_len, ax_scf, ax_n50, ax_gaps) : tuple of Axes
    """
    data = _load_data()
    metrics = _compute_metrics(data)

    inner_gs = gridspec.GridSpecFromSubplotSpec(
        1, 4, subplot_spec=gs_subplot,
        width_ratios=[1.4, 1.0, 1.0, 1.0],
        wspace=0.50,
    )
    ax_len   = fig.add_subplot(inner_gs[0, 0])
    ax_scf   = fig.add_subplot(inner_gs[0, 1])
    ax_n50   = fig.add_subplot(inner_gs[0, 2])
    ax_gaps  = fig.add_subplot(inner_gs[0, 3])

    _draw_panel_c((ax_len, ax_scf, ax_n50, ax_gaps), metrics)

    if show_footnote:
        # Footnote about gap measurement
        fig.text(
            0.5, 0.02,
            'Gap counts: N stretches \u2265100 bp per assembly (NC/T2T standard). '
            'Some references (Dog10K_Boxer_Tasha, UU_Cfam_GSD_1.0) use NCBI 100-bp placeholder gaps; '
            'true gap sizes are unavailable for these.',
            ha='center', va='bottom', fontsize=5.5, color='#444',
        )

    return ax_len, ax_scf, ax_n50, ax_gaps


# ============================================================
# Standalone main
# ============================================================
def main():
    fig = plt.figure(figsize=(FIG_W_DOUBLE, 2.5))
    outer_gs = gridspec.GridSpec(1, 1, figure=fig,
                                  left=0.10, right=0.98,
                                  top=0.85, bottom=0.18)
    make_panel_c(fig, outer_gs[0, 0], show_footnote=True)

    out = OUTPUT_DIR / 'fig1' / 'panels' / 'fig2ros_d_v7.pdf'
    save_fig(fig, out)
    plt.close(fig)


if __name__ == '__main__':
    main()
