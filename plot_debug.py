"""
plot_debug.py
=============
Matplotlib sanity-check plot for the Emperum galaxy generator.

Shows:
  • Disk boundary circle
  • Forbidden core circle (shaded)
  • Spiral arm centrelines (reference only, not enforced as hard structure)
  • Nodes coloured by arm_dist or radius (toggle with --color_by)
  • L-way edges (optional; can be slow for N > 1000 – use --no_edges)

Usage
-----
    # Default: use ./output/, colour by arm_dist, show edges
    python plot_debug.py

    # Skip edges (much faster at N=5000)
    python plot_debug.py --no_edges

    # Colour nodes by radius instead
    python plot_debug.py --no_edges --color_by r

    # Uniform node colour (no gradient)
    python plot_debug.py --color_by none --node_color "#ffffff"

    # Custom gradient: arm stars bright cyan, inter-arm very dark
    python plot_debug.py --gradient_low_color "#00ffff" --gradient_high_color "#050005"

    # Custom edge appearance
    python plot_debug.py --edge_color "#44aaff" --edge_width 0.8 --edge_alpha 0.5

    # Save to PNG instead of opening an interactive window
    python plot_debug.py --no_edges --save galaxy.png

    # Point at a different output directory
    python plot_debug.py --out_dir my_run --no_edges

Spiral-arm parameters must match those used during generation so that the
plotted centrelines align with the actual node distribution.  The script
reads --r_disk, --r_core, --n_arms, --arm_b, and --r_arm_start; if you used
non-default values in run_generate.py, pass the same values here.
"""

from __future__ import annotations

import argparse
import os

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection

from galaxygen import arm_centerline_points


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plot_debug.py",
        description="Debug visualisation for the Emperum galaxy generator.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Output location
    p.add_argument("--out_dir",  default="output",
                   help="Directory containing nodes.csv and edges.csv.")
    p.add_argument("--save",     default=None, metavar="FILE",
                   help="Save figure to FILE (png/pdf/svg) instead of displaying.")

    # Cosmetic toggles
    p.add_argument("--no_edges", action="store_true",
                   help="Skip drawing edges (much faster for large graphs).")
    p.add_argument("--color_by", choices=["arm_dist", "r", "none"],
                   default="arm_dist",
                   help="Node colouring scheme: arm_dist (gradient by arm proximity), "
                        "r (gradient by radius), or none (uniform colour).")

    # Node appearance
    p.add_argument("--node_size",  type=float, default=1.5,
                   help="Scatter marker size.")
    p.add_argument("--node_color", default="#aaccff",
                   help="Uniform node colour used when --color_by none.")
    p.add_argument("--gradient_low_color",  default="#ffe8c0",
                   help="Gradient colour at low data values "
                        "(close-to-arm for arm_dist; inner disk for r).")
    p.add_argument("--gradient_high_color", default="#0d000f",
                   help="Gradient colour at high data values "
                        "(far-from-arm for arm_dist; outer disk for r).")

    # Edge appearance
    p.add_argument("--edge_alpha", type=float, default=0.35,
                   help="Edge line alpha (0=invisible, 1=solid).")
    p.add_argument("--edge_color", default="#2244aa",
                   help="Edge line colour (any matplotlib colour string).")
    p.add_argument("--edge_width", type=float, default=0.4,
                   help="Edge line width in points.")

    # Spatial parameters – must match the generation run
    p.add_argument("--r_disk",      type=float, default=100.0)
    p.add_argument("--r_core",      type=float, default=15.0)
    p.add_argument("--n_arms",      type=int,   default=4)
    p.add_argument("--arm_b",       type=float, default=0.35)
    p.add_argument("--r_arm_start", type=float, default=3.0)

    return p


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def draw_galaxy(args: argparse.Namespace) -> plt.Figure:
    """Load CSV files and draw the galaxy plot.

    Parameters
    ----------
    args : parsed argparse Namespace

    Returns
    -------
    matplotlib Figure
    """
    nodes_path = os.path.join(args.out_dir, "nodes.csv")
    edges_path = os.path.join(args.out_dir, "edges.csv")

    if not os.path.exists(nodes_path):
        raise FileNotFoundError(
            f"nodes.csv not found in '{args.out_dir}'.  "
            "Run run_generate.py first."
        )

    nodes = pd.read_csv(nodes_path)
    edges = pd.read_csv(edges_path) if os.path.exists(edges_path) else pd.DataFrame()

    # ── Figure setup ─────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect("equal")

    BG = "#09090f"
    ax.set_facecolor(BG)
    fig.patch.set_facecolor(BG)

    margin = args.r_disk * 1.08
    ax.set_xlim(-margin, margin)
    ax.set_ylim(-margin, margin)

    # ── Disk boundary ─────────────────────────────────────────────────────
    disk_circle = plt.Circle(
        (0, 0), args.r_disk,
        fill=False, edgecolor="#3a3a5c", linewidth=1.0, linestyle="--", zorder=2,
    )
    ax.add_patch(disk_circle)

    # ── Forbidden core ────────────────────────────────────────────────────
    core_circle = plt.Circle(
        (0, 0), args.r_core,
        facecolor="#1a0505", edgecolor="#cc3333",
        linewidth=1.5, linestyle="-", zorder=3,
    )
    ax.add_patch(core_circle)

    # ── Spiral arm centrelines ────────────────────────────────────────────
    arm_lines = arm_centerline_points(
        args.n_arms, args.arm_b, args.r_arm_start, args.r_disk
    )
    for arm_xy in arm_lines:
        ax.plot(
            arm_xy[:, 0], arm_xy[:, 1],
            color="#2a3a4a", linewidth=1.0, alpha=0.6, zorder=4,
        )

    # ── Edges ─────────────────────────────────────────────────────────────
    if not args.no_edges and len(edges) > 0:
        xy  = nodes[["x", "y"]].values
        src = edges["source"].values.astype(int)
        tgt = edges["target"].values.astype(int)
        segs = [[xy[s], xy[t]] for s, t in zip(src, tgt)]
        lc = LineCollection(
            segs,
            colors=args.edge_color,
            linewidths=args.edge_width,
            alpha=args.edge_alpha,
            zorder=5,
        )
        ax.add_collection(lc)

    # ── Nodes ─────────────────────────────────────────────────────────────
    x = nodes["x"].values
    y = nodes["y"].values

    # Build gradient colormap from the two user-supplied endpoint colours.
    # Default: low values (close-to-arm) → warm bright, high → very dark,
    # so spiral-arm stars stand out prominently against the background.
    grad_cmap = LinearSegmentedColormap.from_list(
        "user_gradient", [args.gradient_low_color, args.gradient_high_color]
    )

    if args.color_by == "arm_dist" and "arm_dist" in nodes.columns:
        c      = nodes["arm_dist"].values
        cmap   = grad_cmap
        clabel = "Arm distance"
    elif args.color_by == "r" and "r" in nodes.columns:
        c      = nodes["r"].values
        cmap   = grad_cmap
        clabel = "Radius"
    else:
        c      = args.node_color
        cmap   = None
        clabel = None

    sc = ax.scatter(
        x, y,
        c=c, cmap=cmap,
        s=args.node_size,
        alpha=0.75,
        linewidths=0,
        zorder=6,
    )

    if clabel and cmap:
        cbar = plt.colorbar(sc, ax=ax, pad=0.01, fraction=0.03, shrink=0.85)
        cbar.set_label(clabel, color="white", fontsize=9)
        cbar.ax.yaxis.set_tick_params(color="white", labelsize=7)
        plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color="white")

    # ── Decorations ───────────────────────────────────────────────────────
    title = (
        f"Emperum Galaxy  —  "
        f"{len(nodes):,} systems  |  {len(edges):,} L-ways"
        + ("  (edges hidden)" if args.no_edges else "")
    )
    ax.set_title(title, color="white", fontsize=11, pad=10)

    for spine in ax.spines.values():
        spine.set_edgecolor("#2a2a3a")
    ax.tick_params(colors="#555566", labelsize=7)

    legend_patches = [
        mpatches.Patch(facecolor="#3a3a5c", edgecolor="#3a3a5c",
                       label=f"Disk boundary (r={args.r_disk})"),
        mpatches.Patch(facecolor="#1a0505", edgecolor="#cc3333",
                       label=f"Forbidden core (r={args.r_core})"),
        mpatches.Patch(facecolor="#2a3a4a", label="Arm centrelines"),
    ]
    if not args.no_edges:
        legend_patches.append(
            mpatches.Patch(facecolor=args.edge_color, label="L-way edges")
        )

    ax.legend(
        handles=legend_patches,
        loc="upper right",
        fontsize=8,
        facecolor="#111122",
        edgecolor="#333355",
        labelcolor="white",
    )

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = build_parser()
    args   = parser.parse_args()

    fig = draw_galaxy(args)

    if args.save:
        fig.savefig(args.save, dpi=150, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        print(f"Saved figure to {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
