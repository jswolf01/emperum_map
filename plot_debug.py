"""
plot_debug.py
=============
Matplotlib sanity-check plot for the Emperum galaxy generator.

Shows:
  • Disk boundary circle
  • Forbidden core circle (shaded)
  • Spiral arm centrelines (reference only, not enforced as hard structure)
  • Nodes coloured by arm_dist, radius, pop, admin_lvl, admin_dist, or hierarchy
  • L-way edges (optional; can be slow for N > 1000 – use --no_edges)
  • In hierarchy mode, both nodes and edges are tinted by their admin domain

Usage
-----
    # Default: use ./output/, colour by arm_dist, show edges
    python plot_debug.py

    # Skip edges (much faster at N=5000)
    python plot_debug.py --no_edges

    # Colour nodes by population
    python plot_debug.py --no_edges --color_by pop

    # Colour nodes by admin level
    python plot_debug.py --no_edges --color_by admin_lvl

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

    # Save as SVG (vector, scales to any size)
    python plot_debug.py --no_edges --svg galaxy.svg

    # --svg with no filename defaults to galaxy.svg
    python plot_debug.py --no_edges --svg

    # Point at a different output directory
    python plot_debug.py --out_dir my_run --no_edges

Spiral-arm parameters must match those used during generation so that the
plotted centrelines align with the actual node distribution.  The script
reads --r_disk, --r_core, --n_arms, --arm_b, and --r_arm_start; if you used
non-default values in run_generate.py, pass the same values here.
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Optional

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
    p.add_argument("--svg",      nargs="?", const="galaxy.svg", default=None,
                   metavar="FILE",
                   help="Save figure as SVG (vector format).  "
                        "FILE defaults to 'galaxy.svg' when omitted.  "
                        "Overrides --save when both are given.")

    # Cosmetic toggles
    p.add_argument("--no_edges", action="store_true",
                   help="Skip drawing edges (much faster for large graphs).")
    p.add_argument("--color_by",
                   choices=["arm_dist", "r", "none", "pop", "admin_lvl", "admin_dist",
                            "hierarchy", "is_choke"],
                   default="arm_dist",
                   help="Node colouring scheme.")

    # Node appearance
    p.add_argument("--node_size",  type=float, default=1.5,
                   help="Scatter marker size.")
    p.add_argument("--node_color", default="#aaccff",
                   help="Uniform node colour used when --color_by none.")
    p.add_argument("--gradient_low_color",  default="#ffe8c0",
                   help="Gradient colour at low data values.")
    p.add_argument("--gradient_high_color", default="#0d000f",
                   help="Gradient colour at high data values.")

    # Edge appearance
    p.add_argument("--edge_alpha", type=float, default=0.35,
                   help="Edge line alpha (0=invisible, 1=solid).")
    p.add_argument("--edge_color", default="#2244aa",
                   help="Edge line colour (any matplotlib colour string).")
    p.add_argument("--edge_width", type=float, default=0.4,
                   help="Edge line width in points.")

    # Spatial parameters – auto-loaded from params.json when present;
    # explicit CLI values always take precedence.
    p.add_argument("--r_disk",      type=float, default=None)
    p.add_argument("--r_core",      type=float, default=None)
    p.add_argument("--n_arms",      type=int,   default=None)
    p.add_argument("--arm_b",       type=float, default=None)
    p.add_argument("--r_arm_start", type=float, default=None)

    return p


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def _hierarchy_color(idx: int) -> tuple:
    """Return an RGB triple for hierarchy index *idx*.

    Uses golden-ratio hue spacing in HSV space so that consecutive indices
    remain perceptually distinct even for large numbers of hierarchies.
    """
    golden_ratio = 0.618033988749895
    h = (0.10 + idx * golden_ratio) % 1.0
    s, v = 0.82, 0.88
    hi  = int(h * 6) % 6
    f   = h * 6 - int(h * 6)
    p   = v * (1.0 - s)
    q   = v * (1.0 - s * f)
    t   = v * (1.0 - s * (1.0 - f))
    rgb = [(v, t, p), (q, v, p), (p, v, t), (p, q, v), (t, p, v), (v, p, q)][hi]
    return rgb


def draw_galaxy(args: argparse.Namespace) -> plt.Figure:
    """Load CSV files and draw the galaxy plot.

    Parameters
    ----------
    args : parsed argparse Namespace.  May optionally carry:
        • selected_nodes  – iterable of integer node IDs to highlight
        • found_nodes     – iterable of integer node IDs to highlight (search)

    Returns
    -------
    matplotlib Figure
    """
    # ── Load saved generation parameters ─────────────────────────────────
    _SPATIAL_DEFAULTS = {
        "r_disk": 100.0, "r_core": 15.0, "n_arms": 4,
        "arm_b": 0.35, "r_arm_start": 3.0,
    }
    params_path = os.path.join(args.out_dir, "params.json")
    if os.path.exists(params_path):
        with open(params_path) as _f:
            _saved = json.load(_f)
        for _key, _fallback in _SPATIAL_DEFAULTS.items():
            if getattr(args, _key, None) is None:
                setattr(args, _key, _saved.get(_key, _fallback))
    else:
        for _key, _fallback in _SPATIAL_DEFAULTS.items():
            if getattr(args, _key, None) is None:
                setattr(args, _key, _fallback)

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
    # Equal data aspect keeps the galaxy circular and edges aligned to nodes
    # regardless of canvas size.  "datalim" lets the axes fill available space
    # (no whitespace) and adjusts the shown data range instead of the box size.
    ax.set_aspect("equal", adjustable="datalim")

    BG = "#09090f"
    ax.set_facecolor(BG)
    fig.patch.set_facecolor(BG)

    margin = args.r_disk * 1.08
    ax.set_xlim(-margin, margin)
    ax.set_ylim(-margin, margin)
    # Disable automatic limit expansion so scatter/collections never shift limits.
    ax.autoscale(False)

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

    # ── Hidden-node mask (reduced view / filter) ─────────────────────────
    _hidden_set = set(getattr(args, "hidden_nodes", None) or [])
    if _hidden_set:
        _vis_mask = np.ones(len(nodes), dtype=bool)
        for _hi in _hidden_set:
            if 0 <= _hi < len(nodes):
                _vis_mask[_hi] = False
    else:
        _vis_mask = None   # all visible

    # ── Hierarchy pre-computation ─────────────────────────────────────────
    # Resolved early so it can colour both edges and nodes in hierarchy mode.
    color_by = getattr(args, "color_by", "arm_dist")

    _HIER_UNASSIGNED = (0.06, 0.06, 0.08)   # dark near-black for unassigned nodes
    _HIER_CROSS_EDGE = (0.08, 0.08, 0.12, 0.12)  # very dim for cross-hierarchy edges
    _h_vals: Optional[np.ndarray]   = None   # int64 hierarchy array (length = n nodes)
    _hmap:   dict                   = {}     # hierarchy_id -> RGB triple
    _hier_node_rgba: Optional[np.ndarray] = None  # (n, 4) float32 for scatter

    if color_by == "hierarchy" and "hierarchy" in nodes.columns:
        _h_vals    = nodes["hierarchy"].values.astype(np.int64)
        _unique_h  = sorted(int(v) for v in np.unique(_h_vals) if v >= 0)
        _hmap      = {h: _hierarchy_color(i) for i, h in enumerate(_unique_h)}
        _hier_node_rgba = np.array(
            [(*_hmap.get(int(h), _HIER_UNASSIGNED), 0.85) if h >= 0
             else (*_HIER_UNASSIGNED, 0.30)
             for h in _h_vals],
            dtype=np.float32,
        )

    # ── Node visibility mask (reduced view / filter) ───────────────────────
    # vis[i] = True means node i is rendered; False means hidden.
    _node_mask = getattr(args, "node_mask", None)
    if _node_mask is not None and len(_node_mask) == len(nodes):
        vis = np.asarray(_node_mask, dtype=bool)
    else:
        vis = np.ones(len(nodes), dtype=bool)

    # Fold hidden_nodes mask into vis so a single boolean array drives all rendering.
    if _vis_mask is not None:
        vis &= _vis_mask

    # ── Edges ─────────────────────────────────────────────────────────────
    if not getattr(args, "no_edges", False) and len(edges) > 0:
        xy  = nodes[["x", "y"]].values
        src = edges["source"].values.astype(int)
        tgt = edges["target"].values.astype(int)

        # Filter edges: hide any edge touching a hidden node.
        edge_vis = vis[src] & vis[tgt]
        src_v = src[edge_vis]
        tgt_v = tgt[edge_vis]

        if len(src_v) > 0:
            segs = [[xy[s], xy[t]] for s, t in zip(src_v, tgt_v)]

            if _h_vals is not None:
                # Hierarchy mode: colour edges that share a hierarchy; dim the rest.
                ec = [
                    (*_hmap.get(int(_h_vals[s]), _HIER_UNASSIGNED), 0.55)
                    if (_h_vals[s] >= 0 and int(_h_vals[s]) == int(_h_vals[t]))
                    else _HIER_CROSS_EDGE
                    for s, t in zip(src_v, tgt_v)
                ]
                lc = LineCollection(segs, colors=ec,
                                    linewidths=args.edge_width, zorder=5)
            else:
                lc = LineCollection(
                    segs,
                    colors=args.edge_color,
                    linewidths=args.edge_width,
                    alpha=args.edge_alpha,
                    zorder=5,
                )
            ax.add_collection(lc)

    # ── Nodes ─────────────────────────────────────────────────────────────
    x = nodes["x"].values[vis]
    y = nodes["y"].values[vis]

    # Build user gradient colormap — applied to ALL gradient color modes so
    # the gradient low/high controls in the GUI always take effect.
    grad_cmap = LinearSegmentedColormap.from_list(
        "user_gradient", [args.gradient_low_color, args.gradient_high_color]
    )

    sc     = None
    clabel = None
    cmap   = None

    if _hier_node_rgba is not None:
        # Hierarchy mode: categorical RGBA colours, no gradient colorbar.
        sc = ax.scatter(
            x, y,
            c=_hier_node_rgba[vis],
            s=args.node_size,
            linewidths=0,
            zorder=6,
        )
    else:
        # Colormaps for worldbuilding attributes
        _POP_CMAP   = LinearSegmentedColormap.from_list(
            "pop", ["#0a0a1a", "#1a3a6a", "#2266cc", "#44aaff", "#ffffff"])
        _ADMIN_CMAP = LinearSegmentedColormap.from_list(
            "admin", ["#111111", "#331100", "#884400", "#cc6600", "#ff9900", "#ffff44"])
        _DIST_CMAP  = LinearSegmentedColormap.from_list(
            "adist", ["#ffff44", "#66aa22", "#116633", "#003322", "#000a0a"])

        vmin_val: Optional[float] = None
        vmax_val: Optional[float] = None

        if color_by == "arm_dist" and "arm_dist" in nodes.columns:
            c, cmap, clabel = nodes["arm_dist"].values[vis], grad_cmap, "Arm distance"
        elif color_by == "r" and "r" in nodes.columns:
            c, cmap, clabel = nodes["r"].values[vis], grad_cmap, "Radius"
        elif color_by == "pop" and "pop" in nodes.columns:
            c, cmap, clabel = nodes["pop"].values[vis], _POP_CMAP, "Population"
            # Anchor to the full semantic range so every colour level is correct
            # regardless of whether all values 0-100 are present in this dataset.
            vmin_val, vmax_val = 0.0, 100.0
        elif color_by == "admin_lvl" and "admin_lvl" in nodes.columns:
            c, cmap, clabel = nodes["admin_lvl"].values[vis], _ADMIN_CMAP, "Admin level"
            # Anchor: 0 = no admin centre (darkest), 5 = top-tier capital (brightest).
            # Without this, matplotlib auto-normalises to whichever levels happen to
            # exist and a level-4 node would incorrectly get the level-5 colour.
            vmin_val, vmax_val = 0.0, 5.0
        elif color_by == "admin_dist" and "admin_dist" in nodes.columns:
            raw = nodes["admin_dist"].values.astype(float)
            valid_vals = raw[raw >= 0]
            raw_max = float(valid_vals.max()) + 1 if len(valid_vals) > 0 else 1.0
            raw[raw < 0] = raw_max   # map unreachable (-1) to just beyond the max
            c, cmap, clabel = raw[vis], _DIST_CMAP, "Admin distance (hops)"
            vmin_val, vmax_val = 0.0, raw_max
        elif color_by == "is_choke" and "is_choke" in nodes.columns:
            c = nodes["is_choke"].values[vis].astype(float)
            cmap, clabel = "Reds", "Chokepoint"
            vmin_val, vmax_val = 0.0, 1.0
        else:
            c = args.node_color
            cmap = None

        sc = ax.scatter(
            x, y,
            c=c, cmap=cmap,
            vmin=vmin_val, vmax=vmax_val,
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
    # Re-lock limits after colorbar; colorbar resizes ax but must not shift data coords.
    ax.set_xlim(-margin, margin)
    ax.set_ylim(-margin, margin)

    # ── Search / selection highlights ─────────────────────────────────────
    found_nodes = getattr(args, "found_nodes", None)
    if found_nodes is not None and len(found_nodes) > 0:
        found_vis = [fn for fn in found_nodes if fn < len(nodes) and vis[fn]]
        if found_vis:
            fx = nodes["x"].values[found_vis]
            fy = nodes["y"].values[found_vis]
            ax.scatter(fx, fy, s=max(args.node_size * 8, 30),
                       facecolors="none", edgecolors="#00ff88",
                       linewidths=1.5, zorder=9, label="Search results")

    selected_nodes = getattr(args, "selected_nodes", None)
    if selected_nodes is not None and len(selected_nodes) > 0:
        sel_vis = [sn for sn in selected_nodes if sn < len(nodes) and vis[sn]]
        if sel_vis:
            sx = nodes["x"].values[sel_vis]
            sy = nodes["y"].values[sel_vis]
            ax.scatter(sx, sy, s=max(args.node_size * 12, 50),
                       facecolors="none", edgecolors="#ffff00",
                       linewidths=2.0, zorder=10, label="Selected")

    # ── Decorations ───────────────────────────────────────────────────────
    no_edges = getattr(args, "no_edges", False)
    n_visible = int(vis.sum())
    n_hidden  = len(nodes) - n_visible
    _filter_note = (
        f"  [{n_visible:,} of {len(nodes):,} shown]"
        if n_hidden > 0 else ""
    )
    title = (
        f"Emperum Galaxy  —  "
        f"{len(nodes):,} systems  |  {len(edges):,} L-ways"
        + ("  (edges hidden)" if no_edges else "")
        + _filter_note
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
    if not no_edges:
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

    return fig


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = build_parser()
    args   = parser.parse_args()

    fig = draw_galaxy(args)

    if args.svg:
        fig.savefig(args.svg, format="svg", bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        print(f"Saved figure to {args.svg}")
    elif args.save:
        fig.savefig(args.save, dpi=150, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        print(f"Saved figure to {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
