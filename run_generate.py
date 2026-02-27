"""
run_generate.py
===============
CLI entrypoint for the Emperum procedural galaxy generator.

All parameters are optional; unspecified parameters fall back to the defaults
defined in ``GalaxyConfig``.

Quick start
-----------
    python run_generate.py

With custom parameters (matching the default preset)::

    python run_generate.py \\
        --n_nodes 5000 \\
        --r_disk 100 \\
        --r_core 15 \\
        --l_max 9 \\
        --target_degree 4 \\
        --n_arms 4 \\
        --arm_sigma 5.5 \\
        --arm_b 0.35 \\
        --r_scale 38 \\
        --seed 7 \\
        --out_dir output

Then visualise the result::

    python plot_debug.py
"""

import argparse
import json
import os
import sys

from galaxygen import GalaxyConfig, GalaxyGenerator


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_generate.py",
        description=(
            "Procedural galaxy generator for the Emperum worldbuilding project.\n"
            "Produces nodes.csv, edges.csv, and (optionally) graph.gexf in OUT_DIR."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ── Spatial parameters ────────────────────────────────────────────────
    p.add_argument(
        "--n_nodes", type=int, default=5_000,
        metavar="N",
        help="Number of star-system nodes to generate.",
    )
    p.add_argument(
        "--r_disk", type=float, default=100.0,
        metavar="R",
        help="Outer disk boundary radius.",
    )
    p.add_argument(
        "--r_core", type=float, default=15.0,
        metavar="R",
        help="Forbidden core radius (no nodes; no edges may cross this circle).",
    )

    # ── Edge parameters ───────────────────────────────────────────────────
    p.add_argument(
        "--l_max", type=float, default=9.0,
        metavar="L",
        help="Maximum edge (L-way) length.",
    )
    p.add_argument(
        "--target_degree", type=float, default=4.0,
        metavar="D",
        help="Target average node degree (approximate).",
    )
    p.add_argument(
        "--node_connection_chance", type=float, default=1.0,
        metavar="P",
        help=(
            "Probability [0, 1] that any given node is eligible to have edge "
            "connections.  Nodes that are not selected remain fully isolated.  "
            "1.0 = all nodes may connect (default); 0.3 = ~30%% of nodes connect."
        ),
    )

    # ── Spiral-arm parameters ─────────────────────────────────────────────
    p.add_argument(
        "--n_arms", type=int, default=4,
        metavar="N",
        help="Number of logarithmic spiral arms.",
    )
    p.add_argument(
        "--arm_b", type=float, default=0.35,
        metavar="B",
        help="Spiral tightness parameter (larger = more open / fewer windings).",
    )
    p.add_argument(
        "--arm_sigma", type=float, default=5.5,
        metavar="S",
        help="Arm Gaussian half-width (same units as radii).",
    )
    p.add_argument(
        "--arm_base", type=float, default=0.15,
        metavar="P",
        help="Baseline acceptance probability far from any arm (0 < arm_base < 1).",
    )

    # ── Radial density ────────────────────────────────────────────────────
    p.add_argument(
        "--r_scale", type=float, default=38.0,
        metavar="R",
        help="Exponential scale length for radial density falloff.",
    )

    # ── Sampling tuning ───────────────────────────────────────────────────
    p.add_argument(
        "--boost", type=float, default=6.0,
        metavar="B",
        help="Acceptance-probability multiplier (increase if sampling is slow).",
    )

    # ── Chokepoints ───────────────────────────────────────────────────────
    p.add_argument(
        "--choke_count", type=int, default=0,
        metavar="N",
        help=(
            "Chokepoint count.  0 = none (default); -1 = random-chance mode "
            "(see --choke_chance); N>0 = exactly N chokepoints."
        ),
    )
    p.add_argument(
        "--choke_chance", type=float, default=0.3,
        metavar="P",
        help=(
            "Probability [0, 1] that any chokepoints spawn when "
            "--choke_count=-1.  Ignored otherwise."
        ),
    )
    p.add_argument(
        "--choke_sep", type=int, default=-1,
        metavar="H",
        help=(
            "Minimum hop-distance between chokepoints.  "
            "-1 = auto (max(3, diameter // 5))."
        ),
    )

    # ── Reproducibility ───────────────────────────────────────────────────
    p.add_argument(
        "--seed", type=int, default=7,
        metavar="S",
        help="Random seed for reproducible output.",
    )

    # ── Output ───────────────────────────────────────────────────────────
    p.add_argument(
        "--out_dir", type=str, default="output",
        metavar="DIR",
        help="Directory to write output files (created if absent).",
    )
    p.add_argument(
        "--no_gexf", action="store_true",
        help="Skip GEXF export (useful when networkx is not installed).",
    )

    return p


def main() -> None:
    parser = build_parser()
    args   = parser.parse_args()

    cfg = GalaxyConfig(
        n_nodes                = args.n_nodes,
        r_disk                 = args.r_disk,
        r_core                 = args.r_core,
        l_max                  = args.l_max,
        target_degree          = args.target_degree,
        node_connection_chance = args.node_connection_chance,
        n_arms                 = args.n_arms,
        arm_b                  = args.arm_b,
        arm_sigma              = args.arm_sigma,
        arm_base               = args.arm_base,
        r_scale                = args.r_scale,
        boost                  = args.boost,
        seed                   = args.seed,
        out_dir                = args.out_dir,
        write_gexf             = not args.no_gexf,
        choke_count            = args.choke_count,
        choke_chance           = args.choke_chance,
        choke_sep              = args.choke_sep,
    )

    # Print config so the user can confirm parameters before waiting
    print("Configuration")
    print("─" * 40)
    for field in cfg.__dataclass_fields__:
        print(f"  {field:<22} = {getattr(cfg, field)}")
    print()

    gen = GalaxyGenerator(cfg)
    nodes_df, edges_df = gen.run()

    # Persist generation parameters so plot_debug.py can read them automatically
    params_path = os.path.join(cfg.out_dir, "params.json")
    params = {field: getattr(cfg, field) for field in cfg.__dataclass_fields__}
    with open(params_path, "w") as f:
        json.dump(params, f, indent=2)
    print(f"Wrote {params_path}")

    # Remind user of next steps
    print(
        f"\nNext steps:\n"
        f"  • Debug plot : python plot_debug.py --out_dir {cfg.out_dir}\n"
        f"  • Gephi      : import {cfg.out_dir}/graph.gexf  "
        f"(or see README.md for CSV import steps)"
    )


if __name__ == "__main__":
    main()
