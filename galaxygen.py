"""
galaxygen.py
============
Core procedural galaxy generator for the Emperum worldbuilding project.

Generates N star systems (nodes) arranged in a 2-D face-on galactic disk with
spiral-arm density waves and an empty forbidden core, then connects them with
L-way hyperspace edges subject to hard geometric constraints.

Constraints enforced
--------------------
1. Forbidden core:  no node inside r < R_CORE; no edge whose segment
   passes through the circle of radius R_CORE centred at (0, 0).
2. Edge length cap: Euclidean distance ≤ L_MAX.
3. Radial gradient: higher node density near the centre (exponential falloff).
4. Spiral arms:     logarithmic density-wave spiral arms with tunable count,
   tightness, thickness, and inter-arm baseline.

Usage (importable)
------------------
    from galaxygen import GalaxyConfig, GalaxyGenerator
    cfg = GalaxyConfig(n_nodes=5000, seed=7)
    gen = GalaxyGenerator(cfg)
    nodes_df, edges_df = gen.run()

Usage (script, uses all defaults)
----------------------------------
    python galaxygen.py
"""

from __future__ import annotations

import dataclasses
import math
import os
import time
from typing import Tuple

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


# ---------------------------------------------------------------------------
# Configuration dataclass
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class GalaxyConfig:
    """All tunable parameters for galaxy generation.

    Spatial units
    -------------
    All distances share an arbitrary unit (think "hundred light-years" or
    whatever you need for worldbuilding).  The generator is scale-agnostic.

    Notes on BOOST
    --------------
    BOOST is an internal multiplier that raises the raw acceptance probability
    P_rad * P_arm into a useful acceptance-rate range.  If you tighten the
    radial or arm preferences substantially you may need a higher BOOST so
    that the rejection sampler converges in a reasonable number of batches.
    """

    # ---- node count ----
    n_nodes: int = 5_000

    # ---- spatial extents ----
    r_disk: float = 100.0   # outer disk boundary radius
    r_core: float = 15.0    # forbidden core radius (no nodes, no edge crossings)

    # ---- edge parameters ----
    l_max: float = 9.0          # maximum edge length (hard limit)
    target_degree: float = 4.0  # desired average node degree (approximate)

    # ---- spiral-arm parameters ----
    n_arms: int = 4             # number of logarithmic spiral arms
    arm_b: float = 0.35         # spiral tightness (larger → more open)
    arm_sigma: float = 5.5      # arm Gaussian half-width (same units as r)
    arm_base: float = 0.15      # baseline acceptance far from any arm ∈ (0, 1)
    r_arm_start: float = 3.0    # spiral r-value at winding angle φ = 0

    # ---- radial density ----
    r_scale: float = 38.0       # exponential scale length; smaller → denser centre

    # ---- rejection sampling tuning ----
    boost: float = 6.0           # acceptance-probability multiplier
    max_sample_iters: int = 30   # maximum batch attempts before giving up

    # ---- reproducibility ----
    seed: int = 7

    # ---- output ----
    out_dir: str = "output"
    write_gexf: bool = True      # attempt GEXF export (requires networkx)


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def segment_min_dist_sq_to_origin(p1s: np.ndarray, p2s: np.ndarray) -> np.ndarray:
    """Minimum squared distance from the origin to each segment p1→p2.

    Handles both a single pair (shape ``(2,)``) and a batch ``(N, 2)``.

    Parameters
    ----------
    p1s, p2s : ndarray, shape ``(2,)`` or ``(N, 2)``

    Returns
    -------
    ndarray of shape ``()`` or ``(N,)`` – squared distances.
    """
    d = p2s - p1s                                  # direction vector(s)
    dot_dd  = np.sum(d  * d,  axis=-1)             # |d|²
    dot_p1d = np.sum(p1s * d, axis=-1)             # p1 · d

    # Clamp t ∈ [0, 1] to stay on the segment; guard against zero-length segs
    denom = np.where(dot_dd > 1e-12, dot_dd, 1e-12)
    t = np.clip(-dot_p1d / denom, 0.0, 1.0)

    if p1s.ndim == 1:
        closest = p1s + t * d
    else:
        closest = p1s + t[:, None] * d

    return np.sum(closest * closest, axis=-1)


def segments_cross_circle(
    p1s: np.ndarray, p2s: np.ndarray, radius: float
) -> np.ndarray:
    """Boolean array: True where segment p1→p2 intersects the circle at the origin.

    Parameters
    ----------
    p1s, p2s : ndarray, shape ``(N, 2)``
    radius   : float

    Returns
    -------
    ndarray of shape ``(N,)``, dtype bool.
    """
    return segment_min_dist_sq_to_origin(p1s, p2s) < radius * radius


# ---------------------------------------------------------------------------
# Spiral-arm helpers
# ---------------------------------------------------------------------------

def arm_distances_batch(
    r: np.ndarray,
    theta: np.ndarray,
    n_arms: int,
    arm_b: float,
    r_arm_start: float,
    r_disk: float,
) -> np.ndarray:
    """Approximate distance of each point (r, θ) to the nearest spiral arm.

    The approximation is ``|r – r_arm(θ)|`` where ``r_arm(θ)`` is the
    logarithmic-spiral radius at the *same* azimuthal angle.  Multiple
    windings are tested and the minimum distance is returned.

    This is intentionally approximate (cheaper than exact shortest Euclidean
    distance to a curve) and sufficient for density-shaping purposes.

    Parameters
    ----------
    r, theta     : 1-D arrays of length N
    n_arms       : number of spiral arms
    arm_b        : logarithmic-spiral tightness parameter
    r_arm_start  : reference radius (arm radius at winding angle φ = 0)
    r_disk       : outer disk radius (used to bound the winding search)

    Returns
    -------
    arm_dist : 1-D array of length N (same units as r)
    """
    min_dist = np.full(len(r), np.inf)

    # How many full windings can fit within [r_arm_start, r_disk]?
    #   r_disk = r_arm_start * exp(b * phi_max)  →  phi_max = log(r_disk/r_arm_start) / b
    if arm_b > 0 and r_arm_start > 0:
        phi_max = math.log(max(r_disk / r_arm_start, 1.01)) / arm_b
        n_windings = int(math.ceil(phi_max / (2.0 * math.pi))) + 1
    else:
        n_windings = 5

    for k in range(n_arms):
        # Angular offset: arm k starts 2πk/N_ARMS ahead of arm 0
        theta_k = theta - 2.0 * math.pi * k / n_arms

        # Sweep winding numbers that can plausibly hit a radius near r
        for n in range(-2, n_windings + 2):
            phi = theta_k + 2.0 * math.pi * n
            r_arm = r_arm_start * np.exp(arm_b * phi)
            # Only compare windings whose radii overlap the disk
            valid = (r_arm > 0.0) & (r_arm <= r_disk * 1.5)
            dist  = np.where(valid, np.abs(r - r_arm), np.inf)
            np.minimum(min_dist, dist, out=min_dist)

    return min_dist


def arm_centerline_points(
    n_arms: int,
    arm_b: float,
    r_arm_start: float,
    r_disk: float,
    n_pts: int = 400,
) -> list[np.ndarray]:
    """Return (x, y) arrays tracing each spiral arm centerline, for plotting.

    Parameters
    ----------
    n_arms, arm_b, r_arm_start, r_disk : matching GalaxyConfig fields
    n_pts : number of sample points per arm

    Returns
    -------
    List of ``n_arms`` arrays, each of shape ``(n_pts, 2)``.
    """
    if arm_b <= 0 or r_arm_start <= 0:
        return []

    # Winding angle needed to reach r_disk from r_arm_start
    phi_max = math.log(max(r_disk / r_arm_start, 1.01)) / arm_b
    phi_vals = np.linspace(0.0, phi_max, n_pts)
    r_vals   = r_arm_start * np.exp(arm_b * phi_vals)

    arms_xy = []
    for k in range(n_arms):
        angle = phi_vals + 2.0 * math.pi * k / n_arms
        x = r_vals * np.cos(angle)
        y = r_vals * np.sin(angle)
        arms_xy.append(np.column_stack([x, y]))

    return arms_xy


# ---------------------------------------------------------------------------
# Main generator class
# ---------------------------------------------------------------------------

class GalaxyGenerator:
    """Procedural generator for a face-on galactic disk with L-way edges.

    Parameters
    ----------
    cfg : GalaxyConfig
        All tunable parameters.  Default values produce a plausible 5 000-node
        galaxy in a few seconds on a standard laptop.
    """

    def __init__(self, cfg: GalaxyConfig) -> None:
        self.cfg = cfg
        self._rng = np.random.default_rng(cfg.seed)

    # ------------------------------------------------------------------
    # Stage A: node-position sampling
    # ------------------------------------------------------------------

    def _sample_nodes(
        self,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Sample N node positions via batch rejection sampling.

        Candidates are drawn uniformly in the annulus [R_CORE, R_DISK] (using
        the area-correct inverse-CDF method), then accepted with probability::

            P = clip(P_rad * P_arm * BOOST, 0, 1)

        where::

            P_rad(r) = exp(-r / R_SCALE)
            P_arm(d) = ARM_BASE + (1 - ARM_BASE) * exp(-(d / ARM_SIGMA)²)
            d        = approximate distance to nearest arm centreline

        Returns
        -------
        xy    : (N, 2) Cartesian positions
        r     : (N,)  radii
        theta : (N,)  azimuths ∈ [0, 2π)
        arm_d : (N,)  approximate arm-centreline distances
        """
        cfg = self.cfg
        collected_xy    = []
        collected_r     = []
        collected_theta = []
        collected_arm_d = []
        total_collected = 0

        # Batch size: at least 4× the target to give the sampler breathing room
        batch_size = max(20_000, cfg.n_nodes * 5)

        for iteration in range(cfg.max_sample_iters):
            # ---- draw candidates ----
            # Uniform-area sampling in annulus: r = sqrt(U * (R²_disk – R²_core) + R²_core)
            r2_lo = cfg.r_core ** 2
            r2_hi = cfg.r_disk ** 2
            r2    = self._rng.uniform(r2_lo, r2_hi, batch_size)
            r     = np.sqrt(r2)
            theta = self._rng.uniform(0.0, 2.0 * math.pi, batch_size)

            # ---- acceptance probabilities ----
            p_rad = np.exp(-r / cfg.r_scale)

            arm_d = arm_distances_batch(
                r, theta,
                cfg.n_arms, cfg.arm_b, cfg.r_arm_start, cfg.r_disk,
            )
            p_arm = cfg.arm_base + (1.0 - cfg.arm_base) * np.exp(
                -(arm_d / cfg.arm_sigma) ** 2
            )

            p = np.clip(p_rad * p_arm * cfg.boost, 0.0, 1.0)

            # ---- accept / reject ----
            accept    = self._rng.random(batch_size) < p
            r_acc     = r[accept]
            theta_acc = theta[accept]
            arm_d_acc = arm_d[accept]

            x = r_acc * np.cos(theta_acc)
            y = r_acc * np.sin(theta_acc)

            collected_xy.append(np.column_stack([x, y]))
            collected_r.append(r_acc)
            collected_theta.append(theta_acc)
            collected_arm_d.append(arm_d_acc)
            total_collected += len(r_acc)

            if total_collected >= cfg.n_nodes:
                break
        else:
            if total_collected < cfg.n_nodes:
                raise RuntimeError(
                    f"Sampled only {total_collected}/{cfg.n_nodes} nodes after "
                    f"{cfg.max_sample_iters} batches.  Try increasing BOOST "
                    f"(currently {cfg.boost}), or loosening radial / arm constraints."
                )

        xy    = np.vstack(collected_xy)   [:cfg.n_nodes]
        r     = np.concatenate(collected_r)    [:cfg.n_nodes]
        theta = np.concatenate(collected_theta)[:cfg.n_nodes]
        arm_d = np.concatenate(collected_arm_d)[:cfg.n_nodes]

        return xy, r, theta, arm_d

    # ------------------------------------------------------------------
    # Stage B: edge generation
    # ------------------------------------------------------------------

    def _build_edges(self, xy: np.ndarray) -> pd.DataFrame:
        """Build L-way edges with KD-tree acceleration and core filtering.

        Algorithm
        ---------
        1. Find all candidate node pairs within L_MAX using ``cKDTree.query_pairs``.
        2. Batch-test every candidate for forbidden-core intersection.
        3. Shuffle the valid candidates, then take the first
           ``ceil(N × TARGET_DEGREE / 2)`` as final undirected edges.

        If fewer valid candidates exist than the target, all valid edges are
        returned and a warning is printed.

        Parameters
        ----------
        xy : (N, 2) array of node Cartesian positions

        Returns
        -------
        DataFrame with columns: source, target, length, weight
        """
        cfg = self.cfg
        n = len(xy)
        target_edges = int(math.ceil(n * cfg.target_degree / 2.0))

        print(f"  Building KD-tree for {n} nodes …")
        tree = cKDTree(xy)

        print(f"  Querying all pairs within L_MAX={cfg.l_max} …")
        # output_type='ndarray' (scipy ≥ 1.6) avoids an expensive set→list→array
        # conversion; fall back gracefully for older scipy versions.
        try:
            pairs = tree.query_pairs(cfg.l_max, output_type="ndarray")
        except TypeError:
            # Older scipy returns a set of (i, j) tuples
            raw = list(tree.query_pairs(cfg.l_max))
            pairs = np.array(raw, dtype=np.int64).reshape(-1, 2)

        if len(pairs) == 0:
            print("  WARNING: no candidate pairs within L_MAX. "
                  "Try increasing --l_max or decreasing --n_nodes.")
            return pd.DataFrame(columns=["source", "target", "length", "weight"])

        print(f"  {len(pairs):,} candidate pairs before core filter …")

        # ---- core-crossing filter (vectorised) ----
        p1s = xy[pairs[:, 0]]
        p2s = xy[pairs[:, 1]]
        crosses    = segments_cross_circle(p1s, p2s, cfg.r_core)
        valid_mask = ~crosses
        valid_pairs = pairs[valid_mask]

        n_removed = int(crosses.sum())
        print(f"  {len(valid_pairs):,} valid pairs after core filter "
              f"({n_removed:,} removed as core-crossing).")

        if len(valid_pairs) == 0:
            print("  WARNING: all candidate pairs cross the core.")
            return pd.DataFrame(columns=["source", "target", "length", "weight"])

        # ---- shuffle and select ----
        if len(valid_pairs) < target_edges:
            print(f"  WARNING: only {len(valid_pairs):,} valid pairs available; "
                  f"target was {target_edges:,}.  "
                  "Try increasing --l_max or --target_degree.")
            selected = valid_pairs
        else:
            perm     = self._rng.permutation(len(valid_pairs))
            selected = valid_pairs[perm[:target_edges]]

        # ---- compute edge attributes ----
        lengths = np.linalg.norm(xy[selected[:, 0]] - xy[selected[:, 1]], axis=1)
        weights = 1.0 / np.maximum(lengths, 1e-9)   # weight = 1/length

        df = pd.DataFrame({
            "source": selected[:, 0].astype(np.int64),
            "target": selected[:, 1].astype(np.int64),
            "length": lengths,
            "weight": weights,
        })
        return df

    # ------------------------------------------------------------------
    # Acceptance tests
    # ------------------------------------------------------------------

    def _run_checks(
        self,
        nodes: pd.DataFrame,
        edges: pd.DataFrame,
        xy: np.ndarray,
    ) -> None:
        """Print acceptance test results to stdout."""
        cfg = self.cfg
        sep = "─" * 52

        print(f"\n{sep}")
        print("  ACCEPTANCE TESTS")
        print(sep)

        # Node count
        ok = len(nodes) == cfg.n_nodes
        print(f"  Node count : {len(nodes):>6,}  (target {cfg.n_nodes:,})  "
              f"{'✓' if ok else '✗ FAIL'}")

        # Radial bounds
        r = nodes["r"].values
        ok_min = r.min() >= cfg.r_core - 1e-6
        ok_max = r.max() <= cfg.r_disk + 1e-6
        print(f"  Min r      : {r.min():>9.4f}  >= {cfg.r_core}  "
              f"{'✓' if ok_min else '✗ FAIL'}")
        print(f"  Max r      : {r.max():>9.4f}  <= {cfg.r_disk}  "
              f"{'✓' if ok_max else '✗ FAIL'}")

        if len(edges) == 0:
            print("  (no edges to check)")
        else:
            # Edge length
            ok_len = edges["length"].max() <= cfg.l_max + 1e-6
            print(f"  Max length : {edges['length'].max():>9.4f}  <= {cfg.l_max}  "
                  f"{'✓' if ok_len else '✗ FAIL'}")

            # Core-crossing check
            n_check = len(edges)
            if n_check <= 50_000:
                sample_edges = edges
                label = f"all {n_check:,}"
            else:
                sample_edges = edges.sample(200, random_state=cfg.seed)
                label = "200 random sample"

            p1s = xy[sample_edges["source"].values]
            p2s = xy[sample_edges["target"].values]
            crosses = segments_cross_circle(p1s, p2s, cfg.r_core)
            ok_core = crosses.sum() == 0
            print(f"  Core cross : {crosses.sum():>3} crossings in {label}  "
                  f"{'✓' if ok_core else '✗ FAIL'}")

            # Degree stats
            src = edges["source"].values
            tgt = edges["target"].values
            deg = np.bincount(
                np.concatenate([src, tgt]), minlength=cfg.n_nodes
            )
            print(f"\n  Degree distribution (target avg ≈ {cfg.target_degree}):")
            print(f"    min={deg.min()}  "
                  f"median={np.median(deg):.1f}  "
                  f"max={deg.max()}  "
                  f"avg={deg.mean():.3f}")
            print(f"    edges={len(edges):,}  "
                  f"(target ≈ {int(cfg.n_nodes * cfg.target_degree / 2):,})")

            # Connected components via path-compressed union-find
            parent = np.arange(cfg.n_nodes, dtype=np.int64)

            def find(x: int) -> int:
                while parent[x] != x:
                    parent[x] = parent[parent[x]]   # path halving
                    x = parent[x]
                return x

            for s, t in zip(src, tgt):
                ps, pt = find(int(s)), find(int(t))
                if ps != pt:
                    parent[ps] = pt

            n_components = len({find(i) for i in range(cfg.n_nodes)})
            print(f"    connected components: {n_components:,}")

        print(sep + "\n")

    # ------------------------------------------------------------------
    # GEXF export
    # ------------------------------------------------------------------

    def _write_gexf(
        self, nodes: pd.DataFrame, edges: pd.DataFrame
    ) -> None:
        """Export a GEXF file for Gephi (requires networkx ≥ 2.0)."""
        try:
            import networkx as nx
        except ImportError:
            print("  networkx not found – skipping GEXF export.  "
                  "Install with: pip install networkx")
            return

        G = nx.Graph()

        # Nodes – cast to plain Python types so networkx serialises cleanly
        for row in nodes.itertuples(index=False):
            G.add_node(
                int(row.id),
                x=float(row.x),
                y=float(row.y),
                r=float(row.r),
                theta=float(row.theta),
                arm_dist=float(row.arm_dist),
            )

        # Edges
        for row in edges.itertuples(index=False):
            G.add_edge(
                int(row.source),
                int(row.target),
                length=float(row.length),
                weight=float(row.weight),
            )

        gexf_path = os.path.join(self.cfg.out_dir, "graph.gexf")
        nx.write_gexf(G, gexf_path)
        print(f"  Wrote {gexf_path}")

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Execute all generation stages; write outputs; return DataFrames.

        Stages
        ------
        A – sample node positions with radial + spiral-arm preferences.
        B – build L-way edges with KD-tree acceleration and core filtering.

        Returns
        -------
        nodes_df : DataFrame  (id, x, y, r, theta, arm_dist)
        edges_df : DataFrame  (source, target, length, weight)
        """
        cfg = self.cfg
        os.makedirs(cfg.out_dir, exist_ok=True)
        t_start = time.perf_counter()

        # ── Stage A ──────────────────────────────────────────────────
        print("Stage A: sampling node positions …")
        t0 = time.perf_counter()
        xy, r, theta, arm_d = self._sample_nodes()
        print(f"  {len(xy):,} nodes sampled in {time.perf_counter() - t0:.2f}s")

        nodes = pd.DataFrame({
            "id":       np.arange(len(xy), dtype=np.int64),
            "x":        xy[:, 0],
            "y":        xy[:, 1],
            "r":        r,
            "theta":    theta,
            "arm_dist": arm_d,
        })

        # ── Stage B ──────────────────────────────────────────────────
        print("\nStage B: building edges …")
        t0 = time.perf_counter()
        edges = self._build_edges(xy)
        print(f"  {len(edges):,} edges built in {time.perf_counter() - t0:.2f}s")

        # ── Acceptance tests ─────────────────────────────────────────
        self._run_checks(nodes, edges, xy)

        # ── Write outputs ─────────────────────────────────────────────
        nodes_path = os.path.join(cfg.out_dir, "nodes.csv")
        edges_path = os.path.join(cfg.out_dir, "edges.csv")
        nodes.to_csv(nodes_path, index=False)
        edges.to_csv(edges_path, index=False)
        print(f"Wrote {nodes_path}")
        print(f"Wrote {edges_path}")

        if cfg.write_gexf:
            self._write_gexf(nodes, edges)

        elapsed = time.perf_counter() - t_start
        print(f"\nTotal time: {elapsed:.2f}s")

        return nodes, edges


# ---------------------------------------------------------------------------
# Script entry point (uses all GalaxyConfig defaults)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    gen = GalaxyGenerator(GalaxyConfig())
    gen.run()
