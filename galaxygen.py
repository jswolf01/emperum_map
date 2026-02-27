"""
galaxygen.py
============
Core procedural galaxy generator for the Emperum worldbuilding project.

Generates N star systems (nodes) arranged in a 2-D face-on galactic disk with
spiral-arm density waves and an empty forbidden core, then connects them with
L-way hyperspace edges subject to hard geometric constraints.

After spatial generation, worldbuilding attributes are assigned to every node:
  • uid        – unique 7-digit system registry number (zero-padded string)
  • name       – free-text name (blank by default)
  • pop        – population index 0-100 (0 = uninhabited / isolated)
  • admin_lvl  – administrative hierarchy level 0-5 (0 = none)
  • admin_dist – minimum hop-count to the nearest admin node (−1 if none exist)
  • hierarchy  – integer index of the admin-level-1 domain this node belongs to
                 (−1 if unpopulated or unreachable from any lvl-1 root)

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
from collections import deque
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
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
    node_connection_chance: float = 1.0  # probability [0, 1] that a node is eligible for any edges

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

    # ---- population parameters ----
    # pop_core_dispersal: 0 = no radial gradient (uniform distribution across
    #   disk); 1 = moderate core-concentration (inner systems tend higher pop);
    #   values > 1 = stronger core bias.
    pop_core_dispersal: float = 1.0
    # pop_dispersal: spatial clustering strength.  0 = independent assignment;
    #   1 = moderate neighbour-averaging; values > 1 = stronger clustering.
    pop_dispersal: float = 1.0

    # ---- admin hierarchy – exact counts (−1 = compute from ratios) ----
    admin_count_1: int = -1   # number of lvl-1 (top) administrative centres
    admin_count_2: int = -1
    admin_count_3: int = -1
    admin_count_4: int = -1
    admin_count_5: int = -1   # number of lvl-5 (smallest regional) centres

    # ---- admin hierarchy – ratios between successive levels ----
    # admin_ratio_XY = (# lvl-X nodes) / (# lvl-Y nodes)  where X > Y
    admin_ratio_21: float = 2.0   # lvl-2 : lvl-1
    admin_ratio_32: float = 3.0   # lvl-3 : lvl-2
    admin_ratio_43: float = 4.0   # lvl-4 : lvl-3
    admin_ratio_54: float = 5.0   # lvl-5 : lvl-4

    # ---- admin hierarchy – spatial spacing ----
    # admin_coverage_scale: multiplier on the auto-computed coverage radius.
    #   1.0 = each centre covers ~n_connected/n_L nodes (log-distance rule).
    #   < 1 = denser placement, smaller regions; > 1 = sparser, larger regions.
    admin_coverage_scale: float = 1.0
    # admin_sep_scale: sep_radius = max(1, round(coverage_radius × sep_scale)).
    #   Minimum hop-distance enforced between two centres of the SAME level.
    #   1.0 = centres just touching coverage zones; 1.5 = moderate overlap buffer;
    #   2.0+ = strong separation, fewer centres will fit.
    admin_sep_scale: float = 1.5

    # ---- chokepoints ----
    # Chokepoints are nodes with unusually high betweenness centrality —
    # they provide exclusive access to large sectors of the galaxy.
    #
    # choke_count:  0  = disabled (default, no chokepoints)
    #              -1  = random-chance mode: roll choke_chance to decide if
    #                    any chokepoints spawn, then pick a random count
    #             N>0  = designate exactly N chokepoints
    choke_count: int = 0
    # choke_chance: probability [0, 1] that chokepoints appear when
    #   choke_count == -1.  Ignored for other choke_count values.
    choke_chance: float = 0.3
    # choke_sep: minimum hop-distance enforced between any two chokepoints.
    #   -1 = auto (derived from graph diameter: max(3, diameter // 5))
    choke_sep: int = -1


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
# Standalone attribute helpers (callable from GUI without full regeneration)
# ---------------------------------------------------------------------------

def compute_degree(n_nodes: int, edges: pd.DataFrame) -> np.ndarray:
    """Return per-node degree array of length *n_nodes*."""
    degree = np.zeros(n_nodes, dtype=np.int64)
    if len(edges) > 0:
        src = edges["source"].values.astype(np.int64)
        tgt = edges["target"].values.astype(np.int64)
        np.add.at(degree, src, 1)
        np.add.at(degree, tgt, 1)
    return degree


def compute_hierarchy(nodes: pd.DataFrame, edges: pd.DataFrame) -> np.ndarray:
    """Assign each populated node to a discrete administrative hierarchy.

    The hierarchy chain runs from a node to the nearest node at the lowest
    extant admin level (e.g. admin_lvl 5), then hop-by-hop up through each
    extant level until reaching an admin-level-1 root.  "Nearest" means
    fewest graph hops, not Euclidean distance.

    Each admin-level-1 node is the root of one discrete hierarchy.  Every
    populated non-admin node that can reach a chain back to an admin-level-1
    node is assigned the integer index of that root (0-based, ordered by
    ascending node ID of lvl-1 nodes).

    Parameters
    ----------
    nodes : DataFrame  – must contain 'admin_lvl' and 'pop' columns
    edges : DataFrame  – must contain 'source' and 'target' columns

    Returns
    -------
    hierarchy : int64 array of length len(nodes)
        Non-negative integer = hierarchy index shared by all nodes in the
        same administrative domain.  −1 = unassigned (unpopulated or
        unreachable from any admin-level-1 root through the chain).
    """
    n = len(nodes)
    hierarchy = np.full(n, -1, dtype=np.int64)

    admin_lvl_arr = nodes["admin_lvl"].values
    pop_arr       = nodes["pop"].values

    # Determine which admin levels actually have at least one node.
    extant = sorted(lvl for lvl in range(1, 6) if (admin_lvl_arr == lvl).any())

    # Hierarchy requires at least one admin-level-1 root node.
    if not extant or extant[0] != 1:
        return hierarchy

    # Build adjacency list.
    adj: list[list[int]] = [[] for _ in range(n)]
    if len(edges) > 0:
        src = edges["source"].values.astype(int)
        tgt = edges["target"].values.astype(int)
        for s, t in zip(src, tgt):
            adj[s].append(t)
            adj[t].append(s)

    # Node-index lists per extant admin level.
    level_nodes: dict[int, np.ndarray] = {
        lvl: np.where(admin_lvl_arr == lvl)[0] for lvl in extant
    }

    # Step 1 – each admin-level-1 node seeds its own hierarchy ID.
    for hier_id, node_id in enumerate(level_nodes[1]):
        hierarchy[int(node_id)] = hier_id

    # Step 2 – propagate hierarchy downward through consecutive extant levels.
    # For the transition parent_lvl → child_lvl, BFS outward from every
    # parent-level node that already has a hierarchy assignment and stamp
    # the first-contact child-level nodes with that hierarchy.
    for i in range(1, len(extant)):
        parent_lvl = extant[i - 1]
        child_lvl  = extant[i]
        child_set  = {int(x) for x in level_nodes[child_lvl]}

        visited: dict[int, int] = {}
        q: deque[int] = deque()
        for p in level_nodes[parent_lvl]:
            p = int(p)
            if hierarchy[p] >= 0:
                visited[p] = int(hierarchy[p])
                q.append(p)

        while q:
            node = q.popleft()
            for nb in adj[node]:
                if nb not in visited:
                    visited[nb] = visited[node]
                    if nb in child_set:
                        hierarchy[nb] = visited[nb]
                    q.append(nb)

    # Step 3 – assign populated non-admin nodes to the nearest node at the
    # lowest extant admin level (e.g. admin_lvl 5 if present).
    lowest_lvl = extant[-1]
    is_target  = (pop_arr > 0) & (admin_lvl_arr == 0)

    visited3: dict[int, int] = {}
    q3: deque[int] = deque()
    for ln in level_nodes[lowest_lvl]:
        ln = int(ln)
        if hierarchy[ln] >= 0:
            visited3[ln] = int(hierarchy[ln])
            q3.append(ln)

    while q3:
        node = q3.popleft()
        for nb in adj[node]:
            if nb not in visited3:
                visited3[nb] = visited3[node]
                if is_target[nb]:
                    hierarchy[nb] = visited3[nb]
                q3.append(nb)

    return hierarchy


def compute_admin_dist(nodes: pd.DataFrame, edges: pd.DataFrame) -> np.ndarray:
    """Compute minimum hop-distance to nearest admin node for every node.

    Uses multi-source BFS starting from all nodes where admin_lvl > 0.
    Nodes unreachable from any admin node receive distance −1.

    Parameters
    ----------
    nodes : DataFrame  – must contain 'admin_lvl' column
    edges : DataFrame  – must contain 'source' and 'target' columns

    Returns
    -------
    dist : int64 array of length len(nodes)
    """
    n = len(nodes)
    admin_mask = nodes["admin_lvl"].values > 0
    admin_indices = np.where(admin_mask)[0]

    dist = np.full(n, -1, dtype=np.int64)

    if len(admin_indices) == 0:
        return dist

    # Build adjacency list
    adj: list[list[int]] = [[] for _ in range(n)]
    if len(edges) > 0:
        src = edges["source"].values.astype(int)
        tgt = edges["target"].values.astype(int)
        for s, t in zip(src, tgt):
            adj[s].append(t)
            adj[t].append(s)

    # Multi-source BFS
    queue: deque[int] = deque()
    for idx in admin_indices:
        dist[idx] = 0
        queue.append(int(idx))

    while queue:
        node = queue.popleft()
        for neighbor in adj[node]:
            if dist[neighbor] == -1:
                dist[neighbor] = dist[node] + 1
                queue.append(neighbor)

    return dist


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

        # ---- node-eligibility filter ----
        # Each node is independently eligible to participate in edges with
        # probability node_connection_chance.  Pairs where either endpoint is
        # ineligible are dropped so that ineligible nodes remain isolated.
        if cfg.node_connection_chance < 1.0:
            connectable = self._rng.random(n) < cfg.node_connection_chance
            both_ok = connectable[valid_pairs[:, 0]] & connectable[valid_pairs[:, 1]]
            valid_pairs = valid_pairs[both_ok]
            n_eligible = int(connectable.sum())
            print(f"  {n_eligible:,} / {n} nodes eligible "
                  f"({cfg.node_connection_chance:.1%} chance); "
                  f"{len(valid_pairs):,} pairs remain after eligibility filter.")
            if len(valid_pairs) == 0:
                print("  WARNING: no valid pairs after eligibility filter.")
                return pd.DataFrame(columns=["source", "target", "length", "weight"])

        # ---- shuffle then spanning-forest-first selection ----
        # Shuffling first preserves randomness.  Then a union-find spanning
        # forest ensures all nodes that have at least one valid pair end up in
        # a single connected component (no isolated islands).  Extra edges fill
        # up to target_edges after the spanning forest is built.
        perm     = self._rng.permutation(len(valid_pairs))
        shuffled = valid_pairs[perm]

        parent = np.arange(n, dtype=np.int64)

        def _find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]   # path halving
                x = parent[x]
            return x

        span_idx: list[int] = []
        extra_idx: list[int] = []
        for i in range(len(shuffled)):
            s, t = int(shuffled[i, 0]), int(shuffled[i, 1])
            ps, pt = _find(s), _find(t)
            if ps != pt:
                parent[ps] = pt
                span_idx.append(i)
            else:
                extra_idx.append(i)

        n_span   = len(span_idx)
        n_extra  = max(0, target_edges - n_span)
        all_idx  = span_idx + extra_idx[:n_extra]

        if not all_idx:
            print("  WARNING: no valid pairs could be selected.")
            return pd.DataFrame(columns=["source", "target", "length", "weight"])

        if len(valid_pairs) < target_edges:
            print(f"  WARNING: only {len(valid_pairs):,} valid pairs available; "
                  f"target was {target_edges:,}.  "
                  "Try increasing --l_max or --target_degree.")

        selected = shuffled[np.array(all_idx, dtype=np.int64)]
        n_components = len({_find(int(i)) for i in range(n)})
        print(f"  Connected components after spanning forest: {n_components:,} "
              f"(spanning forest used {n_span:,} edges)")

        # ---- bridge isolated components into a single connected graph ----
        # The spanning forest guarantees that every node with at least one
        # valid pair is reachable from every other such node *within its
        # component*.  But low node counts, low connection-chance, or low
        # l_max can leave multiple disconnected components.  We add the
        # minimum set of "bridge" edges required to merge them all.
        bridge_pairs: list[tuple[int, int]] = []
        if n_components > 1:
            bridge_pairs = self._bridge_components(xy, parent, n, tree)
            n_after = len({_find(int(i)) for i in range(n)})
            print(f"  After bridging: {n_after:,} component(s) "
                  f"({len(bridge_pairs):,} bridge edge(s) added).")

        # Combine spanning-forest + extra + bridge edges
        all_pairs = list(selected)
        for bp in bridge_pairs:
            all_pairs.append(np.array(bp, dtype=np.int64))

        if not all_pairs:
            print("  WARNING: no valid pairs could be selected.")
            return pd.DataFrame(columns=["source", "target", "length", "weight"])

        selected_final = np.array(all_pairs, dtype=np.int64)
        # ---- ensure connectivity ----
        # After random selection, connected nodes may form multiple isolated
        # islands.  We add the minimum-length bridge edges (drawn from the
        # valid-pairs pool) needed to merge every component into one, using a
        # Kruskal-style scan over length-sorted valid pairs.
        if len(selected) > 0:
            _par: list = list(range(n))

            def _find(x: int) -> int:
                root = x
                while _par[root] != root:
                    root = _par[root]
                while _par[x] != root:          # path compression
                    _par[x], x = root, _par[x]
                return root

            def _union(a: int, b: int) -> bool:
                ra, rb = _find(a), _find(b)
                if ra == rb:
                    return False
                _par[ra] = rb
                return True

            for _s, _t in selected:
                _union(int(_s), int(_t))

            _touched = set(int(v) for v in selected.reshape(-1))
            _n_comp = len({_find(i) for i in _touched})

            if _n_comp > 1:
                # Build set of already-selected pairs for O(1) skip.
                _sel_set = {(min(int(_s), int(_t)), max(int(_s), int(_t)))
                            for _s, _t in selected}
                # Sort all valid pairs by Euclidean length (ascending).
                _all_lens = np.linalg.norm(
                    xy[valid_pairs[:, 0]] - xy[valid_pairs[:, 1]], axis=1)
                _order = np.argsort(_all_lens)
                _bridges: list = []
                for _idx in _order:
                    _s = int(valid_pairs[_idx, 0])
                    _t = int(valid_pairs[_idx, 1])
                    _key = (min(_s, _t), max(_s, _t))
                    if _key in _sel_set:
                        continue
                    if _union(_s, _t):
                        _bridges.append([_s, _t])
                        _sel_set.add(_key)
                        _n_comp -= 1
                        if _n_comp <= 1:
                            break
                if _bridges:
                    _bridge_arr = np.array(_bridges, dtype=np.int64)
                    selected = np.concatenate([selected, _bridge_arr], axis=0)
                    print(f"  Added {len(_bridges):,} bridge edge(s) to merge "
                          f"isolated network islands.")

        # ---- compute edge attributes ----
        lengths = np.linalg.norm(
            xy[selected_final[:, 0]] - xy[selected_final[:, 1]], axis=1
        )
        weights = 1.0 / np.maximum(lengths, 1e-9)   # weight = 1/length

        df = pd.DataFrame({
            "source": selected_final[:, 0].astype(np.int64),
            "target": selected_final[:, 1].astype(np.int64),
            "length": lengths,
            "weight": weights,
        })
        return df

    # ------------------------------------------------------------------
    # Connectivity helper
    # ------------------------------------------------------------------

    def _bridge_components(
        self,
        xy: np.ndarray,
        parent: np.ndarray,
        n: int,
        tree: "cKDTree",
    ) -> list[tuple[int, int]]:
        """Add minimum bridge edges so the graph has a single connected component.

        Uses a Prim-like approach: repeatedly find the shortest valid cross-
        component edge and add it, merging the two components.  Repeats until
        one component remains.

        Bridge edges are allowed to exceed L_MAX (they are a last resort to
        guarantee connectivity).  Core-crossing edges are avoided if possible;
        if no non-crossing bridge exists, a crossing edge is used rather than
        leaving the graph disconnected.

        Parameters
        ----------
        xy     : (N, 2) array of node positions.
        parent : union-find parent array (modified in place as components merge).
        n      : number of nodes.
        tree   : cKDTree of all node positions.

        Returns
        -------
        List of (source, target) pairs for the bridge edges.
        """
        def _find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        bridge_list: list[tuple[int, int]] = []
        # How many nearest neighbours to examine per node per iteration.
        k_search = min(50, n)

        while True:
            # Gather one representative per component.
            rep_of: dict[int, int] = {}
            for i in range(n):
                r = _find(i)
                if r not in rep_of:
                    rep_of[r] = i

            if len(rep_of) <= 1:
                break   # single component – done

            reps = list(rep_of.values())
            rep_xy = xy[reps]

            # Query k nearest neighbours for every representative.
            dists, idxs = tree.query(rep_xy, k=k_search)

            best_dist  = float("inf")
            best_pair: tuple[int, int] | None = None
            best_crosses = True   # prefer non-crossing edges

            for ri, (row_d, row_idx) in enumerate(zip(dists, idxs)):
                src      = reps[ri]
                src_root = _find(src)
                for d, tgt in zip(row_d, row_idx):
                    tgt = int(tgt)
                    if _find(tgt) == src_root:
                        continue   # same component – skip
                    crosses = bool(
                        segments_cross_circle(
                            xy[src : src + 1], xy[tgt : tgt + 1],
                            self.cfg.r_core,
                        )[0]
                    )
                    # Prefer non-crossing; among equal crossing status prefer shorter.
                    if (not crosses and (best_crosses or d < best_dist)) or \
                       (crosses and best_crosses and d < best_dist):
                        best_dist    = d
                        best_pair    = (src, tgt)
                        best_crosses = crosses
                    if best_pair and not best_crosses:
                        break   # already have a valid non-crossing edge

            if best_pair is None:
                print("  WARNING: could not find any bridge edge; "
                      "some components remain disconnected.")
                break

            bridge_list.append(best_pair)
            ps = _find(best_pair[0])
            pt = _find(best_pair[1])
            parent[ps] = pt

        return bridge_list

    # ------------------------------------------------------------------
    # Stage C: worldbuilding attribute assignment
    # ------------------------------------------------------------------

    def _generate_uids(self) -> np.ndarray:
        """Generate unique 7-digit system registry IDs (as zero-padded strings).

        Samples integers in [1_000_000, 9_999_999] with collision rejection.
        For n ≤ 20_000 out of 9_000_000 possible values, the birthday-paradox
        collision rate is < 0.02 %, so a single batch almost always suffices.
        """
        n = self.cfg.n_nodes
        seen: set[int] = set()
        result: list[str] = []
        while len(result) < n:
            batch = self._rng.integers(1_000_000, 10_000_000, size=n + 50)
            for v in batch:
                vi = int(v)
                if vi not in seen:
                    seen.add(vi)
                    result.append(f"{vi:07d}")
                    if len(result) >= n:
                        break
        return np.array(result[:n])

    def _assign_population(
        self,
        nodes: pd.DataFrame,
        edges: pd.DataFrame,
        degree: np.ndarray,
    ) -> np.ndarray:
        """Assign population values 0-100 to each node.

        Algorithm
        ---------
        1. Isolated nodes (degree = 0) receive pop = 0.
        2. Connected nodes start with independent samples from Normal(50, 18).
        3. A radial bonus shifts inner nodes toward higher values:
              radial_bonus = (1 − r/r_disk) × 40 × pop_core_dispersal
        4. Spatial clustering: blends each node's score with a weighted average
           of its direct neighbours' scores.  Strength = min(0.6, 0.3 × pop_dispersal).
        5. Scores are clamped to [1, 100] and rounded to integers.
        """
        cfg = self.cfg
        n = len(nodes)
        r = nodes["r"].values
        connected = degree > 0

        # Step 2: base bell-curve scores
        raw = self._rng.normal(50.0, 18.0, n)

        # Step 3: radial bonus
        r_norm = np.clip(r / max(cfg.r_disk, 1e-9), 0.0, 1.0)
        radial_bonus = (1.0 - r_norm) * 40.0 * cfg.pop_core_dispersal
        raw = raw + radial_bonus

        # Step 4: spatial clustering via sparse adjacency matrix-vector product
        if cfg.pop_dispersal > 0 and len(edges) > 0:
            src = edges["source"].values.astype(np.int64)
            tgt = edges["target"].values.astype(np.int64)
            # Build symmetric sparse adjacency (undirected)
            data = np.ones(2 * len(src), dtype=np.float64)
            row  = np.concatenate([src, tgt])
            col  = np.concatenate([tgt, src])
            A = csr_matrix((data, (row, col)), shape=(n, n))
            deg_vec = np.array(A.sum(axis=1)).flatten()
            neighbor_sum = A.dot(raw)
            with np.errstate(divide="ignore", invalid="ignore"):
                neighbor_avg = np.where(deg_vec > 0, neighbor_sum / deg_vec, raw)
            w = min(0.6, 0.3 * cfg.pop_dispersal)
            raw = raw * (1.0 - w) + neighbor_avg * w

        # Step 5: finalise
        pop = np.round(np.clip(raw, 1.0, 100.0)).astype(np.int64)
        pop[~connected] = 0   # isolated nodes are uninhabited

        return pop

    def _resolve_admin_counts(self, n_connected: int) -> Tuple[int, int, int, int, int]:
        """Return (n5, n4, n3, n2, n1) admin node counts for each level.

        Exact counts from GalaxyConfig take precedence over ratio-based
        calculation.  The total is capped at n_connected.
        """
        cfg = self.cfg

        # If all explicit counts are non-negative, use them directly
        explicit = [cfg.admin_count_1, cfg.admin_count_2, cfg.admin_count_3,
                    cfg.admin_count_4, cfg.admin_count_5]
        if all(c >= 0 for c in explicit):
            n1, n2, n3, n4, n5 = explicit
        else:
            # Ratio-based calculation targeting ~5% of connected nodes.
            # Overrides are applied in cascade order: n1 override happens before
            # n2 is computed from n1, so n2=ratio uses the effective n1 value.
            target_total = max(1, int(n_connected * 0.05))
            r21, r32, r43, r54 = (cfg.admin_ratio_21, cfg.admin_ratio_32,
                                   cfg.admin_ratio_43, cfg.admin_ratio_54)
            multiplier = 1.0 + r21 + r21*r32 + r21*r32*r43 + r21*r32*r43*r54
            n1 = max(1, int(round(target_total / multiplier)))
            if cfg.admin_count_1 >= 0:
                n1 = cfg.admin_count_1
            n2 = max(1, int(round(n1 * r21)))
            if cfg.admin_count_2 >= 0:
                n2 = cfg.admin_count_2
            n3 = max(1, int(round(n2 * r32)))
            if cfg.admin_count_3 >= 0:
                n3 = cfg.admin_count_3
            n4 = max(1, int(round(n3 * r43)))
            if cfg.admin_count_4 >= 0:
                n4 = cfg.admin_count_4
            n5 = max(1, int(round(n4 * r54)))
            if cfg.admin_count_5 >= 0:
                n5 = cfg.admin_count_5

        # Cap total at available connected nodes.
        # Scale all levels proportionally, but always keep at least one
        # level-1 node so the hierarchy system has a root to propagate from.
        total = n1 + n2 + n3 + n4 + n5
        if total > n_connected:
            scale = n_connected / max(total, 1)
            n5 = int(n5 * scale)
            n4 = int(n4 * scale)
            n3 = int(n3 * scale)
            n2 = int(n2 * scale)
            n1 = max(1, int(n1 * scale))

        return n5, n4, n3, n2, n1

    # ------------------------------------------------------------------
    # Admin placement helpers
    # ------------------------------------------------------------------

    def _build_adj_list(self, n: int, edges: pd.DataFrame) -> list[list[int]]:
        """Build an undirected adjacency list from *edges*."""
        adj: list[list[int]] = [[] for _ in range(n)]
        if len(edges) > 0:
            src = edges["source"].values.astype(int)
            tgt = edges["target"].values.astype(int)
            for s, t in zip(src, tgt):
                adj[s].append(t)
                adj[t].append(s)
        return adj

    def _estimate_graph_spread(
        self,
        adj: list[list[int]],
        connected_idx: np.ndarray,
    ) -> Tuple[float, float]:
        """Estimate graph diameter and average degree via BFS sampling.

        Returns
        -------
        diameter_estimate : approximate longest BFS distance found
        avg_degree        : mean degree of connected nodes
        """
        n_conn = len(connected_idx)
        if n_conn == 0:
            return 1.0, 1.0

        avg_deg = sum(len(adj[i]) for i in connected_idx) / n_conn

        # BFS from up to 7 random connected nodes; take the maximum eccentricity.
        sample_size = min(7, n_conn)
        samples = self._rng.choice(connected_idx, size=sample_size, replace=False)

        max_dist = 1
        for start in samples:
            visited: dict[int, int] = {int(start): 0}
            q: deque[int] = deque([int(start)])
            while q:
                node = q.popleft()
                d    = visited[node]
                for nb in adj[node]:
                    if nb not in visited:
                        visited[nb] = d + 1
                        q.append(nb)
            if visited:
                max_dist = max(max_dist, max(visited.values()))

        return float(max_dist), float(avg_deg)

    def _bfs_neighborhood(
        self, adj: list[list[int]], start: int, max_hops: int
    ) -> dict[int, int]:
        """BFS from *start* up to *max_hops* hops. Returns {node: distance}."""
        visited: dict[int, int] = {start: 0}
        q: deque[int] = deque([start])
        while q:
            node = q.popleft()
            d    = visited[node]
            if d >= max_hops:
                continue
            for nb in adj[node]:
                if nb not in visited:
                    visited[nb] = d + 1
                    q.append(nb)
        return visited

    def _greedy_spaced_placement(
        self,
        adj: list[list[int]],
        scores: np.ndarray,
        count: int,
        sep_radius: int,
        excluded: set[int],
    ) -> list[int]:
        """Greedy placement of *count* centres with ≥ sep_radius hop separation.

        Algorithm
        ---------
        Candidates (connected, not yet assigned to any level) are visited in
        descending score order.  A candidate is placed only if it is at least
        *sep_radius* hops from every already-placed centre of this level.
        After placement, all nodes within (sep_radius − 1) hops are added to
        the *too_close* exclusion set so subsequent candidates skip them.

        The BFS for each placed centre only explores the local neighbourhood
        up to sep_radius hops, so the per-level cost is O(count × D^sep_radius)
        where D is the average degree — very fast even for N = 20 000.

        Falls back to progressively relaxed sep_radius if the target count
        cannot be reached with the original constraint.
        """
        if count <= 0:
            return []

        placed: list[int] = []

        for current_sep in range(sep_radius, 0, -1):
            placed = []
            too_close: set[int] = set(excluded)

            # Sort by score desc once; iterate linearly.
            candidates = sorted(
                [i for i in range(len(scores))
                 if scores[i] > 0 and i not in excluded],
                key=lambda i: -scores[i],
            )

            for c in candidates:
                if len(placed) >= count:
                    break
                if c in too_close:
                    continue

                placed.append(c)

                # Exclude all nodes within (current_sep - 1) hops.
                if current_sep > 1:
                    nearby = self._bfs_neighborhood(adj, c, current_sep - 1)
                    too_close.update(nearby.keys())
                else:
                    too_close.add(c)

            if len(placed) >= count:
                break

        return placed

    def _assign_admin_levels(
        self,
        nodes: pd.DataFrame,
        edges: pd.DataFrame,
        pop: np.ndarray,
        degree: np.ndarray,
    ) -> np.ndarray:
        """Assign administrative hierarchy levels 0–5 with spatial separation.

        Algorithm
        ---------
        1.  Build a graph adjacency list and estimate the network's diameter
            and average degree from BFS samples.
        2.  For each level L, auto-compute a *coverage radius* — the BFS
            radius a single centre of that level should ideally service — as::

                r_cover = log(n_connected / n_L) / log(avg_degree) × coverage_scale

            This captures the intuition that if you have n_L centres spread
            across n_connected nodes, each centre "owns" n_connected / n_L
            nodes, and the BFS radius of a ball that size is log(n/n_L)/log(D).

        3.  The *separation radius* (minimum hops between two same-level
            centres) is sep_radius = max(1, round(r_cover × sep_scale)).

        4.  Greedy placement with the sep_radius constraint, processing
            from the most-important level (1) to the least (5).  Each node
            is assigned at most one level; higher levels are excluded from
            subsequent levels' candidate pools.

        Fallback: if the required count cannot be placed under the initial
        sep_radius, the algorithm halves the radius and retries, so it
        always produces *some* result.
        """
        cfg = self.cfg
        n   = len(nodes)
        admin_lvl = np.zeros(n, dtype=np.int64)

        connected_mask = degree > 0
        connected_idx  = np.where(connected_mask)[0]
        n_connected    = len(connected_idx)
        if n_connected == 0:
            return admin_lvl

        # Build adjacency structure and graph spread estimates.
        adj = self._build_adj_list(n, edges)
        diameter, avg_deg = self._estimate_graph_spread(adj, connected_idx)
        eff_deg = max(avg_deg, 1.5)   # guard against degree-0/1 graphs

        # Scores: pop (60 %) + normalised degree (40 %).
        max_deg = max(int(degree.max()), 1)
        scores  = np.zeros(n, dtype=np.float64)
        scores[connected_mask] = (
            pop[connected_mask] * 0.6
            + (degree[connected_mask] / max_deg) * 100.0 * 0.4
        )

        n5, n4, n3, n2, n1 = self._resolve_admin_counts(n_connected)

        def _sep_for(target_count: int) -> int:
            if target_count <= 0:
                return 1
            nodes_per = max(n_connected / target_count, 2.0)
            r_cover   = math.log(nodes_per) / math.log(eff_deg)
            r_scaled  = max(1.0, r_cover * cfg.admin_coverage_scale)
            sep       = max(1, int(round(r_scaled * cfg.admin_sep_scale)))
            return sep

        sep1 = _sep_for(n1)
        sep2 = _sep_for(n2)
        sep3 = _sep_for(n3)
        sep4 = _sep_for(n4)
        sep5 = _sep_for(n5)

        print(f"  Admin sep radii (hops):  "
              f"lvl1={sep1}  lvl2={sep2}  lvl3={sep3}  "
              f"lvl4={sep4}  lvl5={sep5}  "
              f"(graph diam≈{int(diameter)}, avg_deg≈{avg_deg:.1f})")

        # Greedy placement: highest priority first; each assigned node is
        # excluded from consideration for lower-priority levels.
        assigned: set[int] = set()
        for lvl, count, sep in [
            (1, n1, sep1), (2, n2, sep2), (3, n3, sep3),
            (4, n4, sep4), (5, n5, sep5),
        ]:
            if count <= 0:
                continue
            placed = self._greedy_spaced_placement(adj, scores, count, sep, assigned)
            for node in placed:
                admin_lvl[node] = lvl
                assigned.add(node)

        return admin_lvl

    def _estimate_betweenness(
        self,
        adj: list[list[int]],
        n: int,
        connected_idx: np.ndarray,
    ) -> np.ndarray:
        """Approximate betweenness centrality via sampled Brandes BFS.

        Runs the full Brandes single-source algorithm on a random sample of
        ``k = min(n_conn, max(30, 2 × √n_conn))`` source nodes.  The results
        are scaled by ``n_conn / k`` so the magnitude approximates true
        betweenness up to that sampling factor.

        Returns
        -------
        btw : np.ndarray, shape (n,)
            Approximate betweenness scores.  Isolated nodes stay at 0.
        """
        n_conn  = len(connected_idx)
        k       = min(n_conn, max(30, int(n_conn ** 0.5) * 2))
        sources = self._rng.choice(connected_idx, size=k, replace=False)

        btw = np.zeros(n, dtype=np.float64)

        for s_raw in sources:
            s = int(s_raw)

            sigma: list[int]   = [0]  * n   # shortest-path counts
            dist:  list[int]   = [-1] * n   # hop distances (-1 = unvisited)
            delta: list[float] = [0.0] * n  # pair-dependency accumulator
            pred:  list[list[int]] = [[] for _ in range(n)]

            sigma[s] = 1
            dist[s]  = 0
            order: list[int] = []
            q = deque([s])

            # Forward BFS — build shortest-path DAG
            while q:
                u = q.popleft()
                order.append(u)
                for v in adj[u]:
                    if dist[v] < 0:          # first visit
                        dist[v] = dist[u] + 1
                        q.append(v)
                    if dist[v] == dist[u] + 1:   # on a shortest path
                        sigma[v] += sigma[u]
                        pred[v].append(u)

            # Backward accumulation (Brandes)
            for v in reversed(order):
                for u in pred[v]:
                    if sigma[v] > 0:
                        delta[u] += (sigma[u] / sigma[v]) * (1.0 + delta[v])
                if v != s:
                    btw[v] += delta[v]

        # Scale to compensate for sampling fraction
        if k < n_conn:
            btw *= (n_conn / k)

        return btw

    def _assign_chokepoints(
        self,
        nodes: pd.DataFrame,
        edges: pd.DataFrame,
        degree: np.ndarray,
    ) -> np.ndarray:
        """Designate chokepoint nodes based on betweenness centrality.

        A chokepoint is a node that lies on a disproportionate number of
        shortest paths through the network — it provides *exclusive access*
        to large sectors of the galaxy.  Chokepoints are spread out so that
        no two are within ``choke_sep`` hops of each other.

        Parameters
        ----------
        choke_count == 0  → all False (disabled, default)
        choke_count == -1 → random-chance mode: roll ``choke_chance``; if it
                            succeeds pick a random count ∈ [1, n_conn//300+1]
        choke_count > 0   → exactly that many chokepoints

        Returns
        -------
        is_choke : np.ndarray[bool], shape (n,)
        """
        cfg      = self.cfg
        n        = len(nodes)
        is_choke = np.zeros(n, dtype=bool)

        if cfg.choke_count == 0:
            return is_choke

        connected_mask = degree > 0
        connected_idx  = np.where(connected_mask)[0]
        n_connected    = len(connected_idx)
        if n_connected < 3:
            return is_choke

        # Determine target count ────────────────────────────────────────────
        if cfg.choke_count == -1:
            if self._rng.random() >= cfg.choke_chance:
                print("  Chokepoints: random roll failed — none placed.")
                return is_choke
            max_count = max(1, n_connected // 300)
            count = int(self._rng.integers(1, max_count + 2))
        else:
            count = int(cfg.choke_count)

        count = max(1, min(count, n_connected))

        # Build graph structures ────────────────────────────────────────────
        adj = self._build_adj_list(n, edges)

        if cfg.choke_sep == -1:
            diameter, _ = self._estimate_graph_spread(adj, connected_idx)
            sep_radius  = max(3, int(diameter // 5))
        else:
            sep_radius = max(1, int(cfg.choke_sep))

        # Score by approximate betweenness centrality ───────────────────────
        scores = self._estimate_betweenness(adj, n, connected_idx)

        # Greedy spaced placement ───────────────────────────────────────────
        placed = self._greedy_spaced_placement(
            adj, scores, count, sep_radius, excluded=set()
        )

        is_choke[placed] = True
        print(f"  {len(placed)} chokepoint(s) designated "
              f"(target {count}, sep ≥ {sep_radius} hops).")
        return is_choke

    def assign_attributes(
        self,
        nodes: pd.DataFrame,
        edges: pd.DataFrame,
    ) -> pd.DataFrame:
        """Compute and attach all worldbuilding attributes to *nodes*.

        This method can be called independently of ``run()`` to re-assign
        attributes (e.g. when the user tweaks population or admin parameters
        in the GUI without regenerating the spatial layout).

        Parameters
        ----------
        nodes : DataFrame  – spatial columns already present
        edges : DataFrame  – FTL lane connections

        Returns
        -------
        nodes : DataFrame with uid, name, pop, admin_lvl, admin_dist,
                hierarchy, is_choke columns added (or replaced if already
                present).
        """
        n = len(nodes)
        degree = compute_degree(n, edges)

        print("Stage C: assigning worldbuilding attributes …")

        # UIDs – only generate if not already present
        if "uid" not in nodes.columns:
            nodes = nodes.copy()
            nodes["uid"] = self._generate_uids()
        if "name" not in nodes.columns:
            nodes = nodes.copy() if "uid" in nodes.columns else nodes
            nodes["name"] = ""

        pop       = self._assign_population(nodes, edges, degree)
        admin_lvl = self._assign_admin_levels(nodes, edges, pop, degree)

        nodes = nodes.copy()
        nodes["pop"]       = pop
        nodes["admin_lvl"] = admin_lvl

        # admin_dist requires admin_lvl column
        nodes["admin_dist"] = compute_admin_dist(nodes, edges)

        # hierarchy requires admin_lvl, pop, and edges
        nodes["hierarchy"] = compute_hierarchy(nodes, edges)

        # is_choke — betweenness-based exclusive-access bottlenecks
        nodes["is_choke"] = self._assign_chokepoints(nodes, edges, degree)

        n_hier = int((nodes["hierarchy"] >= 0).sum())
        n_choke = int(nodes["is_choke"].sum())
        print(f"  {int((pop > 0).sum()):,} inhabited systems, "
              f"{int((admin_lvl > 0).sum()):,} admin centres, "
              f"{n_hier:,} nodes placed in hierarchies, "
              f"{n_choke} chokepoint(s).")
        return nodes

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
        attr_cols = ["x", "y", "r", "theta", "arm_dist",
                     "uid", "name", "pop", "admin_lvl", "admin_dist", "hierarchy"]
        for row in nodes.itertuples(index=False):
            attrs = {col: getattr(row, col)
                     for col in attr_cols if hasattr(row, col)}
            # Ensure numeric types for networkx
            for k in ["x", "y", "r", "theta", "arm_dist"]:
                if k in attrs:
                    attrs[k] = float(attrs[k])
            for k in ["pop", "admin_lvl", "admin_dist", "hierarchy"]:
                if k in attrs:
                    attrs[k] = int(attrs[k])
            G.add_node(int(row.id), **attrs)

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
        C – assign worldbuilding attributes (uid, name, pop, admin_lvl,
            admin_dist).

        Returns
        -------
        nodes_df : DataFrame  (id, x, y, r, theta, arm_dist,
                               uid, name, pop, admin_lvl, admin_dist, hierarchy)
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

        # ── Stage C ──────────────────────────────────────────────────
        print()
        t0 = time.perf_counter()
        nodes = self.assign_attributes(nodes, edges)
        print(f"  Attributes assigned in {time.perf_counter() - t0:.2f}s")

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
