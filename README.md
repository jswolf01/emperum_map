# Emperum Galaxy Generator

Procedural generator for a face-on galactic disk with spiral arms, a forbidden
core, and L-way hyperspace edges — for the **Emperum** sci-fi worldbuilding
project.

Outputs two CSVs (importable directly into Gephi / QGIS) and an optional GEXF
file.  Everything is seeded and fully reproducible.

---

## Project layout

```
emperum_map/
├── galaxygen.py      Core generator: GalaxyConfig dataclass + GalaxyGenerator
├── run_generate.py   CLI entrypoint (argparse wrapper)
├── plot_debug.py     Matplotlib sanity-check plot
├── README.md         This file
└── output/           Created at runtime
    ├── nodes.csv
    ├── edges.csv
    └── graph.gexf    (if networkx is installed)
```

---

## Installation

Python 3.9+ is required.  Install dependencies with pip:

```bash
pip install numpy pandas scipy matplotlib networkx
```

`networkx` is only needed for GEXF export; all other features work without it.

---

## Quick start

```bash
# Generate with default parameters (5 000 nodes, seed 7)
python run_generate.py

# Then view the debug plot (skip edges for speed at N=5000)
python plot_debug.py --no_edges

# Save the plot to a PNG
python plot_debug.py --no_edges --save galaxy.png
```

---

## CLI reference — `run_generate.py`

| Flag | Default | Description |
|------|---------|-------------|
| `--n_nodes N` | 5000 | Number of star-system nodes |
| `--r_disk R` | 100.0 | Outer disk boundary radius |
| `--r_core R` | 15.0 | Forbidden core radius |
| `--l_max L` | 9.0 | Maximum edge length |
| `--target_degree D` | 4.0 | Target average node degree (approximate) |
| `--n_arms N` | 4 | Number of logarithmic spiral arms |
| `--arm_b B` | 0.35 | Spiral tightness (larger = more open) |
| `--arm_sigma S` | 5.5 | Arm Gaussian half-width |
| `--arm_base P` | 0.15 | Baseline acceptance away from arms (0–1) |
| `--r_scale R` | 38.0 | Radial density scale length |
| `--boost B` | 6.0 | Acceptance-probability multiplier |
| `--seed S` | 7 | Random seed |
| `--out_dir DIR` | output | Output directory |
| `--no_gexf` | (off) | Skip GEXF export |

Example with the exact canonical parameters:

```bash
python run_generate.py \
    --n_nodes 5000 \
    --r_disk 100 \
    --r_core 15 \
    --l_max 9 \
    --target_degree 4 \
    --n_arms 4 \
    --arm_sigma 5.5 \
    --arm_b 0.35 \
    --r_scale 38 \
    --seed 7 \
    --out_dir output
```

---

## CLI reference — `plot_debug.py`

| Flag | Default | Description |
|------|---------|-------------|
| `--out_dir DIR` | output | Directory with nodes.csv / edges.csv |
| `--no_edges` | (off) | Skip drawing edges (recommended for N > 1000) |
| `--color_by` | arm_dist | Node colour: `arm_dist`, `r`, or `none` |
| `--node_size S` | 1.5 | Scatter marker size |
| `--edge_alpha A` | 0.35 | Edge line transparency |
| `--save FILE` | (show) | Save to file instead of interactive window |
| `--r_disk`, `--r_core`, `--n_arms`, `--arm_b`, `--r_arm_start` | (matching defaults) | Must match generation run |

---

## Output files

### `output/nodes.csv`

| Column | Type | Description |
|--------|------|-------------|
| `id` | int | Unique node identifier |
| `x`, `y` | float | Cartesian position |
| `r` | float | Distance from galactic centre |
| `theta` | float | Azimuth (radians, 0–2π) |
| `arm_dist` | float | Approximate distance to nearest spiral arm centreline |

### `output/edges.csv`

| Column | Type | Description |
|--------|------|-------------|
| `source` | int | Node id of one endpoint |
| `target` | int | Node id of other endpoint |
| `length` | float | Euclidean distance between endpoints |
| `weight` | float | 1 / length (useful as edge weight in graph analysis) |

### `output/graph.gexf`

NetworkX-generated GEXF 1.2 file containing all nodes with their spatial
attributes and all edges with `length` and `weight`.  Ready to open directly
in Gephi.

---

## Gephi import steps

### Option A — GEXF (easiest)

1. Open Gephi.
2. **File → Open…** → select `output/graph.gexf`.
3. Accept the dialog (Undirected graph, merge strategy doesn't matter here).
4. Switch to the **Overview** tab.
5. The nodes already carry `x` / `y` coordinates from generation.  In the
   **Layout** panel on the left, choose **"No layout (use node positions)"** or
   simply leave the layout as-is.
6. To apply ForceAtlas2 *on top of* the spatial layout (for visual clarity),
   select **ForceAtlas2**, enable **"Prevent overlap"**, run for ~1000
   iterations, then stop.
7. Run statistics:
   - **Statistics → Network Diameter** (gives betweenness centrality).
   - **Statistics → Modularity** (community detection).
   - **Statistics → Avg. Path Length** (connectivity quality).

### Option B — CSV (manual import)

1. **File → Import Spreadsheet…**
2. Select `output/edges.csv`.
   - Separator: **Comma**
   - Import as: **Edges table**
   - Graph type: **Undirected**
   - Columns: `source` = Source, `target` = Target; import `length` and
     `weight` as edge attributes (Double).
   - Click **Next → Finish**.

3. In the same dialog (or repeat for nodes):
   Select `output/nodes.csv`.
   - Import as: **Nodes table**
   - ID column: `id`
   - Import `x`, `y`, `r`, `theta`, `arm_dist` as node attributes (Double).
   - Click **Next → Finish**.

4. Go to **Data Laboratory** and confirm nodes and edges loaded correctly.

5. Switch to **Overview** tab.

6. Assign node positions from attributes:
   - **Tools → Plugins → (or built-in)**: some Gephi versions let you set
     X/Y directly.  Alternatively, run the **GeoLayout** plugin or use the
     **"Set from attributes"** option if available.
   - Simplest alternative: accept the default layout, run ForceAtlas2 for
     ~500 iterations; the global structure (hubs, corridors) will emerge.

7. Run statistics as in Option A.

### Recommended Gephi colour mappings

| Attribute | Panel | Use for |
|-----------|-------|---------|
| `arm_dist` | Appearance → Nodes → Colour | Shows arm structure |
| `r` | Appearance → Nodes → Colour | Shows radial gradient |
| degree | Appearance → Nodes → Size | Shows hubs |
| Betweenness centrality | Appearance → Nodes → Colour | Reveals chokepoints |
| Modularity class | Appearance → Nodes → Colour | Community structure |

---

## Using as a Python library

```python
from galaxygen import GalaxyConfig, GalaxyGenerator

cfg = GalaxyConfig(
    n_nodes=1000,
    r_disk=100.0,
    r_core=15.0,
    l_max=9.0,
    target_degree=4.0,
    n_arms=4,
    seed=42,
)

gen = GalaxyGenerator(cfg)
nodes_df, edges_df = gen.run()

# nodes_df: pandas DataFrame with id, x, y, r, theta, arm_dist
# edges_df: pandas DataFrame with source, target, length, weight
```

---

## Algorithm overview

### Stage A — Node sampling

Candidate positions are drawn uniformly in the annulus
`[R_CORE, R_DISK]` (area-correct inverse-CDF sampling), then accepted with
probability:

```
P = clip(P_rad × P_arm × BOOST, 0, 1)

P_rad(r) = exp(−r / R_SCALE)
P_arm(d) = ARM_BASE + (1 − ARM_BASE) × exp(−(d / ARM_SIGMA)²)
d        = approximate distance to nearest spiral arm centreline
```

`d` is approximated as `min_k |r − r_arm_k(θ)|`, where `r_arm_k(θ)` is the
logarithmic-spiral radius of arm `k` at azimuth `θ` (multiple windings are
tested).  This is cheaper than exact Euclidean shortest-distance to a curve
and sufficient for density shaping.

### Stage B — Edge generation

1. Build a `cKDTree` from all node positions.
2. Query all pairs within `L_MAX` (O(N log N + E) with the KD-tree).
3. Batch-test every candidate for forbidden-core intersection using a
   vectorised closest-point-on-segment calculation.
4. Shuffle valid candidates; take the first `⌈N × TARGET_DEGREE / 2⌉` edges.

### Forbidden-core intersection test

The segment `P1 → P2` crosses the forbidden core iff the minimum Euclidean
distance from the origin to the segment is less than `R_CORE`:

```
t*     = clamp(−(P1·d) / |d|², 0, 1)   where d = P2 − P1
closest = P1 + t* · d
crosses = |closest|² < R_CORE²
```

This is fully vectorised over all candidate edges at once.

---

## Acceptance tests (printed after generation)

After each run the generator prints:

- Node count == N
- Min r >= R_CORE (tolerance 1e-6)
- Max r <= R_DISK (tolerance 1e-6)
- Max edge length <= L_MAX (tolerance 1e-6)
- Core-crossing count (0 crossings expected)
- Degree distribution: min / median / avg / max
- Number of connected components

---

## Performance notes

| N | Typical wall time (laptop) |
|---|--------------------------|
| 1 000 | < 1 s |
| 5 000 | 2–6 s |
| 10 000 | 5–15 s |
| 50 000 | 30–90 s |

Performance is dominated by `cKDTree.query_pairs` and the batch core-crossing
filter, both of which scale sub-quadratically.

If sampling is slow (acceptance rate too low), increase `--boost`.  If there
are not enough valid edges, increase `--l_max`.

---

## Known limitations and honest caveats

- **Static snapshot**: Logarithmic spiral arms are density waves that rotate
  at a different angular speed than stars.  This generator produces a single
  frozen snapshot, not a dynamically consistent model.

- **Approximate arm distance**: `arm_dist` measures `|r − r_arm(θ)|`, not the
  true shortest Euclidean distance from a point to a spiral curve.  The error
  is small near the arm but can be significant at large arm pitch angles.

- **Edge selection is non-local**: The shuffle-and-select strategy for edges
  does not guarantee uniform degree distribution.  Some nodes may be isolated
  (degree 0) if no valid neighbours exist within `L_MAX`.  Use
  `--target_degree` and `--l_max` together to tune this.

- **Gephi layout distortion**: If you run ForceAtlas2 in Gephi the node
  positions will drift from the generated `x`, `y` values.  Treat the spatial
  layout as an *input reference*, not the final visualisation.

- **2-D projection only**: The galaxy is modelled as a face-on disk.  No
  thickness or inclination is simulated.

- **L-ways are not physical paths**: Connectivity in the Emperum setting is
  defined by barycentric FTL topology, not Euclidean proximity; the edge set
  is a graph-theoretic structure, not a travel-time model.
