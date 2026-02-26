# Emperum Galaxy Generator

Procedural generator for a face-on galactic disk with spiral arms, a forbidden
core, and L-way hyperspace edges — for the **Emperum** sci-fi worldbuilding
project.

Outputs two CSVs (importable directly into Gephi / QGIS) and an optional GEXF
file.  Everything is seeded and fully reproducible.

---

## Plain-English Overview — What Does This Tool Do?

This tool creates a **fictional galaxy map** for sci-fi worldbuilding.  Here is
what it produces, in plain terms:

- **Star systems** — thousands of points scattered across a circular disk.
  More stars are placed near the centre and along spiral arms, like a real
  galaxy.  A hollow "forbidden core" at the centre keeps the very middle empty.

- **Hyperspace lanes (L-ways)** — connections (edges) between nearby star
  systems.  These represent how ships travel between stars.  The tool
  guarantees the whole network is one fully connected map — every system can be
  reached from every other system.

- **Population** — each connected star system is given a population score
  (0–100).  Systems near the galactic centre tend to be more populous.
  Isolated systems with no lanes have population 0.

- **Administrative hierarchy** — some systems are designated as capitals or
  administrative centres at five levels (1 = top-tier capital, 5 = regional
  hub).  These are chosen automatically based on population and connectivity.

- **Admin distance** — for each system, the tool calculates how many hyperspace
  hops it takes to reach the nearest administrative centre.

- **Administrative domains (hierarchies)** — every system is assigned to a
  "domain" belonging to the nearest level-1 capital.  Think of these as
  interstellar nations or sectors.

The interactive GUI lets you tweak all these settings, preview the result,
click individual star systems to inspect and edit their data, search for
systems, and export the map as PNG, SVG, or network files.

---

## Project Layout

```
emperum_map/
├── galaxygen.py      Core generator: GalaxyConfig dataclass + GalaxyGenerator
├── galaxy_gui.py     Interactive GUI (recommended way to use the tool)
├── run_generate.py   CLI entrypoint (command-line use without a GUI)
├── plot_debug.py     Matplotlib visualisation used by the GUI and CLI
├── README.md         This file
└── output/           Created at runtime
    ├── nodes.csv         Star-system data
    ├── edges.csv         Hyperspace lane data
    ├── graph.gexf        Network file for Gephi (if networkx is installed)
    └── params.json       Generation parameters used for this run
```

---

## Installation

Python 3.9+ is required.  Install dependencies:

```bash
pip install numpy pandas scipy matplotlib networkx
```

`networkx` is only needed for GEXF export; all other features work without it.
`tkinter` (for the GUI) is bundled with the standard Python installer.
On Ubuntu/Debian: `sudo apt-get install python3-tk`

---

## Quick Start — GUI (Recommended)

```bash
python galaxy_gui.py
```

1. Click **Generate** to create a galaxy with the current settings.
2. Click **Preview** to render it in the window.
3. Explore the tabs on the left and right to adjust settings.
4. Click any star on the preview to inspect its attributes.
5. Use **Export PNG…** or **Export SVG…** to save the image.

---

## Quick Start — Command Line

```bash
# Generate with default parameters (5 000 nodes, seed 7)
python run_generate.py

# View the result (skip edges for speed at N=5000)
python plot_debug.py --no_edges

# Save the plot to a PNG
python plot_debug.py --no_edges --save galaxy.png
```

---

## GUI — Left Panel

The left panel is split into two tabs:

### Generation Tab

Controls how the galaxy is built.  Changes here require pressing **Generate**
again to take effect.

| Section | What it does |
|---------|-------------|
| **Spatial** | Number of star systems, outer disk size, and forbidden-core radius |
| **Edges (Hyperspace Lanes)** | Maximum lane length, target average connections per system, and connection chance (< 1 creates some isolated systems) |
| **Spiral Arms** | Number of arms, how tightly they wind, their thickness, and how many stars appear between arms |
| **Radial Density** | How sharply star density falls off toward the outer disk |
| **Sampling & Seed** | Boost (acceptance-rate amplifier) and the random seed (same seed = same galaxy every time) |

### Appearance Tab

Controls how the preview looks.  Changes here only require pressing **Preview**
(no regeneration needed).

| Section | What it does |
|---------|-------------|
| **Output** | Where to save CSV/GEXF files and whether to export a GEXF file |
| **Node Appearance** | Color-by mode, node size, flat color (for "none" mode), and gradient low/high colors (for arm_dist and r modes) |
| **Edge Appearance** | Hide edges toggle, edge color, width, and opacity |

**Color-by modes explained:**

| Mode | What you see |
|------|-------------|
| `arm_dist` | Stars near spiral arms get the "low" gradient color; inter-arm stars get the "high" color |
| `r` | Stars near the galactic center get the "low" color; outer stars get the "high" color |
| `pop` | Darker = less populated; brighter = more populated |
| `admin_lvl` | Dark = no admin center; bright yellow = top-tier capital |
| `admin_dist` | Yellow = near an admin center; dark = far away |
| `hierarchy` | Each administrative domain gets a unique color |
| `none` | All stars shown in the flat color you choose |

---

## GUI — Right Panel

The right panel is split into three tabs:

### Attributes Tab

Controls the worldbuilding data assigned to each system.  Press **Apply
Attributes** to re-run these without regenerating the spatial layout.

| Section | What it does |
|---------|-------------|
| **Population** | Core dispersal (how much inner systems are favored) and clustering (how much neighbors influence each other's population) |
| **Admin Levels** | Set exact counts for each admin level (−1 = auto), or use ratios between levels.  Coverage and separation scales control how spread out admin centers are |

### Inspect & Search Tab

**Node Inspector** — click any star on the preview map to fill this panel.
Shows the system's UID, admin distance, and hierarchy index (read-only), plus
editable fields for name, population, and admin level.  Use **Save Node Edits**
to persist changes and **Recalc Admin Dist** to recompute hop distances after
manual edits.

**Search** — find systems by attribute value.  Select an attribute, choose an
operator (`contains`, `=`, `>`, etc.), enter a value, and press **Search**.
Matching systems are highlighted in green on the preview.

### Filter View Tab

The **reduced view** lets you de-clutter the visualization by showing or hiding
star systems based on attribute conditions.

1. Tick **Enable Filter View** to activate it.
2. Set the **Logic** (AND = must match ALL conditions; OR = must match ANY).
3. Set the **Action** (Hide matching OR Show only matching).
4. Use **+ Add Condition** to add one or more filter rules.  Each rule has:
   - **Attribute** — which column to filter on (pop, admin_lvl, r, name, etc.)
   - **Operator** — `=`, `>`, `<`, `>=`, `<=`, `between`, or `contains`
   - **Value** / **To** — the threshold or range (both fields for `between`)
5. Press **Apply Filter** to refresh the preview.

When a node is hidden, all hyperspace lanes connected to it are also hidden.
The map title shows how many systems are currently visible.

Press **Clear All** to remove all conditions and disable the filter.

---

## CLI Reference — `run_generate.py`

| Flag | Default | Description |
|------|---------|-------------|
| `--n_nodes N` | 5000 | Number of star-system nodes |
| `--r_disk R` | 100.0 | Outer disk boundary radius |
| `--r_core R` | 15.0 | Forbidden core radius |
| `--l_max L` | 9.0 | Maximum edge length |
| `--target_degree D` | 4.0 | Target average node degree |
| `--node_connection_chance P` | 1.0 | Probability each node participates in edges |
| `--n_arms N` | 4 | Number of logarithmic spiral arms |
| `--arm_b B` | 0.35 | Spiral tightness (larger = more open) |
| `--arm_sigma S` | 5.5 | Arm Gaussian half-width |
| `--arm_base P` | 0.15 | Baseline density away from arms (0–1) |
| `--r_scale R` | 38.0 | Radial density scale length |
| `--boost B` | 6.0 | Acceptance-probability multiplier |
| `--seed S` | 7 | Random seed |
| `--out_dir DIR` | output | Output directory |
| `--no_gexf` | (off) | Skip GEXF export |

---

## CLI Reference — `plot_debug.py`

| Flag | Default | Description |
|------|---------|-------------|
| `--out_dir DIR` | output | Directory with nodes.csv / edges.csv |
| `--no_edges` | (off) | Skip drawing edges (recommended for N > 1000) |
| `--color_by MODE` | arm_dist | `arm_dist`, `r`, `pop`, `admin_lvl`, `admin_dist`, `hierarchy`, or `none` |
| `--node_size S` | 1.5 | Scatter marker size |
| `--node_color C` | `#aaccff` | Uniform color for `none` mode |
| `--gradient_low_color C` | `#ffe8c0` | Gradient color at low data values |
| `--gradient_high_color C` | `#0d000f` | Gradient color at high data values |
| `--edge_alpha A` | 0.35 | Edge transparency |
| `--edge_color C` | `#2244aa` | Edge color |
| `--edge_width W` | 0.4 | Edge line width |
| `--save FILE` | (show) | Save to PNG/PDF instead of interactive window |
| `--svg [FILE]` | (off) | Save as SVG vector file |
| `--r_disk`, `--r_core`, `--n_arms`, `--arm_b`, `--r_arm_start` | (matching defaults) | Spatial params; auto-loaded from params.json when present |

---

## Output Files

### `output/nodes.csv`

| Column | Type | Description |
|--------|------|-------------|
| `id` | int | Unique node index (0 to N-1) |
| `x`, `y` | float | Cartesian position |
| `r` | float | Distance from galactic centre |
| `theta` | float | Azimuth (radians, 0–2π) |
| `arm_dist` | float | Approximate distance to nearest spiral arm centreline |
| `uid` | str | 7-digit unique system registry number |
| `name` | str | Free-text name (blank by default; editable in the GUI) |
| `pop` | int | Population index 0–100 (0 = uninhabited) |
| `admin_lvl` | int | Administrative level 0–5 (0 = none) |
| `admin_dist` | int | Hop-count to nearest admin centre (−1 if none reachable) |
| `hierarchy` | int | Index of the level-1 domain this system belongs to |

### `output/edges.csv`

| Column | Type | Description |
|--------|------|-------------|
| `source` | int | Node id of one endpoint |
| `target` | int | Node id of other endpoint |
| `length` | float | Euclidean distance between endpoints |
| `weight` | float | 1 / length (useful as edge weight in graph analysis) |

### `output/graph.gexf`

NetworkX-generated GEXF 1.2 file containing all nodes with their spatial and
worldbuilding attributes, and all edges with `length` and `weight`.  Ready to
open directly in Gephi.

---

## Gephi Import Steps

### Option A — GEXF (easiest)

1. Open Gephi.
2. **File → Open…** → select `output/graph.gexf`.
3. Accept the dialog.
4. Switch to the **Overview** tab.
5. In the **Layout** panel, choose **"No layout (use node positions)"** or run
   ForceAtlas2 for ~1000 iterations for visual clarity.
6. Run statistics: **Network Diameter**, **Modularity**, **Avg. Path Length**.

### Option B — CSV (manual import)

1. **File → Import Spreadsheet…** → select `output/edges.csv`.
   - Separator: Comma · Import as: Edges table · Graph type: Undirected
2. Repeat for `output/nodes.csv` (Import as: Nodes table, ID column: `id`).
3. Confirm in **Data Laboratory**, then switch to **Overview**.

---

## Using as a Python Library

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

# nodes_df columns: id, x, y, r, theta, arm_dist, uid, name,
#                   pop, admin_lvl, admin_dist, hierarchy
# edges_df columns: source, target, length, weight
```

---

## Algorithm Overview

### Stage A — Node Sampling

Candidates are drawn uniformly in the annulus `[R_CORE, R_DISK]` and accepted
with probability:

```
P = clip(P_rad × P_arm × BOOST, 0, 1)

P_rad(r) = exp(−r / R_SCALE)
P_arm(d) = ARM_BASE + (1 − ARM_BASE) × exp(−(d / ARM_SIGMA)²)
d        = approximate distance to nearest spiral arm centreline
```

`d` is approximated as `min_k |r − r_arm_k(θ)|` (cheaper than exact Euclidean
distance to the curve and sufficient for density shaping).

### Stage B — Edge Generation

1. Build a `cKDTree` from all node positions.
2. Query all pairs within `L_MAX`.
3. Filter out edges that cross the forbidden core (vectorised).
4. Apply the node-eligibility filter (connection chance).
5. Build a spanning forest with union-find (guarantees one connected component
   among all nodes that have at least one valid pair).
6. **Bridge remaining isolated components** — if multiple disconnected groups
   remain after the spanning forest (e.g., due to a low connection-chance
   "barrier"), the shortest valid cross-component edge is added iteratively
   until the graph is fully connected.  Core-crossing is avoided when possible.
7. Fill remaining edge budget with random extra edges.

### Stage C — Attribute Assignment

1. **UIDs** — unique 7-digit system registry numbers.
2. **Population** — bell curve + radial bonus + neighbour clustering.
3. **Admin levels** — greedy placement at decreasing spatial scales,
   scored by population × connectivity.
4. **Admin distance** — multi-source BFS hop count.
5. **Hierarchy** — each node is assigned to the domain of the nearest
   level-1 capital.

---

## Acceptance Tests (Printed After Generation)

After each run the generator prints:

- Node count == N
- Min r ≥ R_CORE (tolerance 1e-6)
- Max r ≤ R_DISK (tolerance 1e-6)
- Max edge length ≤ L_MAX (bridge edges may exceed this)
- Core-crossing count (0 expected for regular edges)
- Degree distribution: min / median / avg / max
- Number of connected components (1 after bridging)

---

## Performance Notes

| N | Typical wall time (laptop) |
|---|--------------------------|
| 1 000 | < 1 s |
| 5 000 | 2–6 s |
| 10 000 | 5–15 s |
| 50 000 | 30–90 s |

Performance is dominated by `cKDTree.query_pairs` and the batch core-crossing
filter.  If sampling is slow (acceptance rate too low), increase `--boost`.  If
there are not enough valid edges, increase `--l_max`.

---

## Known Limitations and Honest Caveats

- **Static snapshot**: Logarithmic spiral arms are density waves.  This
  generator produces a single frozen snapshot, not a dynamically consistent
  model.

- **Approximate arm distance**: `arm_dist` measures `|r − r_arm(θ)|`, not the
  true Euclidean distance from a point to the spiral curve.  The error is small
  near the arm but can be significant at large pitch angles.

- **Bridge edges may be long**: When the connection-chance filter creates
  barriers, the bridge-building step may add edges that exceed `L_MAX`.  These
  are flagged in the console output.

- **Gradient colors apply only to arm_dist and r modes**: The custom
  gradient-low / gradient-high colors in the GUI control the color scale for
  `arm_dist` and `r` color-by modes.  The `pop`, `admin_lvl`, and `admin_dist`
  modes use built-in semantic colormaps that are fixed.

- **Gephi layout distortion**: Running ForceAtlas2 in Gephi will move nodes
  from their generated `x`, `y` positions.

- **2-D projection only**: The galaxy is modelled as a face-on disk with no
  thickness or inclination.

- **L-ways are not physical paths**: Connectivity is a graph-theoretic
  structure for worldbuilding purposes, not a travel-time model.
