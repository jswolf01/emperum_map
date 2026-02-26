# Emperum Galaxy Generator — GUI Guide

Procedural generator for a face-on galactic disk with spiral arms, a forbidden
core, L-way hyperspace edges, and a full worldbuilding attribute stack — for
the **Emperum** sci-fi setting.

Outputs `nodes.csv`, `edges.csv`, and (optionally) `graph.gexf`.

---

## Project layout

```
emperum_map/
├── galaxy_gui.py     Desktop GUI (Tkinter + embedded Matplotlib)
├── galaxygen.py      Core generator: GalaxyConfig + GalaxyGenerator
├── run_generate.py   CLI entry point (argparse wrapper)
├── plot_debug.py     Standalone Matplotlib plot / export tool
├── readme.md         This file — GUI-focused guide
├── README.md         Full technical reference
└── output/           Created at runtime
    ├── nodes.csv
    ├── edges.csv
    ├── params.json
    └── graph.gexf    (if networkx is installed)
```

---

## Installation

Python 3.9+ required.

```bash
pip install numpy pandas scipy matplotlib networkx
```

`networkx` is only needed for GEXF export; all other features work without it.
For the GUI, `tkinter` is also required (bundled with most Python installers;
on Ubuntu/Debian: `sudo apt-get install python3-tk`).

---

## Launching the GUI

```bash
python galaxy_gui.py
```

A three-panel window opens:

| Panel | Contents |
|-------|----------|
| **Left** | Generation parameters — scrollable, collapsible sections |
| **Centre** | Embedded Matplotlib preview with zoom / pan |
| **Right** | Worldbuilding attributes, Node Inspector, Search |

The action bar runs along the bottom: **Generate**, **Preview**,
**Export PNG…**, **Export SVG…**, and a status line.

---

## Typical workflow

1. **Tune parameters** in the left panel (collapse sections you are not using).
2. Click **Generate** — builds the spatial graph and writes all output files.
3. Click **Preview** — renders the galaxy in the centre panel.
4. Use the **Right panel** to adjust population / admin settings, then click
   **Apply Attributes** to re-run worldbuilding assignment without rebuilding
   the spatial layout.
5. Switch `Color by` to `hierarchy` (or any other mode) and click **Preview**
   again — no regeneration needed, just a visual redraw.
6. Click a node to inspect it in the **Node Inspector**.
7. **Export PNG…** or **Export SVG…** to save the current view.

---

## Left panel — generation parameters

### Spatial

| Control | Default | Description |
|---------|---------|-------------|
| Star systems (n_nodes) | 5 000 | Total node count |
| Disk radius (r_disk) | 100.0 | Outer boundary |
| Core radius (r_core) | 15.0 | Forbidden core; no nodes or edge crossings inside |

### Edges (Hyperspace Lanes)

| Control | Default | Description |
|---------|---------|-------------|
| Max lane length (l_max) | 9.0 | Hard length cap per edge |
| Avg degree (target_degree) | 4.0 | Edge budget ≈ N × D / 2 |
| Connection chance | 1.0 | Fraction of nodes eligible for any edges (0 = all isolated, 1 = all eligible) |

### Spiral Arms

| Control | Default | Description |
|---------|---------|-------------|
| Number of arms | 4 | Logarithmic spiral count |
| Arm tightness (arm_b) | 0.35 | Smaller = tighter coil, larger = more open |
| Arm width σ | 5.5 | Gaussian half-width around each centreline |
| Inter-arm density | 0.15 | Baseline acceptance between arms (0 = empty gaps, 1 = uniform) |

### Radial Density

| Control | Default | Description |
|---------|---------|-------------|
| Scale length (r_scale) | 38.0 | Exponential falloff; smaller = denser core |

### Sampling & Reproducibility

| Control | Default | Description |
|---------|---------|-------------|
| Boost multiplier | 6.0 | Scales acceptance probability; raise if generation stalls |
| Random seed | 7 | Any integer gives a fully reproducible galaxy |

### Output

| Control | Default | Description |
|---------|---------|-------------|
| Output directory | `output` | Directory for all written files |
| Export GEXF | ✓ | Writes `graph.gexf` for Gephi (requires networkx) |

### Node Appearance

| Control | Default | Description |
|---------|---------|-------------|
| Color by | `arm_dist` | Data driving node color — see table below |
| Node size | 1.5 | Scatter marker size |
| Flat color | `#aaccff` | Used when `Color by = none` |
| Gradient low | `#ffe8c0` | Color at the low end of the gradient |
| Gradient high | `#0d000f` | Color at the high end of the gradient |

**Color by options:**

| Value | Description |
|-------|-------------|
| `arm_dist` | Distance to nearest spiral arm — arm nodes bright, inter-arm dark |
| `r` | Galactic radius — inner bright, outer dark |
| `pop` | Population index — dark (0) through blue/cyan to white (100) |
| `admin_lvl` | Admin level — dark (none) through brown/orange to yellow (5) |
| `admin_dist` | Hop-count to nearest admin node — yellow (0) fading to dark |
| `hierarchy` | **Administrative domain** — each discrete hierarchy gets its own color; all nodes and edges in a domain share that color |
| `none` | Flat single color |

> When `hierarchy` is selected, **edges** are also colored by domain.
> Intra-domain edges are tinted in the domain color; cross-domain edges are
> dimmed to near-invisible.  This is the primary way to visualize the
> administrative hierarchy structure.

### Edge Appearance

| Control | Default | Description |
|---------|---------|-------------|
| Hide edges | ☐ | Skip edge rendering (strongly recommended for N > 1 000) |
| Edge color | `#2244aa` | Flat color (ignored in `hierarchy` color mode) |
| Edge width | 0.4 | Line width in points |
| Edge opacity | 0.35 | 0 = invisible, 1 = fully opaque |

---

## Right panel — worldbuilding attributes

### Apply Attributes

Re-runs the entire worldbuilding stack (population, admin levels, hierarchy)
on the existing spatial layout without rebuilding the graph.  Use this after
tweaking population or admin parameters without wanting to re-generate
positions and edges.

### Population

Population is assigned only to connected nodes (isolated nodes receive `pop=0`).

| Control | Default | Description |
|---------|---------|-------------|
| Core dispersal | 1.0 | 0 = uniform across disk; higher = inner systems richer |
| Clustering | 1.0 | 0 = independent; higher = neighbours pull each other's pop up |

### Admin Levels

Five administrative levels are placed hierarchically.  Level 1 is the highest
(fewest, most powerful); level 5 is the lowest (most numerous).

**Exact counts** (`-1` = use ratios instead):

| Spinbox | Meaning |
|---------|---------|
| Lvl 1 count | Target number of top-tier capitals |
| Lvl 2 count | … |
| … | … |
| Lvl 5 count | Target number of regional logistics centres |

**Ratios** (used when counts are `-1`):

| Control | Default | Description |
|---------|---------|-------------|
| Lvl 2 : Lvl 1 | 2.0 | How many lvl-2 nodes per lvl-1 node |
| Lvl 3 : Lvl 2 | 3.0 | … |
| Lvl 4 : Lvl 3 | 4.0 | … |
| Lvl 5 : Lvl 4 | 5.0 | … |

**Spatial spacing:**

| Control | Default | Description |
|---------|---------|-------------|
| Coverage scale | 1.0 | Scales the auto-computed coverage radius per level |
| Separation scale | 1.5 | Minimum separation = coverage_radius × this value |

### Node Inspector

Click any node on the preview to populate the inspector.

| Field | Editable | Description |
|-------|----------|-------------|
| System ID (uid) | No | Unique 7-digit registry number |
| Admin dist (hops) | No | Hop-count to nearest admin node |
| Hierarchy index | No | Integer ID of the admin domain this node belongs to |
| Name | **Yes** | Free-text system name |
| Population (0–100) | **Yes** | Population index |
| Admin level (0–5) | **Yes** | Administrative level (0 = none) |

Buttons:

- **Save Node Edits** — writes the editable fields back to `nodes.csv`.
- **Recalc Admin Dist** — recomputes hop-distances from current `admin_lvl`
  values (useful after manual admin edits).
- **Deselect** — clears the selection.

> After manually changing `admin_lvl` values, use **Apply Attributes** (top of
> right panel) to recompute hierarchy assignments across all nodes.

### Search

Search any column in the node table using a chosen attribute, operator, and
value.  Matching nodes are highlighted as green rings on the preview.

| Attribute options | `name`, `uid`, `pop`, `admin_lvl`, `admin_dist`, `hierarchy`, `id` |
|---|---|
| Operator options | `contains`, `=`, `>`, `<`, `>=`, `<=` |

The first matched node is auto-selected in the inspector.

---

## Preview canvas controls

| Gesture | Action |
|---------|--------|
| Scroll wheel | Zoom in / out centred on cursor |
| Right-click drag | Pan |
| Middle-click drag | Pan |
| Left-click | Select nearest node (within ~1.5 % of current view width) |

The Matplotlib navigation toolbar (just below the canvas) provides additional
zoom-to-rect, pan, and home/back/forward controls.

---

## Administrative hierarchy — how it works

After admin levels are assigned, `compute_hierarchy()` traces a chain from
each populated node up to its nearest admin-level-1 root:

```
non-admin node → nearest lvl-5 node → nearest lvl-4 → ... → lvl-1 root
```

"Nearest" is measured in **graph hops**, not Euclidean distance.  If an
intermediate level is absent (e.g. admin_lvl 4 has count 0), it is skipped.

Each admin-level-1 node is the root of exactly one discrete hierarchy.  All
nodes whose chain terminates at the same lvl-1 root share the same integer
`hierarchy` value (0-based, assigned in ascending node-ID order of lvl-1
roots).  Unpopulated nodes and nodes disconnected from every admin-level-1
root receive `hierarchy = -1`.

To visualize: select `hierarchy` from the **Color by** dropdown and click
**Preview**.

---

## Output files

### `nodes.csv`

| Column | Type | Description |
|--------|------|-------------|
| `id` | int | Unique node index (0-based) |
| `x`, `y` | float | Cartesian position |
| `r` | float | Distance from galactic centre |
| `theta` | float | Azimuth in radians (0–2π) |
| `arm_dist` | float | Approximate distance to nearest spiral arm centreline |
| `uid` | str | 7-digit zero-padded system registry number |
| `name` | str | Free-text name (blank by default, editable in inspector) |
| `pop` | int | Population index 0–100 (0 = uninhabited / isolated) |
| `admin_lvl` | int | Administrative level 0–5 (0 = none) |
| `admin_dist` | int | Hop-count to nearest admin node (−1 if unreachable) |
| `hierarchy` | int | Admin-domain index (−1 if unassigned) |

### `edges.csv`

| Column | Type | Description |
|--------|------|-------------|
| `source` | int | Node id of one endpoint |
| `target` | int | Node id of other endpoint |
| `length` | float | Euclidean distance between endpoints |
| `weight` | float | 1 / length |

### `params.json`

All generation parameters saved as JSON.  `plot_debug.py` reads this
automatically so overlay graphics always match the generation run.

### `graph.gexf` *(optional)*

Full graph for Gephi with all node attributes (including `hierarchy`) and edge
attributes.  Written when the GEXF checkbox is enabled and `networkx` is
installed.

---

## Command-line alternative

If you prefer a headless workflow:

```bash
# Generate
python run_generate.py --n_nodes 5000 --seed 7

# Visualize (all color modes available)
python plot_debug.py --color_by hierarchy --no_edges
python plot_debug.py --color_by admin_lvl
python plot_debug.py --color_by pop --save galaxy_pop.png
```

See `README.md` for the full CLI reference.
