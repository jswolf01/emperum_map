# Emperum Galaxy Generator — GUI Guide

Procedural generator for a face-on galactic disk with spiral arms, a forbidden
core, L-way hyperspace edges, and a full worldbuilding attribute stack — for
the **Emperum** sci-fi setting.

Outputs `nodes.csv`, `edges.csv`, and (optionally) `graph.gexf`.

---

## Plain-language overview

If you are new to the project, start here.  No programming knowledge needed.

**What this tool does:**
It creates a fictional star map.  You tell it how many star systems to place and
how they should be arranged (spiral arms, a dense center, an empty forbidden
zone in the middle), and it generates a galaxy for you — complete with travel
lanes between nearby stars.  You can then assign political attributes to those
stars: population, administrative rank, and which empire or faction controls
each system.

**The three panels:**

| Panel | Plain description |
|-------|-------------------|
| **Left** | Sliders and settings that control what the galaxy looks like — number of stars, how many spiral arms, how connected they are, etc. Split into two tabs: **Generation** (shape and structure) and **Appearance** (colors and sizes). |
| **Centre** | The preview — a picture of the galaxy that updates when you click Preview. You can zoom in with the scroll wheel, pan by right-click-dragging, and click on a star to see its details. |
| **Right** | Information about the stars — population levels, who is in charge of what. Split into four tabs: **Attributes** (set population/admin rules), **Inspector** (see/edit a selected star), **Search** (find stars by property), **Filter** (show or hide specific groups of stars). |

**Typical first run:**

1. Open the app: `python galaxy_gui.py`
2. Leave everything on default and click **Generate** (bottom bar).
3. Click **Preview** — you will see a galaxy appear.
4. Click **Apply Attributes** (Attributes tab, right panel) to assign populations and governments.
5. Switch "Color by" (Appearance tab, left panel) to `hierarchy` and click **Preview** again to see political boundaries in color.
6. Click any star — the Inspector tab will fill with its details.
7. Use **Export PNG…** to save an image.

**Common questions:**

- *Why are some stars dark/unlit?*  They have no population (isolated nodes).
- *Why are some regions all one color in hierarchy mode?*  Those stars all report to the same Level-1 capital — they are one administrative empire.
- *What are L-ways?*  The lines between stars — hyperspace travel routes in the Emperum setting.
- *How do I get more connections?*  Increase "Avg degree" or "Max lane length" in the Generation tab.

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
| **Left** | Two tabs — Generation and Appearance |
| **Centre** | Embedded Matplotlib preview with zoom / pan |
| **Right** | Four tabs — Attributes, Inspector, Search, Filter |

The action bar runs along the bottom: **Generate**, **Preview**,
**Export PNG…**, **Export SVG…**, and a status line.

---

## Typical workflow

1. **Tune parameters** in the Generation tab (left panel).
2. Click **Generate** — builds the spatial graph and writes all output files.
3. Click **Preview** — renders the galaxy in the centre panel.
4. Use the **Attributes tab** (right) to adjust population / admin settings,
   then click **Apply Attributes** to re-run worldbuilding without rebuilding
   the spatial layout.
5. Switch `Color by` to `hierarchy` (Appearance tab) and click **Preview**
   again — no regeneration needed, just a visual redraw.
6. Click a node to inspect it in the **Inspector tab**.
7. Use the **Filter tab** to focus on a subset of nodes (see below).
8. **Export PNG…** or **Export SVG…** to save the current view.

---

## Left panel — Generation tab

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

> **Connectivity guarantee:** After the random edge selection, the generator
> automatically adds the shortest possible bridge edges to ensure that all
> connected nodes form a single network — no isolated islands.

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

---

## Left panel — Appearance tab

### Node Appearance

| Control | Default | Description |
|---------|---------|-------------|
| Color by | `arm_dist` | Data driving node color — see table below |
| Node size | 1.5 | Scatter marker size |
| Flat color | `#aaccff` | Used when `Color by = none` |
| Gradient low | `#ffe8c0` | Color at the **low** end of the gradient (low data value) |
| Gradient high | `#0d000f` | Color at the **high** end of the gradient (high data value) |

The **Gradient low** and **Gradient high** colors apply to **all** gradient
color modes — arm_dist, r, pop, admin_lvl, and admin_dist.  Whatever colors you
pick here will be used for any of those modes.

**Color by options:**

| Value | Description |
|-------|-------------|
| `arm_dist` | Distance to nearest spiral arm — arm nodes bright, inter-arm dark |
| `r` | Galactic radius — inner bright, outer dark |
| `pop` | Population index — low → gradient low color, 100 → gradient high color |
| `admin_lvl` | Admin level — low → gradient low color, level 5 → gradient high color |
| `admin_dist` | Hop-count to nearest admin node — 0 → gradient low color, far → gradient high color |
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

## Right panel — Attributes tab

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

---

## Right panel — Inspector tab

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

> After manually changing `admin_lvl` values, use **Apply Attributes** (Attributes
> tab) to recompute hierarchy assignments across all nodes.

---

## Right panel — Search tab

Search any column in the node table using a chosen attribute, operator, and
value.  Matching nodes are highlighted as green rings on the preview.

| Attribute options | `name`, `uid`, `pop`, `admin_lvl`, `admin_dist`, `hierarchy`, `id` |
|---|---|
| Operator options | `contains`, `=`, `>`, `<`, `>=`, `<=` |

The first matched node is auto-selected in the inspector.

---

## Right panel — Filter tab (reduced view)

The filter lets you show or hide a subset of nodes — and automatically hides
their connected edges too.  Useful for decluttering the view when you only want
to see, for example, highly populated systems or a specific administrative tier.

**How to use:**

1. Check **Enable reduced view filter**.
2. Enable one or more condition rows (check the checkbox at the left of each row).
3. For each enabled row, pick an attribute (e.g. `pop`), an operator (e.g. `>`), and a value (e.g. `50`).
4. Choose whether multiple conditions combine with **AND** (all must match) or **OR** (any must match).
5. Choose the mode:
   - **Show only matching** — only matching nodes are drawn; everything else vanishes.
   - **Hide matching** — matching nodes are hidden; everything else remains.
6. Click **Apply Filter** (or just click **Preview**) to redraw.

> The filter title in the preview shows how many nodes are currently hidden.

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

To visualize: select `hierarchy` from the **Color by** dropdown (Appearance
tab) and click **Preview**.

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
