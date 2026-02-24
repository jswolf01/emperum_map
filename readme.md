# Emperum Galaxy Generator

Procedural generator for a face-on galactic disk with spiral arms and
hyperspace-lane edges, for the Emperum worldbuilding project.

Outputs `nodes.csv`, `edges.csv`, and (optionally) `graph.gexf` for use
in Gephi or any graph tool.

---

## Quick start

```bash
# Generate with all defaults (5 000 nodes, 4 spiral arms)
python run_generate.py

# Visualise the result
python plot_debug.py

# Save the plot to a PNG
python plot_debug.py --save galaxy.png
```

---

## `run_generate.py` – generation parameters

### Spatial

| Argument | Default | Description |
|---|---|---|
| `--n_nodes` | `5000` | Number of star-system nodes to generate. |
| `--r_disk` | `100.0` | Outer disk boundary radius (arbitrary units). |
| `--r_core` | `15.0` | Forbidden core radius. No nodes are placed inside it and no edge may cross it. |

### Edges

| Argument | Default | Description |
|---|---|---|
| `--l_max` | `9.0` | Maximum edge (L-way hyperspace lane) length. Pairs farther apart than this are never connected. |
| `--target_degree` | `4.0` | Approximate target average node degree. Controls how many edges are drawn in total (`≈ N × D / 2`). |
| `--node_connection_chance` | `1.0` | Probability `[0, 1]` that any given node is eligible to have edge connections at all. Ineligible nodes remain fully isolated. `1.0` = all nodes may connect; `0.3` = roughly 30 % of nodes connect, 70 % are isolated. Edge density among connected nodes is driven purely by local node density (spiral arm placement has no additional effect). |

### Spiral arms

| Argument | Default | Description |
|---|---|---|
| `--n_arms` | `4` | Number of logarithmic spiral arms. |
| `--arm_b` | `0.35` | Spiral tightness. Larger values produce more open, loosely wound arms. |
| `--arm_sigma` | `5.5` | Gaussian half-width of each arm (same units as radii). Controls how tightly nodes cluster around the arm centreline. |
| `--arm_base` | `0.15` | Baseline acceptance probability far from any arm `(0, 1)`. Higher values fill the inter-arm space more evenly. |

### Radial density

| Argument | Default | Description |
|---|---|---|
| `--r_scale` | `38.0` | Exponential scale length for the radial density falloff. Smaller values concentrate nodes near the centre. |

### Sampling tuning

| Argument | Default | Description |
|---|---|---|
| `--boost` | `6.0` | Acceptance-probability multiplier for the rejection sampler. Increase if sampling is very slow or fails to reach `n_nodes`. |

### Reproducibility & output

| Argument | Default | Description |
|---|---|---|
| `--seed` | `7` | Random seed. Use any integer for a different but reproducible galaxy. |
| `--out_dir` | `output` | Directory to write output files (created if absent). |
| `--no_gexf` | flag | Skip GEXF export. Useful when networkx is not installed. |

### Full example

```bash
python run_generate.py \
    --n_nodes 8000 \
    --r_disk 100 \
    --r_core 15 \
    --l_max 9 \
    --target_degree 4 \
    --node_connection_chance 0.4 \
    --n_arms 4 \
    --arm_b 0.35 \
    --arm_sigma 5.5 \
    --arm_base 0.15 \
    --r_scale 38 \
    --seed 42 \
    --out_dir my_run
```

---

## `plot_debug.py` – visualisation parameters

### Input & output

| Argument | Default | Description |
|---|---|---|
| `--out_dir` | `output` | Directory containing `nodes.csv` and `edges.csv`. |
| `--save` | _(none)_ | Save the figure to a file (`galaxy.png`, `galaxy.pdf`, etc.) instead of opening an interactive window. |

### Node appearance

| Argument | Default | Description |
|---|---|---|
| `--color_by` | `arm_dist` | Node colouring scheme. `arm_dist` = gradient by distance to nearest spiral arm centreline; `r` = gradient by galactic radius; `none` = uniform flat colour. |
| `--node_size` | `1.5` | Scatter marker size (matplotlib `s` units). |
| `--node_color` | `#aaccff` | Flat colour used when `--color_by none`. |
| `--gradient_low_color` | `#ffe8c0` | Gradient colour at **low** data values. For `arm_dist` this applies to nodes closest to arm centrelines; for `r` it applies to the innermost nodes. Default (warm bright) makes arm stars visually prominent. |
| `--gradient_high_color` | `#0d000f` | Gradient colour at **high** data values. For `arm_dist` this applies to inter-arm nodes; for `r` it applies to outer-disk nodes. Default (near-black) lets distant nodes recede into the background. |

> **Gradient direction note:** the defaults are intentionally inverted from
> the classic `plasma` colourmap — arm stars (low `arm_dist`) are rendered
> bright and visible, while inter-arm stars fade to near-black.  To restore
> the old look, swap the two values:
> `--gradient_low_color "#0d000f" --gradient_high_color "#ffe8c0"`

### Edge appearance

| Argument | Default | Description |
|---|---|---|
| `--no_edges` | flag | Skip drawing edges entirely. Much faster for large graphs (N > 1 000). |
| `--edge_color` | `#2244aa` | Edge line colour (any matplotlib colour string or hex). |
| `--edge_width` | `0.4` | Edge line width in points. |
| `--edge_alpha` | `0.35` | Edge opacity (`0` = invisible, `1` = fully opaque). |

### Spatial parameters (must match generation)

These control the overlay graphics (disk circle, core circle, arm
centrelines) and must be set to the same values used in `run_generate.py`.

| Argument | Default |
|---|---|
| `--r_disk` | `100.0` |
| `--r_core` | `15.0` |
| `--n_arms` | `4` |
| `--arm_b` | `0.35` |
| `--r_arm_start` | `3.0` |

### Full example

```bash
python plot_debug.py \
    --out_dir my_run \
    --color_by arm_dist \
    --node_size 2.0 \
    --gradient_low_color "#ffe8c0" \
    --gradient_high_color "#0d000f" \
    --edge_color "#3366cc" \
    --edge_width 0.5 \
    --edge_alpha 0.4 \
    --save my_run/galaxy.png
```

---

## Output files

| File | Description |
|---|---|
| `output/nodes.csv` | One row per star system: `id, x, y, r, theta, arm_dist` |
| `output/edges.csv` | One row per hyperspace lane: `source, target, length, weight` |
| `output/graph.gexf` | Full graph for Gephi (requires `networkx`; skip with `--no_gexf`) |

---

## Dependencies

```
numpy
pandas
scipy
matplotlib
networkx  # optional – only needed for GEXF export
```

Install with:

```bash
pip install numpy pandas scipy matplotlib networkx
```
