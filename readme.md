# Emperum Galaxy Generator

Procedural generator for a face-on galactic disk with spiral arms and
hyperspace-lane edges, for the Emperum worldbuilding project.

Outputs `nodes.csv`, `edges.csv`, and (optionally) `graph.gexf` for use
in Gephi or any graph tool.

---

## GUI (recommended)

```bash
python galaxy_gui.py
```

This opens a desktop window with all parameters visible in a scrollable
left panel.  The workflow inside the GUI is:

1. **Adjust parameters** — every setting has a slider + number input or a
   colour picker; sections can be collapsed to reduce clutter.
2. **Generate** — builds the galaxy and writes `nodes.csv`, `edges.csv`,
   and `params.json` to the chosen output directory.
3. **Preview** — renders the galaxy in the embedded canvas.  You can
   re-render with different visual settings (colours, edge opacity, etc.)
   without re-generating.
4. **Export PNG… / Export SVG…** — saves a fresh render via a file-save
   dialog.

The GUI requires tkinter, which is bundled with the standard Python
installer.  On Ubuntu/Debian it can be installed with:

```bash
sudo apt-get install python3-tk
```

---

## Command-line quick start

```bash
# Generate with all defaults (5 000 nodes, 4 spiral arms)
python run_generate.py

# Visualise the result
python plot_debug.py

# Save the plot to a PNG
python plot_debug.py --save galaxy.png

# Save the plot as a scalable SVG
python plot_debug.py --svg galaxy.svg
```

---

## `run_generate.py` – generation parameters

### Spatial

| Argument | Default | Description |
|---|---|---|
| `--n_nodes` | `5000` | Total number of star-system nodes. More nodes = a denser, richer map but slower generation and larger output files. |
| `--r_disk` | `100.0` | Outer disk boundary radius (arbitrary units — think of it as the galactic edge). No nodes are placed beyond this radius. |
| `--r_core` | `15.0` | Radius of the forbidden galactic core. No nodes are placed inside it, and no hyperspace lane may pass through it. Increase to create a larger empty centre; decrease for a busier core region. |

### Edges (hyperspace lanes)

| Argument | Default | Description |
|---|---|---|
| `--l_max` | `9.0` | Hard maximum length for a single hyperspace lane. Any pair of nodes farther apart than this is **never** connected, regardless of other settings. Increase to allow longer routes and a more connected graph; decrease to force shorter hops and more localised clusters. Must be large enough relative to local node spacing for lanes to exist at all. |
| `--target_degree` | `4.0` | **Edge budget control.** Sets the total number of edges drawn: `total_edges ≈ N × D / 2`. At the default of 4 with 5 000 nodes that is ~10 000 edges. The generator finds all geometrically valid candidate pairs (within `l_max`, not crossing the core), shuffles them, and takes the first `total_edges` at random. This means the **average** node degree will be close to D, but individual nodes vary — nodes in dense spiral arms tend to have more connections than inter-arm nodes because they have more neighbours within `l_max`. Raise to increase overall connectivity; lower for a sparser network. |
| `--node_connection_chance` | `1.0` | Fraction of nodes that are eligible to receive **any** lanes at all. Each node independently passes a random check at this probability; nodes that fail are left completely isolated (no lanes in or out). `1.0` = every node can connect (default). `0.4` = roughly 40 % of nodes connect, 60 % are isolated lone stars. Does not change how connected eligible nodes are — only whether a node participates at all. |

### Spiral arms

| Argument | Default | Description |
|---|---|---|
| `--n_arms` | `4` | Number of logarithmic spiral arms evenly distributed around the disk. Typical real spirals have 2–4 arms. |
| `--arm_b` | `0.35` | Controls how tightly the arms are wound. The arm follows `r = r_arm_start × exp(b × φ)`. Smaller values (e.g. `0.2`) create tightly coiled arms that make many revolutions; larger values (e.g. `0.6`) produce wide-open arms that unwind quickly. |
| `--arm_sigma` | `5.5` | Gaussian half-width of each arm (in the same units as `r_disk`). Controls how broadly nodes cluster around the arm centreline. Small values (e.g. `2`) create narrow, sharply defined arms; large values (e.g. `10`) make wide, diffuse arms that blend into the inter-arm space. |
| `--arm_base` | `0.15` | Baseline node-placement probability anywhere in the disk, even far from all arms. Acts as an inter-arm density floor. `0.0` would leave the space between arms completely empty; `1.0` would spread nodes uniformly and erase the spiral structure. Default `0.15` gives sparse but noticeable inter-arm scatter. |

### Radial density

| Argument | Default | Description |
|---|---|---|
| `--r_scale` | `38.0` | Exponential scale length for how quickly node density drops off with radius: `P_radial(r) = exp(−r / r_scale)`. Smaller values (e.g. `20`) pack nodes tightly near the galactic centre and leave the outer disk sparse. Larger values (e.g. `60`) spread nodes more evenly across the full disk. Should generally stay in the range `[r_core, r_disk]`. |

### Sampling tuning

| Argument | Default | Description |
|---|---|---|
| `--boost` | `6.0` | Internal acceptance-probability multiplier used by the rejection sampler. The raw probability `P_radial × P_arm` can be very small in parts of the parameter space, so `boost` scales it up into a practical acceptance-rate range. You normally do not need to touch this. Increase it (e.g. to `10–15`) only if generation prints a warning about failing to reach `n_nodes` within the iteration limit — which can happen if you set very tight `arm_sigma` or very small `arm_base`. |

### Reproducibility & output

| Argument | Default | Description |
|---|---|---|
| `--seed` | `7` | Random seed. Any integer gives a different but fully reproducible galaxy. |
| `--out_dir` | `output` | Directory to write output files into (created if absent). |
| `--no_gexf` | flag | Skip GEXF export. Use this when `networkx` is not installed, or to speed up generation when you only need the CSV files. |

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
| `--out_dir` | `output` | Directory containing `nodes.csv` and `edges.csv`. Must match `--out_dir` used during generation. |
| `--save` | _(none)_ | Save the figure to a file instead of opening an interactive window. The format is inferred from the extension: `galaxy.png`, `galaxy.pdf`, etc. |
| `--svg` | _(none)_ | Save as SVG (scalable vector format). Pass a filename (`--svg galaxy.svg`) or omit the filename to default to `galaxy.svg`. SVG scales losslessly to any size and is suitable for import into Inkscape, Illustrator, or web pages. Takes precedence over `--save` when both are given. |

### Node appearance

| Argument | Default | Description |
|---|---|---|
| `--color_by` | `arm_dist` | What data drives the node colour gradient. `arm_dist` colours each node by its approximate distance to the nearest spiral arm centreline — arm stars are bright, inter-arm stars are dark. `r` colours by galactic radius — inner stars are one colour, outer stars the other. `none` applies a single flat colour to all nodes (see `--node_color`). |
| `--node_size` | `1.5` | Scatter marker size in matplotlib `s` units. Increase for larger, more visible dots; decrease for finer detail at high node counts. |
| `--node_color` | `#aaccff` | Flat colour used for every node when `--color_by none`. Ignored otherwise. |
| `--gradient_low_color` | `#ffe8c0` | Colour at the **low** end of the gradient. For `arm_dist` this is the colour of nodes **closest** to an arm; for `r` it is the **innermost** nodes. Default warm-white makes spiral arms visually prominent. |
| `--gradient_high_color` | `#0d000f` | Colour at the **high** end of the gradient. For `arm_dist` this is the colour of nodes **farthest** from any arm; for `r` it is the **outermost** nodes. Default near-black lets distant nodes fade into the background. |

> **Gradient direction note:** the defaults intentionally put bright at low values
> (arm stars) and dark at high values (inter-arm stars), which is the opposite of
> a standard heatmap.  To reverse this, swap the two colours:
> `--gradient_low_color "#0d000f" --gradient_high_color "#ffe8c0"`

### Edge appearance

| Argument | Default | Description |
|---|---|---|
| `--no_edges` | flag | Skip drawing edges entirely. Strongly recommended for N > 1 000 — drawing tens of thousands of line segments is slow. Use to inspect node distribution without waiting for edge rendering. |
| `--edge_color` | `#2244aa` | Colour of all edge lines (any matplotlib colour string or hex code). |
| `--edge_width` | `0.4` | Edge line width in points. Increase for bolder lanes; decrease to reduce clutter at high edge counts. |
| `--edge_alpha` | `0.35` | Edge opacity (`0` = fully invisible, `1` = fully opaque). Lower values let the node layer show through a dense edge network. |

### Spatial parameters (auto-loaded from `params.json`)

`run_generate.py` writes a `params.json` file to `--out_dir` alongside the
CSVs. `plot_debug.py` reads it automatically, so the overlay graphics (disk
circle, core circle, arm centrelines) always match the actual generation run
without any manual re-entry. You only need to pass these arguments explicitly
if you want to **override** what is in the saved file.

| Argument | Source |
|---|---|
| `--r_disk` | auto from `params.json` (fallback `100.0`) |
| `--r_core` | auto from `params.json` (fallback `15.0`) |
| `--n_arms` | auto from `params.json` (fallback `4`) |
| `--arm_b` | auto from `params.json` (fallback `0.35`) |
| `--r_arm_start` | auto from `params.json` (fallback `3.0`) |

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
    --svg my_run/galaxy.svg
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
