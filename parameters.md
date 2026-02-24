# Emperum Galaxy Generator — Plain-Language Parameter Guide

This file explains every slider and setting in the GUI in simple terms.
No programming knowledge required.

---

## How to use this guide

Open the GUI with `python galaxy_gui.py`.  The left panel contains all the
settings, organised into collapsible sections (click a section header to
hide or show it).

**Important:** changing a setting in the **top five sections** (Spatial,
Edges, Spiral Arms, Radial Density, Sampling & Reproducibility) requires
you to press **Generate** again before the change takes effect.  Settings
in the bottom two sections (Node Appearance, Edge Appearance) only affect
how the map looks — press **Preview** to see them without regenerating.

---

## Spatial

These settings control the overall size and shape of the galaxy.

### Star systems (n_nodes) — default 5 000

How many star systems (dots) appear on the map.

- **Higher** → more stars, denser map, slower to generate.
- **Lower** → fewer stars, quicker to generate, useful for testing.

### Disk radius (r_disk) — default 100

The size of the whole galaxy.  Think of it as the outer edge: no star
will be placed beyond this boundary.  The actual unit doesn't matter —
it's just the scale everything else is measured against.

- Changing this alone makes the galaxy larger or smaller without changing
  its shape.

### Core radius (r_core) — default 15

The size of the empty zone at the very centre.  No stars are placed
inside it, and no hyperspace lanes may cross through it.

- **Larger** → a bigger empty "hole" in the middle.
- **Smaller** → stars crowd closer to the centre.

---

## Edges (Hyperspace Lanes)

These settings control the lines that connect star systems to each other.

### Max lane length (l_max) — default 9

The longest a single hyperspace lane can be.  Two stars that are farther
apart than this will never be directly connected.

- **Higher** → long-range connections are possible; the network looks more
  web-like.
- **Lower** → only nearby stars connect; you get tight local clusters with
  fewer bridges between them.

### Avg degree (target_degree) — default 4

Roughly how many lanes each star has on average.  This controls the
total number of lanes drawn across the whole map.

- **Higher** → busier, more connected map; most stars have several routes
  to choose from.
- **Lower** → sparser map; some stars may only have one or two connections.

### Connection chance — default 1.0 (= 100 %)

The probability that any given star is allowed to have lanes at all.
Stars that "fail" this check become completely isolated — no lanes in or
out, no matter how close their neighbours are.

- **1.0** → every star can connect (default).
- **0.5** → roughly half the stars will be isolated lone systems.
- **0.0** → no lanes at all (useful for inspecting node placement only).

---

## Spiral Arms

These settings shape the arms that give the galaxy its spiral structure.

### Number of arms (n_arms) — default 4

How many spiral arms the galaxy has.  Arms are spaced evenly around the
centre.

- **2** → a classic two-arm barred spiral.
- **4** → a richer, busier look (default).
- **6–8** → very many arms; starts to look more like a ring than a spiral
  at high values.

### Arm tightness (arm_b) — default 0.35

Controls how quickly the arms wind outward.

- **Lower (e.g. 0.15)** → arms coil tightly; they wrap around the centre
  many times before reaching the edge.
- **Higher (e.g. 0.6)** → arms unwind quickly and sweep outward in broad,
  open curves.

### Arm width σ (arm_sigma) — default 5.5

How wide or narrow each arm is.

- **Lower (e.g. 2)** → thin, sharply defined arms with clear empty space
  between them.
- **Higher (e.g. 12)** → wide, fluffy arms that blend together; the spiral
  pattern becomes subtler.

### Inter-arm density (arm_base) — default 0.15

How many stars appear in the space *between* the arms.

- **0.0** → the space between arms is completely empty; only the arms
  themselves have stars.
- **0.5** → a noticeable number of stars scattered between the arms.
- **1.0** → stars are spread uniformly everywhere; the spiral structure
  disappears.

---

## Radial Density

### Scale length (r_scale) — default 38

Controls how quickly star density falls off as you move away from the
centre.

- **Lower (e.g. 20)** → stars pack tightly near the core; the outer disk
  is quite sparse.
- **Higher (e.g. 70)** → stars are spread more evenly across the whole disk,
  including the outer edge.
- A value somewhere between the core radius and the disk radius usually
  looks most natural.

---

## Sampling & Reproducibility

### Boost multiplier (boost) — default 6

An internal tuning value.  You almost never need to touch this.

The generator places stars by random trial-and-error.  In some
configurations (very narrow arms, very little inter-arm density) it can
struggle to place enough stars within a reasonable number of tries.
Raising `boost` (e.g. to 10–15) helps if you see a warning message
about the generator failing to reach the requested number of stars.

### Random seed — default 7

A starting number for the random generator.  The same seed always
produces the exact same galaxy.

- Change this to get a completely different galaxy with otherwise
  identical settings.
- Write down the seed of any galaxy you like — you can recreate it
  exactly at any time.

### Output directory — default `output`

The folder where the generated files (`nodes.csv`, `edges.csv`,
`params.json`) are saved.  Will be created if it doesn't exist.
Click the **…** button to browse for a folder.

### Export GEXF checkbox — default on

Also saves a `graph.gexf` file, which can be opened directly in
[Gephi](https://gephi.org/) for interactive network exploration.
Uncheck this if you don't use Gephi or if the `networkx` library isn't
installed.

---

## Node Appearance

Changes here only affect the *look* of the map — press **Preview** to
see the result without needing to re-generate.

### Colour by — default `arm_dist`

What determines the colour of each star dot.

| Option | What it shows |
|---|---|
| `arm_dist` | Stars close to a spiral arm are coloured with the **Gradient low** colour; stars far from any arm get the **Gradient high** colour.  This makes the arms visually pop. |
| `r` | Stars near the centre are coloured with **Gradient low**; stars near the outer edge get **Gradient high**. |
| `none` | Every star gets the same **Flat colour**, ignoring position entirely. |

### Node size — default 1.5

How big each dot is.

- **Larger** → easier to see individual stars; starts to look crowded at
  high node counts.
- **Smaller** → finer, more detailed texture; better for large galaxies.

### Flat colour (mode: none) — default `#aaccff` (pale blue)

The colour used for every star when **Colour by** is set to `none`.
Has no effect in the other two modes.

### Gradient low (near arm) — default `#ffe8c0` (warm white)

The colour given to stars that score *low* on the chosen measure:
- In `arm_dist` mode: stars that are very close to an arm centreline.
- In `r` mode: stars that are close to the galactic centre.

### Gradient high (inter-arm) — default `#0d000f` (near-black)

The colour given to stars that score *high* on the chosen measure:
- In `arm_dist` mode: stars far from any arm (deep inter-arm space).
- In `r` mode: stars near the outer edge of the disk.

> **Tip:** swap Gradient low and Gradient high to invert the colour
> direction — e.g. bright outer rim instead of bright arms.

---

## Edge Appearance

Like Node Appearance, changes here only require **Preview** — no need
to regenerate.

### Hide edges checkbox — default off

Tick this to skip drawing the hyperspace lanes entirely.  The map will
show only star dots.

Strongly recommended if your node count is above roughly 1 000 — drawing
tens of thousands of lines is slow and can make the preview unreadable
anyway.

### Edge colour — default `#2244aa` (dark blue)

The colour of all hyperspace lane lines.

### Edge width — default 0.4

How thick the lane lines are.

- **Higher** → bolder, more visible lanes.
- **Lower** → finer lines; less visual clutter.

### Edge opacity — default 0.35

How transparent the lane lines are (`0` = invisible, `1` = fully solid).

- **Lower** → lanes are ghostly; the star layer shows through clearly.
- **Higher** → lanes are bold and prominent.

---

## Zoom and navigation in the Preview

| Action | Effect |
|---|---|
| **Scroll wheel up** | Zoom in, centred on the cursor |
| **Scroll wheel down** | Zoom out |
| **Toolbar — pan hand** | Click and drag to move around |
| **Toolbar — zoom box** | Draw a rectangle to zoom into that area |
| **Toolbar — home button** | Reset to the full galaxy view |
