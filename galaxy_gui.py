"""
galaxy_gui.py
=============
Tkinter GUI front-end for the Emperum galaxy generator.

Layout
------
Left panel   – tabbed:
                 • Generation   – spatial, edges, arms, density, sampling.
                 • Appearance   – output directory, node/edge colours and sizes.
Centre panel – embedded matplotlib preview with zoom/pan and scroll-wheel zoom.
Right panel  – tabbed:
                 • Attributes      – population and admin-level parameters.
                 • Inspect & Search – node inspector (click to select) + search.
                 • Filter View     – reduced-view filter: hide/show nodes by
                                    attribute conditions with AND/OR logic.

Usage
-----
    python galaxy_gui.py

Dependencies
------------
Same as the core generator (numpy, pandas, scipy, matplotlib) plus tkinter,
which is bundled with the standard Python installer.  On Ubuntu/Debian:
    sudo apt-get install python3-tk
"""

from __future__ import annotations

import json
import os
import threading
import types
from typing import Optional

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import tkinter as tk
from tkinter import ttk, colorchooser, filedialog, messagebox

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from galaxygen import GalaxyConfig, GalaxyGenerator, compute_admin_dist
from plot_debug import draw_galaxy


# ---------------------------------------------------------------------------
# Reusable compound widgets
# ---------------------------------------------------------------------------

class SliderEntry(ttk.Frame):
    """Linked horizontal scale + spinbox for a numeric parameter."""

    def __init__(
        self,
        parent,
        label: str,
        var: tk.Variable,
        lo: float,
        hi: float,
        step: float = 1.0,
        label_width: int = 24,
        spin_width: int = 7,
        **kw,
    ):
        super().__init__(parent, **kw)
        self._var  = var
        self._step = step
        self._lo   = lo
        self._hi   = hi
        self._busy = False

        ttk.Label(self, text=label, width=label_width, anchor="w").grid(
            row=0, column=0, sticky="w", padx=(4, 2), pady=1,
        )
        self._scale = ttk.Scale(
            self, orient="horizontal", length=130,
            from_=lo, to=hi, variable=var,
            command=self._on_scale,
        )
        self._scale.grid(row=0, column=1, padx=4)
        self._spin = ttk.Spinbox(
            self, from_=lo, to=hi, increment=step,
            textvariable=var, width=spin_width,
        )
        self._spin.grid(row=0, column=2, padx=(2, 4))
        self._spin.bind("<Return>",   self._clamp)
        self._spin.bind("<FocusOut>", self._clamp)

    def _on_scale(self, _val: str) -> None:
        if self._busy:
            return
        try:
            raw = float(_val)
        except ValueError:
            return
        snapped = round(raw / self._step) * self._step
        snapped = round(snapped, 10)
        if abs(raw - snapped) > 1e-9:
            self._busy = True
            self._var.set(snapped)
            self._busy = False

    def _clamp(self, _evt=None) -> None:
        try:
            val = float(self._spin.get())
        except ValueError:
            val = self._lo
        val = max(self._lo, min(self._hi, val))
        self._var.set(round(val / self._step) * self._step)


# ---------------------------------------------------------------------------

class ColorEntry(ttk.Frame):
    """Colour swatch + hex entry + HSV colour-picker button."""

    def __init__(self, parent, label: str, var: tk.StringVar,
                 label_width: int = 24, **kw):
        super().__init__(parent, **kw)
        self._var = var

        ttk.Label(self, text=label, width=label_width, anchor="w").grid(
            row=0, column=0, sticky="w", padx=(4, 2), pady=1,
        )
        self._swatch = tk.Label(self, width=3, relief="sunken", cursor="hand2")
        self._swatch.grid(row=0, column=1, padx=(2, 2))
        self._swatch.bind("<Button-1>", self._open_picker)

        self._entry = ttk.Entry(self, textvariable=var, width=10)
        self._entry.grid(row=0, column=2, padx=2)
        self._entry.bind("<Return>",   self._refresh_swatch)
        self._entry.bind("<FocusOut>", self._refresh_swatch)

        ttk.Button(self, text="Pick…", width=6,
                   command=self._open_picker).grid(row=0, column=3, padx=(2, 4))

        var.trace_add("write", self._refresh_swatch)
        self._refresh_swatch()

    def _refresh_swatch(self, *_) -> None:
        val = self._var.get().strip()
        try:
            self._swatch.configure(bg=val)
        except tk.TclError:
            self._swatch.configure(bg="#888888")

    def _open_picker(self, _evt=None) -> None:
        current = self._var.get()
        try:
            _rgb, hexval = colorchooser.askcolor(
                color=current, title="Choose colour", parent=self)
        except Exception:
            return
        if hexval:
            self._var.set(hexval.lower())


# ---------------------------------------------------------------------------

class Section(ttk.Frame):
    """Collapsible parameter section with a toggle-button header."""

    def __init__(self, parent, title: str, start_open: bool = True, **kw):
        super().__init__(parent, **kw)
        self._open  = start_open
        self._title = title

        self._btn = ttk.Button(self, text=f"{'▼' if start_open else '▶'}  {title}",
                               command=self._toggle)
        self._btn.pack(fill="x", padx=2, pady=(4, 0))
        ttk.Separator(self, orient="horizontal").pack(fill="x", padx=2)

        self._inner = ttk.Frame(self, padding=(2, 2, 2, 6))
        if start_open:
            self._inner.pack(fill="x", expand=True)

    @property
    def inner(self) -> ttk.Frame:
        return self._inner

    def _toggle(self) -> None:
        if self._open:
            self._inner.pack_forget()
            self._btn.configure(text=f"▶  {self._title}")
        else:
            self._inner.pack(fill="x", expand=True)
            self._btn.configure(text=f"▼  {self._title}")
        self._open = not self._open


# ---------------------------------------------------------------------------
# Main GUI class
# ---------------------------------------------------------------------------

class GalaxyGUI:
    """Top-level GUI application window."""

    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        root.title("Emperum Galaxy Generator")
        root.minsize(1300, 760)

        self._current_fig: Optional[plt.Figure] = None
        self._canvas_widget: Optional[FigureCanvasTkAgg] = None
        self._toolbar_frame: Optional[ttk.Frame] = None
        self._worker: Optional[threading.Thread] = None

        # Live data for node interaction
        self._nodes_df: Optional[pd.DataFrame] = None
        self._edges_df: Optional[pd.DataFrame] = None
        self._node_tree: Optional[cKDTree] = None   # spatial KD-tree for click detection
        self._selected_node: Optional[int] = None   # currently selected node index
        self._found_nodes: list[int] = []            # search results

        self._build_vars()
        self._build_ui()

    # ── Variable definitions ──────────────────────────────────────────────

    def _build_vars(self) -> None:
        iv = tk.IntVar
        dv = tk.DoubleVar
        sv = tk.StringVar
        bv = tk.BooleanVar

        # ── Generation ──────────────────────────────────────────────────
        self.v_n_nodes                = iv(value=5_000)
        self.v_r_disk                 = dv(value=100.0)
        self.v_r_core                 = dv(value=15.0)
        self.v_l_max                  = dv(value=9.0)
        self.v_target_degree          = dv(value=4.0)
        self.v_node_connection_chance = dv(value=1.0)
        self.v_n_arms                 = iv(value=4)
        self.v_arm_b                  = dv(value=0.35)
        self.v_arm_sigma              = dv(value=5.5)
        self.v_arm_base               = dv(value=0.15)
        self.v_r_scale                = dv(value=38.0)
        self.v_boost                  = dv(value=6.0)
        self.v_seed                   = iv(value=7)
        self.v_out_dir                = sv(value="output")
        self.v_write_gexf             = bv(value=True)

        # ── Visualisation ────────────────────────────────────────────────
        self.v_color_by               = sv(value="arm_dist")
        self.v_node_size              = dv(value=1.5)
        self.v_node_color             = sv(value="#aaccff")
        self.v_grad_low               = sv(value="#ffe8c0")
        self.v_grad_high              = sv(value="#0d000f")
        self.v_no_edges               = bv(value=False)
        self.v_edge_color             = sv(value="#2244aa")
        self.v_edge_width             = dv(value=0.4)
        self.v_edge_alpha             = dv(value=0.35)

        # ── Population ───────────────────────────────────────────────────
        self.v_pop_core_dispersal = dv(value=1.0)
        self.v_pop_dispersal      = dv(value=1.0)

        # ── Admin levels – exact counts (-1 = auto) ──────────────────────
        self.v_admin_count_1 = iv(value=-1)
        self.v_admin_count_2 = iv(value=-1)
        self.v_admin_count_3 = iv(value=-1)
        self.v_admin_count_4 = iv(value=-1)
        self.v_admin_count_5 = iv(value=-1)

        # ── Admin levels – ratios ────────────────────────────────────────
        self.v_admin_ratio_21 = dv(value=2.0)
        self.v_admin_ratio_32 = dv(value=3.0)
        self.v_admin_ratio_43 = dv(value=4.0)
        self.v_admin_ratio_54 = dv(value=5.0)

        # ── Admin levels – spatial spacing ───────────────────────────────
        self.v_admin_coverage_scale = dv(value=1.0)
        self.v_admin_sep_scale      = dv(value=1.5)

        # ── Node inspector (editable fields for selected node) ────────────
        self.v_insp_uid       = sv(value="")
        self.v_insp_name      = sv(value="")
        self.v_insp_pop       = iv(value=0)
        self.v_insp_admin_lvl = iv(value=0)
        self.v_insp_admin_dist= sv(value="")
        self.v_insp_hierarchy = sv(value="")

        # ── Search ───────────────────────────────────────────────────────
        self.v_search_attr    = sv(value="name")
        self.v_search_op      = sv(value="contains")
        self.v_search_query   = sv(value="")
        self.v_search_result  = sv(value="")

        # ── Filter (reduced view) ─────────────────────────────────────────
        self.v_filter_enabled = bv(value=False)
        self.v_filter_logic   = sv(value="AND")   # AND | OR
        self.v_filter_action  = sv(value="hide_matching")  # hide_matching | show_only
        self.v_filter_status  = sv(value="")
        # Dynamic list of condition rows – each entry is a dict with keys:
        #   attr_var, op_var, val_var, val2_var, frame, to_label, val2_entry
        self._filter_rows: list[dict] = []
        self._filter_container: tk.Frame | None = None  # set when tab is built
        # ── Reduced-view filter ───────────────────────────────────────────
        self.v_filter_active = bv(value=False)
        self.v_filter_logic  = sv(value="AND")
        self.v_filter_mode   = sv(value="show_matching")
        self._filter_conds: list[dict] = [
            {
                "enabled": bv(value=False),
                "attr":    sv(value=""),
                "op":      sv(value="="),
                "val":     sv(value=""),
            }
            for _ in range(3)
        ]

    # ── UI construction ───────────────────────────────────────────────────

    def _build_ui(self) -> None:
        self._build_action_bar()

        paned = ttk.PanedWindow(self.root, orient="horizontal")
        paned.pack(fill="both", expand=True, padx=6, pady=(6, 0))

        # Left panel: generation parameters
        left_outer = ttk.Frame(paned, width=370)
        left_outer.pack_propagate(False)
        paned.add(left_outer, weight=0)

        # Centre panel: preview
        centre_frame = ttk.Frame(paned)
        paned.add(centre_frame, weight=1)

        # Right panel: worldbuilding attributes + inspector + search
        right_outer = ttk.Frame(paned, width=360)
        right_outer.pack_propagate(False)
        paned.add(right_outer, weight=0)

        self._build_param_panel(left_outer)
        self._build_preview_panel(centre_frame)
        self._build_attr_panel(right_outer)

    # ── Helpers ───────────────────────────────────────────────────────────

    def _make_scrollable(self, parent: ttk.Frame) -> ttk.Frame:
        """Create a vertically-scrollable pane inside *parent*.

        Returns the inner ``ttk.Frame`` into which content should be packed.
        """
        canvas = tk.Canvas(parent, highlightthickness=0, borderwidth=0)
        vscroll = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vscroll.set)
        vscroll.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        inner = ttk.Frame(canvas)
        win_id = canvas.create_window((0, 0), window=inner, anchor="nw")
        inner.bind("<Configure>",
                   lambda _e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>",
                    lambda e: canvas.itemconfigure(win_id, width=e.width))
        canvas.bind("<MouseWheel>",
                    lambda e: canvas.yview_scroll(int(-1 * e.delta / 120), "units"))
        canvas.bind("<Button-4>", lambda _e: canvas.yview_scroll(-1, "units"))
        canvas.bind("<Button-5>", lambda _e: canvas.yview_scroll(1, "units"))
        return inner

    @staticmethod
    def _make_scrollable_tab(tab_frame: ttk.Frame):
        """Create a scrollable canvas inside a tab frame; return the inner content Frame."""
        canvas = tk.Canvas(tab_frame, highlightthickness=0, borderwidth=0)
        vscroll = ttk.Scrollbar(tab_frame, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vscroll.set)
        vscroll.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        inner = ttk.Frame(canvas)
        win_id = canvas.create_window((0, 0), window=inner, anchor="nw")
        inner.bind("<Configure>",
                   lambda _e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>",
                    lambda e: canvas.itemconfigure(win_id, width=e.width))

        def _wheel(evt):
            if evt.delta:
                canvas.yview_scroll(int(-1 * evt.delta / 120), "units")
        canvas.bind("<MouseWheel>", _wheel)
        canvas.bind("<Button-4>", lambda _e: canvas.yview_scroll(-1, "units"))
        canvas.bind("<Button-5>", lambda _e: canvas.yview_scroll(1, "units"))

        return inner

    # ── Left parameter panel ──────────────────────────────────────────────

    def _build_param_panel(self, parent: ttk.Frame) -> None:
        """Left panel — two tabs: Generation parameters and Appearance settings."""
        nb = ttk.Notebook(parent)
        nb.pack(fill="both", expand=True)

        gen_tab = ttk.Frame(nb)
        nb.add(gen_tab, text="Generation")
        app_tab = ttk.Frame(nb)
        nb.add(app_tab, text="Appearance")

        self._build_generation_tab(gen_tab)
        self._build_appearance_tab(app_tab)

    def _build_generation_tab(self, parent: ttk.Frame) -> None:
        """Generation parameters: spatial, edges, arms, density, sampling."""
        inner = self._make_scrollable(parent)
        LW = 24

        # ── Generation tab ────────────────────────────────────────────

        sec = Section(inner, "Spatial")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        SliderEntry(s, "Star systems (n_nodes)",  self.v_n_nodes, 100, 20_000, 100, label_width=LW).pack(fill="x")
        SliderEntry(s, "Disk radius (r_disk)",    self.v_r_disk,  10.0, 500.0, 1.0, label_width=LW).pack(fill="x")
        SliderEntry(s, "Core radius (r_core)",    self.v_r_core,  0.0, 100.0, 0.5, label_width=LW).pack(fill="x")

        sec = Section(inner, "Edges (Hyperspace Lanes)")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        SliderEntry(s, "Max lane length (l_max)",    self.v_l_max,                  1.0, 50.0, 0.5,  label_width=LW).pack(fill="x")
        SliderEntry(s, "Avg degree (target_degree)", self.v_target_degree,           1.0, 20.0, 0.5,  label_width=LW).pack(fill="x")
        SliderEntry(s, "Connection chance",          self.v_node_connection_chance,  0.0,  1.0, 0.05, label_width=LW).pack(fill="x")

        sec = Section(inner, "Spiral Arms")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        SliderEntry(s, "Number of arms (n_arms)",     self.v_n_arms,    1,    8,    1,    label_width=LW).pack(fill="x")
        SliderEntry(s, "Arm tightness (arm_b)",       self.v_arm_b,     0.01, 2.0,  0.01, label_width=LW).pack(fill="x")
        SliderEntry(s, "Arm width σ (arm_sigma)",     self.v_arm_sigma, 0.5,  30.0, 0.5,  label_width=LW).pack(fill="x")
        SliderEntry(s, "Inter-arm density (arm_base)", self.v_arm_base, 0.0,   1.0, 0.01, label_width=LW).pack(fill="x")

        sec = Section(inner, "Radial Density")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        SliderEntry(s, "Scale length (r_scale)", self.v_r_scale, 5.0, 200.0, 1.0, label_width=LW).pack(fill="x")

        # ── Sampling & Seed ───────────────────────────────────────────
        sec = Section(inner, "Sampling & Reproducibility")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        SliderEntry(s, "Boost multiplier (boost)", self.v_boost, 1.0, 20.0, 0.5, label_width=LW).pack(fill="x")
        row = ttk.Frame(s)
        row.pack(fill="x", pady=1)
        ttk.Label(row, text="Random seed", width=LW, anchor="w").pack(side="left", padx=(4, 2))
        ttk.Spinbox(row, from_=0, to=99_999, increment=1,
                    textvariable=self.v_seed, width=8).pack(side="left")

    def _build_appearance_tab(self, parent: ttk.Frame) -> None:
        """Output directory + GEXF export + node/edge appearance."""
        inner = self._make_scrollable(parent)
        LW = 24

        # ── Output ────────────────────────────────────────────────────
        sec = Section(inner, "Output")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        row = ttk.Frame(s)
        row.pack(fill="x", pady=1)
        ttk.Label(row, text="Output directory", width=LW, anchor="w").pack(side="left", padx=(4, 2))
        ttk.Entry(row, textvariable=self.v_out_dir, width=12).pack(side="left")
        ttk.Button(row, text="…", width=3, command=self._browse_out_dir).pack(side="left", padx=2)
        row2 = ttk.Frame(s)
        row2.pack(fill="x", pady=1)
        ttk.Checkbutton(row2, text="Export GEXF (requires networkx)",
                        variable=self.v_write_gexf).pack(side="left", padx=(4, 2))

        # ── Appearance tab ────────────────────────────────────────────

        sec = Section(inner, "Node Appearance")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        row = ttk.Frame(s)
        row.pack(fill="x", pady=1)
        ttk.Label(row, text="Color by", width=LW, anchor="w").pack(side="left", padx=(4, 2))
        ttk.Combobox(
            row, textvariable=self.v_color_by,
            values=["arm_dist", "r", "pop", "admin_lvl", "admin_dist",
                    "hierarchy", "none"],
            state="readonly", width=11,
        ).pack(side="left")
        SliderEntry(s, "Node size", self.v_node_size, 0.1, 20.0, 0.1, label_width=LW).pack(fill="x")
        ColorEntry(s,  "Flat color  (mode: none)", self.v_node_color, label_width=LW).pack(fill="x")
        ColorEntry(s,  "Gradient low  (low value)",  self.v_grad_low,   label_width=LW).pack(fill="x")
        ColorEntry(s,  "Gradient high  (high value)", self.v_grad_high,  label_width=LW).pack(fill="x")
        ttk.Label(s, text=(
            "  Gradient colors apply to arm_dist and r modes.\n"
            "  Other modes use built-in semantic colormaps."
        ), foreground="#888888", wraplength=310, justify="left").pack(anchor="w", padx=6)

        # ── Edge Appearance ───────────────────────────────────────────
        sec = Section(inner, "Edge Appearance")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        ttk.Checkbutton(s, text="Hide edges  (much faster for N > 1 000)",
                        variable=self.v_no_edges).pack(anchor="w", padx=4, pady=2)
        ColorEntry(s, "Edge color",    self.v_edge_color, label_width=LW).pack(fill="x")
        SliderEntry(s, "Edge width",   self.v_edge_width, 0.1, 5.0, 0.1,  label_width=LW).pack(fill="x")
        SliderEntry(s, "Edge opacity", self.v_edge_alpha, 0.0, 1.0, 0.05, label_width=LW).pack(fill="x")

    # ── Right attribute panel ──────────────────────────────────────────────

    def _build_attr_panel(self, parent: ttk.Frame) -> None:
        """Tabbed right panel: Attributes | Inspect & Search | Filter."""
        nb = ttk.Notebook(parent)
        nb.pack(fill="both", expand=True)

        attr_frame   = ttk.Frame(nb)
        insp_frame   = ttk.Frame(nb)
        filter_frame = ttk.Frame(nb)

        nb.add(attr_frame,   text="Attributes")
        nb.add(insp_frame,   text="Inspect & Search")
        nb.add(filter_frame, text="Filter View")

        self._build_attributes_tab(attr_frame)
        self._build_inspect_tab(insp_frame)
        self._build_filter_tab(filter_frame)

    def _build_attributes_tab(self, parent: ttk.Frame) -> None:
        """Population, admin levels, and Apply Attributes button."""
        inner = self._make_scrollable(parent)
        RLW = 22

        # ── Attributes tab ────────────────────────────────────────────

        ttk.Button(inner, text="Apply Attributes",
                   command=self._on_apply_attrs).pack(fill="x", padx=6, pady=(6, 2))
        ttk.Label(inner,
                  text="Re-assigns pop/admin without regenerating spatial layout.",
                  foreground="#888888", wraplength=330,
                  justify="left").pack(padx=6, pady=(0, 4))

        sec = Section(inner, "Population")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        ttk.Label(s, text=(
            "pop=0 for isolated nodes.  Connected nodes follow a bell curve "
            "modified by galactic position and neighbourhood clustering."
        ), wraplength=320, justify="left", foreground="#aaaaaa").pack(padx=4, pady=(2, 6))
        SliderEntry(s, "Core dispersal", self.v_pop_core_dispersal,
                    0.0, 5.0, 0.1, label_width=RLW).pack(fill="x")
        ttk.Label(s, text="  0=uniform  |  1=moderate core-high  |  5=extreme core bias",
                  foreground="#888888").pack(anchor="w", padx=6)
        SliderEntry(s, "Clustering", self.v_pop_dispersal,
                    0.0, 5.0, 0.1, label_width=RLW).pack(fill="x")
        ttk.Label(s, text="  0=no clustering  |  1=moderate  |  5=strong neighbourhood",
                  foreground="#888888").pack(anchor="w", padx=6)

        sec = Section(inner, "Admin Levels")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        ttk.Label(s, text=(
            "Lvl 1 = top-tier capital  |  Lvl 5 = regional logistics centre.  "
            "Nodes selected by high pop + high connectivity.  "
            "Set counts to -1 to use ratios."
        ), wraplength=320, justify="left", foreground="#aaaaaa").pack(padx=4, pady=(2, 6))

        ttk.Label(s, text="── Exact counts  (−1 = use ratios) ──",
                  foreground="#cccccc").pack(anchor="w", padx=4, pady=(4, 0))
        for lvl, var in [(1, self.v_admin_count_1), (2, self.v_admin_count_2),
                         (3, self.v_admin_count_3), (4, self.v_admin_count_4),
                         (5, self.v_admin_count_5)]:
            row = ttk.Frame(s)
            row.pack(fill="x", pady=1)
            ttk.Label(row, text=f"  Lvl {lvl} count", width=RLW, anchor="w").pack(
                side="left", padx=(4, 2))
            ttk.Spinbox(row, from_=-1, to=10_000, increment=1,
                        textvariable=var, width=7).pack(side="left")

        ttk.Label(s, text="── Ratios  (higher ÷ lower level) ──",
                  foreground="#cccccc").pack(anchor="w", padx=4, pady=(8, 0))
        for label, var in [("Lvl 2 : Lvl 1 ratio", self.v_admin_ratio_21),
                            ("Lvl 3 : Lvl 2 ratio", self.v_admin_ratio_32),
                            ("Lvl 4 : Lvl 3 ratio", self.v_admin_ratio_43),
                            ("Lvl 5 : Lvl 4 ratio", self.v_admin_ratio_54)]:
            SliderEntry(s, label, var, 1.0, 20.0, 0.5, label_width=RLW).pack(fill="x")

        ttk.Label(s, text="── Spatial spacing ──",
                  foreground="#cccccc").pack(anchor="w", padx=4, pady=(8, 0))
        SliderEntry(s, "Coverage scale", self.v_admin_coverage_scale,
                    0.1, 5.0, 0.1, label_width=RLW).pack(fill="x")
        ttk.Label(s, text=(
            "  Scales the auto-computed coverage radius for each level.\n"
            "  Lower = denser centres; higher = sparser, larger regions."
        ), foreground="#888888", wraplength=320).pack(anchor="w", padx=6)
        SliderEntry(s, "Separation scale", self.v_admin_sep_scale,
                    0.5, 5.0, 0.1, label_width=RLW).pack(fill="x")
        ttk.Label(s, text=(
            "  Min-separation = coverage_radius × this.  "
            "  Lower = centres can cluster; higher = forced wide spread."
        ), foreground="#888888", wraplength=320).pack(anchor="w", padx=6)

    def _build_inspect_tab(self, parent: ttk.Frame) -> None:
        """Node inspector + search in a scrollable tab."""
        inner = self._make_scrollable(parent)
        RLW = 22

        # ── Node Inspector ────────────────────────────────────────────
        sec = Section(inner, "Node Inspector", start_open=True)
        sec.pack(fill="x", padx=4, pady=3)
        self._inspector_section = sec
        s = sec.inner

        ttk.Label(s, text="Click a node on the preview to inspect it.",
                  foreground="#888888", wraplength=320).pack(padx=4, pady=(2, 4))

        def _insp_row(parent, label, widget_factory, pady=1):
            row = ttk.Frame(parent)
            row.pack(fill="x", pady=pady)
            ttk.Label(row, text=label, width=RLW, anchor="w").pack(
                side="left", padx=(4, 2))
            w = widget_factory(row)
            w.pack(side="left")
            return w

        _insp_row(s, "System ID (uid)",
                  lambda p: ttk.Entry(p, textvariable=self.v_insp_uid,
                                      state="readonly", width=12))
        _insp_row(s, "Admin dist (hops)",
                  lambda p: ttk.Entry(p, textvariable=self.v_insp_admin_dist,
                                      state="readonly", width=8))
        _insp_row(s, "Hierarchy index",
                  lambda p: ttk.Entry(p, textvariable=self.v_insp_hierarchy,
                                      state="readonly", width=8))
        _insp_row(s, "Name",
                  lambda p: ttk.Entry(p, textvariable=self.v_insp_name, width=18))
        _insp_row(s, "Population (0-100)",
                  lambda p: ttk.Spinbox(p, from_=0, to=100, increment=1,
                                        textvariable=self.v_insp_pop, width=6))
        _insp_row(s, "Admin level (0-5)",
                  lambda p: ttk.Spinbox(p, from_=0, to=5, increment=1,
                                        textvariable=self.v_insp_admin_lvl, width=6))

        btn_row = ttk.Frame(s)
        btn_row.pack(fill="x", pady=(4, 2))
        ttk.Button(btn_row, text="Save Node Edits",
                   command=self._on_save_node).pack(side="left", padx=4)
        ttk.Button(btn_row, text="Recalc Admin Dist",
                   command=self._on_recalc_admin_dist).pack(side="left", padx=4)
        ttk.Button(btn_row, text="Deselect",
                   command=self._deselect_node).pack(side="left", padx=4)

        # ── Search ────────────────────────────────────────────────────
        sec2 = Section(inner, "Search")
        sec2.pack(fill="x", padx=4, pady=3)
        s = sec2.inner

        row = ttk.Frame(s)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Attribute", width=RLW, anchor="w").pack(
            side="left", padx=(4, 2))
        ttk.Combobox(
            row, textvariable=self.v_search_attr,
            values=["name", "uid", "pop", "admin_lvl", "admin_dist",
                    "hierarchy", "id"],
            state="readonly", width=12,
        ).pack(side="left")

        row2 = ttk.Frame(s)
        row2.pack(fill="x", pady=2)
        ttk.Label(row2, text="Operator", width=RLW, anchor="w").pack(
            side="left", padx=(4, 2))
        self._search_op_combo = ttk.Combobox(
            row2, textvariable=self.v_search_op,
            values=["contains", "=", ">", "<", ">=", "<="],
            state="readonly", width=10,
        )
        self._search_op_combo.pack(side="left")

        row3 = ttk.Frame(s)
        row3.pack(fill="x", pady=2)
        ttk.Label(row3, text="Value", width=RLW, anchor="w").pack(
            side="left", padx=(4, 2))
        ttk.Entry(row3, textvariable=self.v_search_query, width=16).pack(side="left")

        btn_row = ttk.Frame(s)
        btn_row.pack(fill="x", pady=(4, 2))
        ttk.Button(btn_row, text="Search", command=self._on_search).pack(
            side="left", padx=4)
        ttk.Button(btn_row, text="Clear Search", command=self._on_clear_search).pack(
            side="left", padx=4)
        ttk.Label(s, textvariable=self.v_search_result,
                  foreground="#88cc88", wraplength=320).pack(anchor="w", padx=4, pady=2)

    # ── Filter / Reduced View tab ─────────────────────────────────────────

    def _build_filter_tab(self, parent: ttk.Frame) -> None:
        """Reduced-view filter builder.

        Lets the user define one or more attribute conditions.  When applied,
        nodes that match (or do not match, depending on action) are hidden from
        the preview along with all edges connected to them.
        """
        inner = self._make_scrollable(parent)

        # ── Enable toggle ─────────────────────────────────────────────
        top = ttk.Frame(inner)
        top.pack(fill="x", padx=6, pady=(6, 2))
        ttk.Checkbutton(
            top, text="Enable Filter View",
            variable=self.v_filter_enabled,
        ).pack(side="left")

        # ── Logic & Action ────────────────────────────────────────────
        opts = ttk.Frame(inner)
        opts.pack(fill="x", padx=6, pady=2)

        ttk.Label(opts, text="Logic:", anchor="w").pack(side="left", padx=(0, 2))
        ttk.Combobox(
            opts, textvariable=self.v_filter_logic,
            values=["AND", "OR"], state="readonly", width=5,
        ).pack(side="left", padx=(0, 10))

        ttk.Label(opts, text="Action:", anchor="w").pack(side="left", padx=(0, 2))
        ttk.Radiobutton(
            opts, text="Hide matching",
            variable=self.v_filter_action, value="hide_matching",
        ).pack(side="left")
        ttk.Radiobutton(
            opts, text="Show only matching",
            variable=self.v_filter_action, value="show_only",
        ).pack(side="left", padx=(4, 0))

        # ── Column headers ────────────────────────────────────────────
        hdr = ttk.Frame(inner)
        hdr.pack(fill="x", padx=6, pady=(6, 0))
        ttk.Label(hdr, text="Attribute", width=11, anchor="w").pack(side="left")
        ttk.Label(hdr, text="Op",       width=8,  anchor="w").pack(side="left")
        ttk.Label(hdr, text="Value",    width=9,  anchor="w").pack(side="left")
        ttk.Label(hdr, text="To",       width=8,  anchor="w").pack(side="left")
        ttk.Separator(inner, orient="horizontal").pack(fill="x", padx=6)

        # ── Condition rows container ───────────────────────────────────
        self._filter_container = ttk.Frame(inner)
        self._filter_container.pack(fill="x", padx=6, pady=2)

        # ── Add condition button ───────────────────────────────────────
        ttk.Button(
            inner, text="+ Add Condition",
            command=lambda: self._add_filter_row(self._filter_container),
        ).pack(anchor="w", padx=6, pady=(4, 2))

        # ── Apply / Clear ──────────────────────────────────────────────
        btn_row = ttk.Frame(inner)
        btn_row.pack(fill="x", padx=6, pady=4)
        ttk.Button(btn_row, text="Apply Filter",
                   command=self._apply_filter_and_preview).pack(side="left", padx=(0, 6))
        ttk.Button(btn_row, text="Clear All",
                   command=self._clear_filter).pack(side="left")

        # ── Status ────────────────────────────────────────────────────
        ttk.Label(inner, textvariable=self.v_filter_status,
                  foreground="#88aaff", wraplength=320,
                  justify="left").pack(anchor="w", padx=6, pady=2)

        ttk.Label(inner, text=(
            "Tip: Use AND to narrow down (must meet all conditions).\n"
            "Use OR to broaden (must meet at least one condition).\n"
            "\"between\" range uses the Val and To fields together."
        ), foreground="#888888", wraplength=330,
            justify="left").pack(anchor="w", padx=6, pady=(4, 6))

        # Seed with one blank condition row so the tab isn't empty
        self._add_filter_row(self._filter_container)

    def _add_filter_row(self, container: ttk.Frame) -> None:
        """Append a new condition row to the filter container."""
        _ATTRS = ["pop", "admin_lvl", "admin_dist", "hierarchy",
                  "r", "arm_dist", "name", "uid", "id"]
        _OPS   = ["=", ">", "<", ">=", "<=", "between", "contains"]

        attr_var = tk.StringVar(value="pop")
        op_var   = tk.StringVar(value="=")
        val_var  = tk.StringVar(value="")
        val2_var = tk.StringVar(value="")

        row = ttk.Frame(container)
        row.pack(fill="x", pady=2)

        ttk.Combobox(row, textvariable=attr_var,
                     values=_ATTRS, state="readonly", width=10).pack(side="left", padx=(0, 2))

        op_cb = ttk.Combobox(row, textvariable=op_var,
                             values=_OPS, state="readonly", width=7)
        op_cb.pack(side="left", padx=(0, 2))

        ttk.Entry(row, textvariable=val_var, width=8).pack(side="left", padx=(0, 2))

        to_label = ttk.Label(row, text="to")
        to_label.pack(side="left", padx=(0, 2))

        val2_entry = ttk.Entry(row, textvariable=val2_var, width=8)
        val2_entry.pack(side="left", padx=(0, 4))

        row_info = dict(
            attr_var=attr_var, op_var=op_var,
            val_var=val_var, val2_var=val2_var,
            frame=row, to_label=to_label, val2_entry=val2_entry,
        )

        # Show/hide the "to" field based on operator
        def _op_changed(*_):
            if op_var.get() == "between":
                to_label.pack(side="left", padx=(0, 2))
                val2_entry.pack(side="left", padx=(0, 4))
            else:
                to_label.pack_forget()
                val2_entry.pack_forget()

        op_var.trace_add("write", _op_changed)
        _op_changed()   # apply initial state

        ttk.Button(
            row, text="✕", width=2,
            command=lambda ri=row_info: self._remove_filter_row(ri),
        ).pack(side="left")

        self._filter_rows.append(row_info)

    def _remove_filter_row(self, row_info: dict) -> None:
        """Remove a condition row from the filter list."""
        row_info["frame"].destroy()
        if row_info in self._filter_rows:
            self._filter_rows.remove(row_info)

    def _compute_hidden_nodes(self) -> set:
        """Evaluate the active filter and return a set of node indices to hide.

        Returns an empty set when the filter is disabled or no data is loaded.
        """
        if not self.v_filter_enabled.get() or self._nodes_df is None:
            return set()
        if not self._filter_rows:
            return set()

        n = len(self._nodes_df)
        masks: list[np.ndarray] = []

        for row_info in self._filter_rows:
            attr = row_info["attr_var"].get()
            op   = row_info["op_var"].get()
            val  = row_info["val_var"].get().strip()

            if not attr or not val:
                continue
            if attr not in self._nodes_df.columns:
                continue

            col = self._nodes_df[attr]
            try:
                if op == "contains":
                    mask = col.astype(str).str.contains(val, case=False, na=False)
                elif op == "=":
                    try:
                        mask = col == float(val)
                    except ValueError:
                        mask = col.astype(str) == val
                elif op == ">":
                    mask = col.astype(float) > float(val)
                elif op == "<":
                    mask = col.astype(float) < float(val)
                elif op == ">=":
                    mask = col.astype(float) >= float(val)
                elif op == "<=":
                    mask = col.astype(float) <= float(val)
                elif op == "between":
                    val2 = row_info["val2_var"].get().strip()
                    if not val2:
                        continue
                    lo, hi = float(val), float(val2)
                    mask = (col.astype(float) >= lo) & (col.astype(float) <= hi)
                else:
                    continue
                masks.append(mask.values.astype(bool))
            except Exception:
                continue

        if not masks:
            return set()

        # Combine with AND or OR
        if self.v_filter_logic.get() == "AND":
            combined = np.ones(n, dtype=bool)
            for m in masks:
                combined &= m
        else:  # OR
            combined = np.zeros(n, dtype=bool)
            for m in masks:
                combined |= m

        action = self.v_filter_action.get()
        if action == "show_only":
            hidden_mask = ~combined       # hide non-matching
        else:                             # hide_matching
            hidden_mask = combined        # hide matching

        return set(int(i) for i in np.where(hidden_mask)[0])

    def _apply_filter_and_preview(self) -> None:
        """Evaluate the filter, update status, and refresh the preview."""
        if self._nodes_df is None:
            messagebox.showwarning("No data", "Generate or load data first.")
            return
        hidden = self._compute_hidden_nodes()
        total  = len(self._nodes_df)
        visible = total - len(hidden)
        self.v_filter_status.set(
            f"Filter {'active' if self.v_filter_enabled.get() else 'inactive'}  —  "
            f"{visible:,} of {total:,} nodes visible."
        )
        if not self._is_busy() and self._canvas_widget is not None:
            self._set_busy(True)
            self._worker = threading.Thread(target=self._preview_worker, daemon=True)
            self._worker.start()

    def _clear_filter(self) -> None:
        """Remove all condition rows and disable the filter."""
        for row_info in list(self._filter_rows):
            row_info["frame"].destroy()
        self._filter_rows.clear()
        self.v_filter_enabled.set(False)
        self.v_filter_status.set("")
        if self._filter_container is not None:
            self._add_filter_row(self._filter_container)

    # ── Preview panel (centre) ─────────────────────────────────────────────

    def _build_preview_panel(self, parent: ttk.Frame) -> None:
        self._preview_frame = ttk.Frame(parent)
        self._preview_frame.pack(fill="both", expand=True)

        self._placeholder = ttk.Label(
            self._preview_frame,
            text=(
                "Press  Generate  to create galaxy data,\n"
                "then  Preview  to render it here.\n\n"
                "Click any node to inspect its attributes."
            ),
            anchor="center",
            justify="center",
        )
        self._placeholder.pack(expand=True)

    # ── Action bar ────────────────────────────────────────────────────────

    def _build_action_bar(self) -> None:
        bar = ttk.Frame(self.root)
        bar.pack(side="bottom", fill="x", padx=6, pady=(0, 6))

        self.btn_generate = ttk.Button(bar, text="Generate",
                                       command=self._on_generate, width=12)
        self.btn_generate.pack(side="left", padx=(0, 4))

        self.btn_preview = ttk.Button(bar, text="Preview",
                                      command=self._on_preview, width=12)
        self.btn_preview.pack(side="left", padx=4)

        ttk.Separator(bar, orient="vertical").pack(side="left", fill="y", padx=8, pady=4)

        self.btn_png = ttk.Button(bar, text="Export PNG…",
                                  command=self._on_export_png, width=13)
        self.btn_png.pack(side="left", padx=4)

        self.btn_svg = ttk.Button(bar, text="Export SVG…",
                                  command=self._on_export_svg, width=13)
        self.btn_svg.pack(side="left", padx=4)

        self._status_var = tk.StringVar(value="Ready.")
        ttk.Label(bar, textvariable=self._status_var, anchor="w").pack(
            side="left", padx=12)

        self._progress = ttk.Progressbar(bar, mode="indeterminate", length=110)
        self._progress.pack(side="right", padx=4)

    # ── Helpers ───────────────────────────────────────────────────────────

    def _browse_out_dir(self) -> None:
        path = filedialog.askdirectory(title="Choose output directory")
        if path:
            self.v_out_dir.set(path)

    def _build_config(self) -> GalaxyConfig:
        return GalaxyConfig(
            n_nodes                = self.v_n_nodes.get(),
            r_disk                 = self.v_r_disk.get(),
            r_core                 = self.v_r_core.get(),
            l_max                  = self.v_l_max.get(),
            target_degree          = self.v_target_degree.get(),
            node_connection_chance = self.v_node_connection_chance.get(),
            n_arms                 = self.v_n_arms.get(),
            arm_b                  = self.v_arm_b.get(),
            arm_sigma              = self.v_arm_sigma.get(),
            arm_base               = self.v_arm_base.get(),
            r_scale                = self.v_r_scale.get(),
            boost                  = self.v_boost.get(),
            seed                   = self.v_seed.get(),
            out_dir                = self.v_out_dir.get(),
            write_gexf             = self.v_write_gexf.get(),
            # Population
            pop_core_dispersal     = self.v_pop_core_dispersal.get(),
            pop_dispersal          = self.v_pop_dispersal.get(),
            # Admin exact counts
            admin_count_1          = self.v_admin_count_1.get(),
            admin_count_2          = self.v_admin_count_2.get(),
            admin_count_3          = self.v_admin_count_3.get(),
            admin_count_4          = self.v_admin_count_4.get(),
            admin_count_5          = self.v_admin_count_5.get(),
            # Admin ratios
            admin_ratio_21         = self.v_admin_ratio_21.get(),
            admin_ratio_32         = self.v_admin_ratio_32.get(),
            admin_ratio_43         = self.v_admin_ratio_43.get(),
            admin_ratio_54         = self.v_admin_ratio_54.get(),
            # Admin spatial spacing
            admin_coverage_scale   = self.v_admin_coverage_scale.get(),
            admin_sep_scale        = self.v_admin_sep_scale.get(),
        )

    def _build_plot_args(self,
                         selected: Optional[list[int]] = None,
                         found: Optional[list[int]] = None) -> types.SimpleNamespace:
        ns = types.SimpleNamespace(
            out_dir             = self.v_out_dir.get(),
            color_by            = self.v_color_by.get(),
            node_size           = self.v_node_size.get(),
            node_color          = self.v_node_color.get(),
            gradient_low_color  = self.v_grad_low.get(),
            gradient_high_color = self.v_grad_high.get(),
            no_edges            = self.v_no_edges.get(),
            edge_color          = self.v_edge_color.get(),
            edge_width          = self.v_edge_width.get(),
            edge_alpha          = self.v_edge_alpha.get(),
            r_disk              = self.v_r_disk.get(),
            r_core              = self.v_r_core.get(),
            n_arms              = self.v_n_arms.get(),
            arm_b               = self.v_arm_b.get(),
            r_arm_start         = 3.0,
            selected_nodes      = selected or [],
            found_nodes         = found or [],
            hidden_nodes        = list(self._compute_hidden_nodes()),
            node_mask           = self._compute_node_mask(),
        )
        return ns

    def _compute_node_mask(self) -> Optional[np.ndarray]:
        """Return a boolean visibility array from the active filter conditions.

        Returns None when the filter is inactive or has no valid conditions
        (meaning all nodes are visible).  Returns a bool ndarray of length
        n_nodes otherwise.
        """
        if not self.v_filter_active.get():
            return None
        if self._nodes_df is None:
            return None

        conditions: list[np.ndarray] = []
        for cond in self._filter_conds:
            if not cond["enabled"].get():
                continue
            attr = cond["attr"].get()
            op   = cond["op"].get()
            val  = cond["val"].get().strip()
            if not attr or not val or attr not in self._nodes_df.columns:
                continue
            col = self._nodes_df[attr]
            try:
                if op == "contains":
                    mask = col.astype(str).str.contains(val, case=False, na=False).values
                elif op == "=":
                    try:
                        mask = (col == float(val)).values
                    except ValueError:
                        mask = (col.astype(str) == val).values
                elif op == ">":
                    mask = (col.astype(float) > float(val)).values
                elif op == "<":
                    mask = (col.astype(float) < float(val)).values
                elif op == ">=":
                    mask = (col.astype(float) >= float(val)).values
                elif op == "<=":
                    mask = (col.astype(float) <= float(val)).values
                else:
                    continue
                conditions.append(mask)
            except Exception:
                continue

        if not conditions:
            return None

        combined = conditions[0].copy()
        logic = self.v_filter_logic.get()
        for cond_mask in conditions[1:]:
            if logic == "AND":
                combined = combined & cond_mask
            else:
                combined = combined | cond_mask

        if self.v_filter_mode.get() == "hide_matching":
            combined = ~combined

        return combined

    def _on_apply_filter(self) -> None:
        """Re-render the preview with the current filter applied."""
        if not self._is_busy() and self._canvas_widget is not None:
            self._set_busy(True)
            self._worker = threading.Thread(target=self._preview_worker, daemon=True)
            self._worker.start()

    def _set_busy(self, busy: bool) -> None:
        state = "disabled" if busy else "normal"
        for btn in (self.btn_generate, self.btn_preview,
                    self.btn_png, self.btn_svg):
            btn.configure(state=state)
        if busy:
            self._progress.start(10)
        else:
            self._progress.stop()

    def _status(self, msg: str) -> None:
        self._status_var.set(msg)

    def _is_busy(self) -> bool:
        return self._worker is not None and self._worker.is_alive()

    def _load_data(self) -> bool:
        """Load nodes.csv + edges.csv into memory. Returns True on success."""
        out_dir = self.v_out_dir.get()
        nodes_path = os.path.join(out_dir, "nodes.csv")
        edges_path = os.path.join(out_dir, "edges.csv")
        if not os.path.exists(nodes_path):
            return False
        try:
            self._nodes_df = pd.read_csv(nodes_path)
            # Ensure string types for uid/name even if CSV has them as numeric
            if "uid" in self._nodes_df.columns:
                self._nodes_df["uid"] = self._nodes_df["uid"].astype(str)
            if "name" in self._nodes_df.columns:
                self._nodes_df["name"] = self._nodes_df["name"].fillna("").astype(str)
            self._edges_df = (pd.read_csv(edges_path)
                              if os.path.exists(edges_path) else pd.DataFrame())
            xy = self._nodes_df[["x", "y"]].values
            self._node_tree = cKDTree(xy)
            return True
        except Exception:
            return False

    def _save_nodes_csv(self) -> None:
        """Persist the in-memory nodes_df back to nodes.csv."""
        if self._nodes_df is None:
            return
        out_dir = self.v_out_dir.get()
        nodes_path = os.path.join(out_dir, "nodes.csv")
        self._nodes_df.to_csv(nodes_path, index=False)

    # ── Generate action ───────────────────────────────────────────────────

    def _on_generate(self) -> None:
        if self._is_busy():
            return
        self._set_busy(True)
        self._status("Generating galaxy…")
        self._worker = threading.Thread(target=self._generate_worker, daemon=True)
        self._worker.start()

    def _generate_worker(self) -> None:
        try:
            cfg = self._build_config()
            os.makedirs(cfg.out_dir, exist_ok=True)
            gen = GalaxyGenerator(cfg)
            nodes, edges = gen.run()

            # Persist params.json
            params_path = os.path.join(cfg.out_dir, "params.json")
            params = {f: getattr(cfg, f) for f in cfg.__dataclass_fields__}
            with open(params_path, "w") as fp:
                json.dump(params, fp, indent=2)

            # Load into memory
            self._nodes_df = nodes
            self._edges_df = edges
            xy = nodes[["x", "y"]].values
            self._node_tree = cKDTree(xy)
            self._selected_node = None
            self._found_nodes = []

            n, d = cfg.n_nodes, cfg.out_dir
            self.root.after(0, lambda: self._status(
                f"Done — {n:,} nodes written to '{d}'.  Press Preview to render."))
        except Exception as exc:
            msg = str(exc)
            self.root.after(0, lambda: (
                self._status(f"Generation failed: {msg}"),
                messagebox.showerror("Generation failed", msg),
            ))
        finally:
            self.root.after(0, lambda: self._set_busy(False))

    # ── Apply attributes action ───────────────────────────────────────────

    def _on_apply_attrs(self) -> None:
        """Re-run Stage C on existing spatial data without regenerating layout."""
        out_dir = self.v_out_dir.get()
        if not os.path.exists(os.path.join(out_dir, "nodes.csv")):
            messagebox.showwarning("No data",
                                   "Run Generate first to produce nodes.csv.")
            return
        if self._is_busy():
            return

        # Ensure data is loaded
        if self._nodes_df is None:
            self._load_data()

        self._set_busy(True)
        self._status("Applying attributes…")

        def _worker():
            try:
                cfg = self._build_config()
                gen = GalaxyGenerator(cfg)
                nodes = self._nodes_df.copy()
                edges = self._edges_df if self._edges_df is not None else pd.DataFrame()

                # Drop old attribute columns so assign_attributes regenerates them
                for col in ["uid", "name", "pop", "admin_lvl", "admin_dist",
                            "hierarchy"]:
                    if col in nodes.columns:
                        nodes = nodes.drop(columns=[col])

                nodes = gen.assign_attributes(nodes, edges)
                self._nodes_df = nodes
                self._save_nodes_csv()

                # Update inspector if a node is selected
                sel = self._selected_node
                def _refresh():
                    if sel is not None:
                        self._populate_inspector(sel)
                    self._status("Attributes applied and saved.")
                self.root.after(0, _refresh)
            except Exception as exc:
                msg = str(exc)
                self.root.after(0, lambda: (
                    self._status(f"Apply failed: {msg}"),
                    messagebox.showerror("Apply attributes failed", msg),
                ))
            finally:
                self.root.after(0, lambda: self._set_busy(False))

        threading.Thread(target=_worker, daemon=True).start()

    # ── Preview action ────────────────────────────────────────────────────

    def _on_preview(self) -> None:
        out_dir = self.v_out_dir.get()
        if not os.path.exists(os.path.join(out_dir, "nodes.csv")):
            messagebox.showwarning("No data",
                                   f"nodes.csv not found in '{out_dir}'.\n"
                                   "Run Generate first.")
            return
        if self._is_busy():
            return
        # Load data for node interaction if not already loaded
        if self._nodes_df is None:
            self._load_data()
        self._set_busy(True)
        self._status("Rendering preview…")
        self._worker = threading.Thread(target=self._preview_worker, daemon=True)
        self._worker.start()

    def _preview_worker(self) -> None:
        try:
            args = self._build_plot_args(
                selected=[self._selected_node] if self._selected_node is not None else [],
                found=self._found_nodes,
            )
            fig = draw_galaxy(args)
            self.root.after(0, lambda: self._display_figure(fig))
            self.root.after(0, lambda: self._status("Preview ready.  Click a node to inspect."))
        except Exception as exc:
            msg = str(exc)
            self.root.after(0, lambda: (
                self._status(f"Preview failed: {msg}"),
                messagebox.showerror("Preview failed", msg),
            ))
        finally:
            self.root.after(0, lambda: self._set_busy(False))

    def _display_figure(self, fig: plt.Figure) -> None:
        """Swap the embedded matplotlib canvas to show a new figure."""
        self._placeholder.pack_forget()

        if self._canvas_widget is not None:
            self._canvas_widget.get_tk_widget().destroy()
        if self._toolbar_frame is not None:
            self._toolbar_frame.destroy()
        if self._current_fig is not None:
            plt.close(self._current_fig)

        self._current_fig = fig

        canvas = FigureCanvasTkAgg(fig, master=self._preview_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

        toolbar_frame = ttk.Frame(self._preview_frame)
        toolbar_frame.pack(fill="x")
        NavigationToolbar2Tk(canvas, toolbar_frame).update()

        self._attach_scroll_zoom(canvas)
        self._attach_pan(canvas)
        self._attach_node_click(canvas)

        self._canvas_widget = canvas
        self._toolbar_frame = toolbar_frame

    # ── Scroll-wheel zoom ─────────────────────────────────────────────────

    def _attach_scroll_zoom(self, canvas: FigureCanvasTkAgg) -> None:
        """Scroll-wheel zoom centred on the cursor position.

        Both axes are scaled by the same factor so the galaxy stays
        undistorted regardless of window shape.  The axes tick labels
        update to reflect the new coordinate range.
        """
        FACTOR = 1.2

        def _do_zoom(factor: float, ex: int, ey: int) -> None:
            if not canvas.figure.axes:
                return
            ax = canvas.figure.axes[0]
            # Convert tkinter pixel coords (top-left origin) to
            # matplotlib display coords (bottom-left origin).
            fig_h = canvas.figure.bbox.height
            xd = float(ex)
            yd = float(fig_h) - float(ey)
            if ax.bbox.contains(xd, yd):
                inv = ax.transData.inverted()
                cx, cy = inv.transform((xd, yd))
            else:
                cx = sum(ax.get_xlim()) / 2.0
                cy = sum(ax.get_ylim()) / 2.0
            xl, yl = ax.get_xlim(), ax.get_ylim()
            # Scale both axes by identical factor (preserves aspect ratio).
            x_half = (xl[1] - xl[0]) / 2.0 / factor
            y_half = (yl[1] - yl[0]) / 2.0 / factor
            ax.set_xlim([cx - x_half, cx + x_half])
            ax.set_ylim([cy - y_half, cy + y_half])
            canvas.draw()

        w = canvas.get_tk_widget()
        w.bind("<MouseWheel>",
               lambda e: _do_zoom(FACTOR if e.delta > 0 else 1 / FACTOR, e.x, e.y))
        w.bind("<Button-4>", lambda e: _do_zoom(FACTOR,     e.x, e.y))
        w.bind("<Button-5>", lambda e: _do_zoom(1 / FACTOR, e.x, e.y))

    # ── Drag-to-pan ───────────────────────────────────────────────────────

    def _attach_pan(self, canvas: FigureCanvasTkAgg) -> None:
        """Right-click-drag (or middle-click-drag) pans the galaxy view.

        Tracks the cumulative pixel delta since the last motion event and
        converts it to data-space displacement using the live axis scale.
        Left-click is reserved for node selection; the NavigationToolbar's
        own Pan/Zoom tools also continue to work normally.
        """
        state: dict = {"active": False, "last": None}

        def _start(e: tk.Event) -> None:
            state["active"] = True
            state["last"]   = (e.x, e.y)

        def _move(e: tk.Event) -> None:
            if not state["active"] or state["last"] is None:
                return
            if not canvas.figure.axes:
                return
            ax   = canvas.figure.axes[0]
            bbox = ax.get_window_extent()
            if bbox.width == 0 or bbox.height == 0:
                return
            xl, yl    = ax.get_xlim(), ax.get_ylim()
            dx_pix    = e.x - state["last"][0]
            dy_pix    = e.y - state["last"][1]
            x_scale   = (xl[1] - xl[0]) / bbox.width
            y_scale   = (yl[1] - yl[0]) / bbox.height
            # Screen y is inverted relative to data y.
            shift_x   = -dx_pix * x_scale
            shift_y   =  dy_pix * y_scale
            ax.set_xlim([xl[0] + shift_x, xl[1] + shift_x])
            ax.set_ylim([yl[0] + shift_y, yl[1] + shift_y])
            state["last"] = (e.x, e.y)
            canvas.draw()

        def _stop(e: tk.Event) -> None:
            state["active"] = False
            state["last"]   = None

        w = canvas.get_tk_widget()
        # Middle-click drag
        w.bind("<Button-2>",        _start)
        w.bind("<B2-Motion>",       _move)
        w.bind("<ButtonRelease-2>", _stop)
        # Right-click drag (more accessible on systems without middle button)
        w.bind("<Button-3>",        _start)
        w.bind("<B3-Motion>",       _move)
        w.bind("<ButtonRelease-3>", _stop)

    # ── Node click / selection ────────────────────────────────────────────

    def _attach_node_click(self, canvas: FigureCanvasTkAgg) -> None:
        """Left-click selects the nearest node within a proximity threshold."""
        def _on_click(event):
            # Only handle left-clicks; right/middle are reserved for pan.
            if event.button != 1:
                return
            if event.inaxes is None:
                return
            if self._node_tree is None or self._nodes_df is None:
                return
            ax = event.inaxes
            click_xy = np.array([[event.xdata, event.ydata]])
            dist, idx = self._node_tree.query(click_xy, k=1)
            idx  = int(idx[0])
            dist = float(dist[0])
            # Threshold: ~1.5% of current axis width
            xlim      = ax.get_xlim()
            threshold = (xlim[1] - xlim[0]) * 0.015
            if dist <= threshold:
                self._selected_node = idx
                self.root.after(0, lambda: self._on_node_selected(idx))

        canvas.mpl_connect("button_press_event", _on_click)

    def _on_node_selected(self, idx: int) -> None:
        """Called on the main thread when the user clicks a node."""
        self._populate_inspector(idx)
        # Refresh preview with selection highlight (non-blocking quick redraw)
        if self._is_busy():
            return
        self._set_busy(True)
        self._status(f"Node {idx} selected.")
        self._worker = threading.Thread(target=self._preview_worker, daemon=True)
        self._worker.start()

    def _populate_inspector(self, idx: int) -> None:
        """Fill the inspector panel with data from node *idx*."""
        if self._nodes_df is None or idx >= len(self._nodes_df):
            return
        row = self._nodes_df.iloc[idx]
        self.v_insp_uid.set(str(row.get("uid", "")) if "uid" in self._nodes_df.columns else "")
        self.v_insp_name.set(str(row.get("name", "")) if "name" in self._nodes_df.columns else "")
        self.v_insp_pop.set(int(row.get("pop", 0)) if "pop" in self._nodes_df.columns else 0)
        self.v_insp_admin_lvl.set(int(row.get("admin_lvl", 0)) if "admin_lvl" in self._nodes_df.columns else 0)
        try:
            adist = int(row.get("admin_dist", -1)) if "admin_dist" in self._nodes_df.columns else -1
        except (ValueError, TypeError):
            adist = -1
        self.v_insp_admin_dist.set(str(adist) if adist >= 0 else "N/A")
        try:
            hier = int(row.get("hierarchy", -1)) if "hierarchy" in self._nodes_df.columns else -1
        except (ValueError, TypeError):
            hier = -1
        self.v_insp_hierarchy.set(str(hier) if hier >= 0 else "N/A")

    def _on_save_node(self) -> None:
        """Save inspector edits back to the in-memory DataFrame and CSV."""
        idx = self._selected_node
        if idx is None or self._nodes_df is None:
            messagebox.showinfo("No node selected", "Click a node first.")
            return
        if "name" in self._nodes_df.columns:
            self._nodes_df.at[idx, "name"] = self.v_insp_name.get()
        if "pop" in self._nodes_df.columns:
            self._nodes_df.at[idx, "pop"] = max(0, min(100, self.v_insp_pop.get()))
        if "admin_lvl" in self._nodes_df.columns:
            self._nodes_df.at[idx, "admin_lvl"] = max(0, min(5, self.v_insp_admin_lvl.get()))
        self._save_nodes_csv()
        self._status(f"Node {idx} saved.")

    def _on_recalc_admin_dist(self) -> None:
        """Recalculate admin_dist for all nodes after manual edits."""
        if self._nodes_df is None:
            messagebox.showinfo("No data", "Generate or load data first.")
            return
        if "admin_lvl" not in self._nodes_df.columns:
            messagebox.showinfo("No admin data", "Apply attributes first.")
            return
        edges = self._edges_df if self._edges_df is not None else pd.DataFrame()
        new_dist = compute_admin_dist(self._nodes_df, edges)
        self._nodes_df["admin_dist"] = new_dist
        self._save_nodes_csv()
        if self._selected_node is not None:
            self._populate_inspector(self._selected_node)
        self._status("Admin distances recalculated and saved.")

    def _deselect_node(self) -> None:
        self._selected_node = None
        self.v_insp_uid.set("")
        self.v_insp_name.set("")
        self.v_insp_pop.set(0)
        self.v_insp_admin_lvl.set(0)
        self.v_insp_admin_dist.set("")
        self.v_insp_hierarchy.set("")
        self._status("Deselected.")
        if not self._is_busy() and self._canvas_widget is not None:
            self._set_busy(True)
            self._worker = threading.Thread(target=self._preview_worker, daemon=True)
            self._worker.start()

    # ── Search ────────────────────────────────────────────────────────────

    def _on_search(self) -> None:
        if self._nodes_df is None:
            messagebox.showwarning("No data", "Generate or load data first.")
            return

        attr  = self.v_search_attr.get()
        op    = self.v_search_op.get()
        query = self.v_search_query.get().strip()

        if not query:
            self.v_search_result.set("Enter a search value.")
            return
        if attr not in self._nodes_df.columns:
            self.v_search_result.set(f"Column '{attr}' not in data.")
            return

        col = self._nodes_df[attr]

        try:
            if op == "contains":
                mask = col.astype(str).str.contains(query, case=False, na=False)
            elif op == "=":
                try:
                    mask = col == float(query)
                except ValueError:
                    mask = col.astype(str) == query
            elif op == ">":
                mask = col.astype(float) > float(query)
            elif op == "<":
                mask = col.astype(float) < float(query)
            elif op == ">=":
                mask = col.astype(float) >= float(query)
            elif op == "<=":
                mask = col.astype(float) <= float(query)
            else:
                mask = pd.Series([False] * len(col))
        except Exception as exc:
            self.v_search_result.set(f"Search error: {exc}")
            return

        hits = list(self._nodes_df.index[mask])
        self._found_nodes = hits
        self.v_search_result.set(f"{len(hits)} node(s) found.")

        # Auto-select first result in inspector
        if hits:
            self._selected_node = hits[0]
            self._populate_inspector(hits[0])

        # Refresh preview with highlights
        if not self._is_busy() and self._canvas_widget is not None:
            self._set_busy(True)
            self._worker = threading.Thread(target=self._preview_worker, daemon=True)
            self._worker.start()

    def _on_clear_search(self) -> None:
        self._found_nodes = []
        self.v_search_result.set("")
        self.v_search_query.set("")
        if not self._is_busy() and self._canvas_widget is not None:
            self._set_busy(True)
            self._worker = threading.Thread(target=self._preview_worker, daemon=True)
            self._worker.start()

    # ── Export actions ────────────────────────────────────────────────────

    def _on_export_png(self) -> None:
        self._export("png")

    def _on_export_svg(self) -> None:
        self._export("svg")

    def _export(self, fmt: str) -> None:
        out_dir = self.v_out_dir.get()
        if not os.path.exists(os.path.join(out_dir, "nodes.csv")):
            messagebox.showwarning("No data", "Run Generate first.")
            return
        if self._is_busy():
            return

        filetypes = [(f"{fmt.upper()} files", f"*.{fmt}"), ("All files", "*.*")]
        path = filedialog.asksaveasfilename(
            defaultextension=f".{fmt}",
            filetypes=filetypes,
            initialfile=f"galaxy.{fmt}",
            title=f"Export as {fmt.upper()}",
        )
        if not path:
            return

        self._set_busy(True)
        self._status(f"Exporting {fmt.upper()}…")

        def _worker() -> None:
            try:
                args = self._build_plot_args()
                fig  = draw_galaxy(args)
                save_kw: dict = dict(bbox_inches="tight",
                                     facecolor=fig.get_facecolor())
                if fmt == "png":
                    save_kw["dpi"] = 150
                fig.savefig(path, format=fmt, **save_kw)
                plt.close(fig)
                self.root.after(0, lambda: self._status(f"Saved → {path}"))
            except Exception as exc:
                msg = str(exc)
                self.root.after(0, lambda: (
                    self._status(f"Export failed: {msg}"),
                    messagebox.showerror("Export failed", msg),
                ))
            finally:
                self.root.after(0, lambda: self._set_busy(False))

        threading.Thread(target=_worker, daemon=True).start()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    root = tk.Tk()
    GalaxyGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
