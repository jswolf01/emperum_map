"""
galaxy_gui.py
=============
Tkinter GUI front-end for the Emperum galaxy generator.

Wraps all generation and visualisation functionality in a desktop interface
with collapsible parameter sections, sensible input widgets (sliders,
spinboxes, HSV colour pickers), a live embedded preview, and PNG/SVG export.

Existing CLI tools (run_generate.py, plot_debug.py) are unchanged and remain
fully functional alongside this file.

Usage
-----
    python galaxy_gui.py

Dependencies
------------
Same as the core generator (numpy, pandas, scipy, matplotlib) plus tkinter,
which is bundled with the standard Python installer.  On Ubuntu/Debian it can
be installed with:  sudo apt-get install python3-tk
"""

from __future__ import annotations

import json
import os
import threading
import types
from typing import Optional

import tkinter as tk
from tkinter import ttk, colorchooser, filedialog, messagebox

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from galaxygen import GalaxyConfig, GalaxyGenerator
from plot_debug import draw_galaxy


# ---------------------------------------------------------------------------
# Reusable compound widgets
# ---------------------------------------------------------------------------

class SliderEntry(ttk.Frame):
    """Linked horizontal scale + spinbox for a numeric parameter.

    The slider and spinbox share a single tkinter Variable and stay in sync.
    Values are snapped to the nearest ``step`` increment on slider drag.
    """

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
        self._busy = False  # reentrance guard

        ttk.Label(self, text=label, width=label_width, anchor="w").grid(
            row=0, column=0, sticky="w", padx=(4, 2), pady=1,
        )

        self._scale = ttk.Scale(
            self, orient="horizontal", length=150,
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

    # ------------------------------------------------------------------
    def _on_scale(self, _val: str) -> None:
        """Snap the slider value to the nearest step."""
        if self._busy:
            return
        try:
            raw = float(_val)
        except ValueError:
            return
        snapped = round(raw / self._step) * self._step
        # Guard against floating-point drift producing tiny differences
        snapped = round(snapped, 10)
        if abs(raw - snapped) > 1e-9:
            self._busy = True
            self._var.set(snapped)
            self._busy = False

    def _clamp(self, _evt=None) -> None:
        """Clamp a manually typed spinbox value into [lo, hi]."""
        try:
            val = float(self._spin.get())
        except ValueError:
            val = self._lo
        val = max(self._lo, min(self._hi, val))
        self._var.set(round(val / self._step) * self._step)


# ---------------------------------------------------------------------------

class ColorEntry(ttk.Frame):
    """Colour swatch + hex text entry + HSV colour-picker button.

    The swatch reflects the current hex value in real time.
    The 'Pick…' button opens the system colour-picker dialog (which includes
    an HSV/wheel selector on most platforms).
    """

    def __init__(self, parent, label: str, var: tk.StringVar,
                 label_width: int = 24, **kw):
        super().__init__(parent, **kw)
        self._var = var

        ttk.Label(self, text=label, width=label_width, anchor="w").grid(
            row=0, column=0, sticky="w", padx=(4, 2), pady=1,
        )

        # Coloured swatch (clickable shortcut to picker)
        self._swatch = tk.Label(
            self, width=3, relief="sunken", cursor="hand2",
        )
        self._swatch.grid(row=0, column=1, padx=(2, 2))
        self._swatch.bind("<Button-1>", self._open_picker)

        self._entry = ttk.Entry(self, textvariable=var, width=10)
        self._entry.grid(row=0, column=2, padx=2)
        self._entry.bind("<Return>",   self._refresh_swatch)
        self._entry.bind("<FocusOut>", self._refresh_swatch)

        ttk.Button(
            self, text="Pick…", width=6, command=self._open_picker,
        ).grid(row=0, column=3, padx=(2, 4))

        var.trace_add("write", self._refresh_swatch)
        self._refresh_swatch()

    # ------------------------------------------------------------------
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
                color=current, title="Choose colour",
                parent=self,
            )
        except Exception:
            return
        if hexval:
            self._var.set(hexval.lower())


# ---------------------------------------------------------------------------

class Section(ttk.Frame):
    """Collapsible parameter section with a toggle-button header."""

    def __init__(self, parent, title: str, **kw):
        super().__init__(parent, **kw)
        self._open  = True
        self._title = title

        self._btn = ttk.Button(
            self, text=f"▼  {title}", command=self._toggle,
        )
        self._btn.pack(fill="x", padx=2, pady=(4, 0))

        ttk.Separator(self, orient="horizontal").pack(fill="x", padx=2)

        self._inner = ttk.Frame(self, padding=(2, 2, 2, 6))
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
        root.minsize(1020, 700)

        self._current_fig: Optional[plt.Figure]  = None
        self._canvas_widget: Optional[FigureCanvasTkAgg] = None
        self._toolbar_frame: Optional[ttk.Frame] = None
        self._worker: Optional[threading.Thread] = None

        self._build_vars()
        self._build_ui()

    # ── Variable definitions ──────────────────────────────────────────────

    def _build_vars(self) -> None:
        """Create all tkinter control variables initialised to GalaxyConfig defaults."""
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

    # ── UI construction ───────────────────────────────────────────────────

    def _build_ui(self) -> None:
        # Action bar is packed first so it always reserves space at the bottom.
        # The paned window is packed second with expand=True to fill what remains.
        self._build_action_bar()

        # Top-level paned layout: left (params) | right (preview)
        paned = ttk.PanedWindow(self.root, orient="horizontal")
        paned.pack(fill="both", expand=True, padx=6, pady=(6, 0))

        left_outer = ttk.Frame(paned, width=390)
        left_outer.pack_propagate(False)
        paned.add(left_outer, weight=0)

        right_frame = ttk.Frame(paned)
        paned.add(right_frame, weight=1)

        self._build_param_panel(left_outer)
        self._build_preview_panel(right_frame)

    # ── Parameter panel ───────────────────────────────────────────────────

    def _build_param_panel(self, parent: ttk.Frame) -> None:
        """Scrollable left panel containing all collapsible parameter sections."""
        scroll_canvas = tk.Canvas(parent, highlightthickness=0, borderwidth=0)
        vscroll = ttk.Scrollbar(parent, orient="vertical",
                                command=scroll_canvas.yview)
        scroll_canvas.configure(yscrollcommand=vscroll.set)

        vscroll.pack(side="right", fill="y")
        scroll_canvas.pack(side="left", fill="both", expand=True)

        inner = ttk.Frame(scroll_canvas)
        win_id = scroll_canvas.create_window((0, 0), window=inner, anchor="nw")

        inner.bind(
            "<Configure>",
            lambda _e: scroll_canvas.configure(
                scrollregion=scroll_canvas.bbox("all")
            ),
        )
        scroll_canvas.bind(
            "<Configure>",
            lambda e: scroll_canvas.itemconfigure(win_id, width=e.width),
        )

        # Mouse-wheel scrolling (cross-platform)
        def _wheel(evt):
            if evt.delta:
                scroll_canvas.yview_scroll(int(-1 * evt.delta / 120), "units")
        scroll_canvas.bind_all("<MouseWheel>", _wheel)
        scroll_canvas.bind_all(
            "<Button-4>", lambda _e: scroll_canvas.yview_scroll(-1, "units"))
        scroll_canvas.bind_all(
            "<Button-5>", lambda _e: scroll_canvas.yview_scroll(1, "units"))

        LW = 26  # consistent label width across all rows

        # ── Spatial ───────────────────────────────────────────────────
        sec = Section(inner, "Spatial")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        SliderEntry(s, "Star systems (n_nodes)",
                    self.v_n_nodes, 100, 20_000, 100,
                    label_width=LW).pack(fill="x")
        SliderEntry(s, "Disk radius (r_disk)",
                    self.v_r_disk, 10.0, 500.0, 1.0,
                    label_width=LW).pack(fill="x")
        SliderEntry(s, "Core radius (r_core)",
                    self.v_r_core, 0.0, 100.0, 0.5,
                    label_width=LW).pack(fill="x")

        # ── Edges ─────────────────────────────────────────────────────
        sec = Section(inner, "Edges (Hyperspace Lanes)")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        SliderEntry(s, "Max lane length (l_max)",
                    self.v_l_max, 1.0, 50.0, 0.5,
                    label_width=LW).pack(fill="x")
        SliderEntry(s, "Avg degree (target_degree)",
                    self.v_target_degree, 1.0, 20.0, 0.5,
                    label_width=LW).pack(fill="x")
        SliderEntry(s, "Connection chance",
                    self.v_node_connection_chance, 0.0, 1.0, 0.05,
                    label_width=LW).pack(fill="x")

        # ── Spiral Arms ───────────────────────────────────────────────
        sec = Section(inner, "Spiral Arms")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        SliderEntry(s, "Number of arms (n_arms)",
                    self.v_n_arms, 1, 8, 1,
                    label_width=LW).pack(fill="x")
        SliderEntry(s, "Arm tightness (arm_b)",
                    self.v_arm_b, 0.01, 2.0, 0.01,
                    label_width=LW).pack(fill="x")
        SliderEntry(s, "Arm width σ (arm_sigma)",
                    self.v_arm_sigma, 0.5, 30.0, 0.5,
                    label_width=LW).pack(fill="x")
        SliderEntry(s, "Inter-arm density (arm_base)",
                    self.v_arm_base, 0.0, 1.0, 0.01,
                    label_width=LW).pack(fill="x")

        # ── Radial Density ────────────────────────────────────────────
        sec = Section(inner, "Radial Density")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        SliderEntry(s, "Scale length (r_scale)",
                    self.v_r_scale, 5.0, 200.0, 1.0,
                    label_width=LW).pack(fill="x")

        # ── Sampling & Seed ───────────────────────────────────────────
        sec = Section(inner, "Sampling & Reproducibility")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner
        SliderEntry(s, "Boost multiplier (boost)",
                    self.v_boost, 1.0, 20.0, 0.5,
                    label_width=LW).pack(fill="x")

        row = ttk.Frame(s)
        row.pack(fill="x", pady=1)
        ttk.Label(row, text="Random seed", width=LW, anchor="w").pack(
            side="left", padx=(4, 2))
        ttk.Spinbox(row, from_=0, to=99_999, increment=1,
                    textvariable=self.v_seed, width=8).pack(side="left")

        # ── Output ────────────────────────────────────────────────────
        sec = Section(inner, "Output")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner

        row = ttk.Frame(s)
        row.pack(fill="x", pady=1)
        ttk.Label(row, text="Output directory", width=LW, anchor="w").pack(
            side="left", padx=(4, 2))
        ttk.Entry(row, textvariable=self.v_out_dir, width=13).pack(side="left")
        ttk.Button(row, text="…", width=3,
                   command=self._browse_out_dir).pack(side="left", padx=2)

        row2 = ttk.Frame(s)
        row2.pack(fill="x", pady=1)
        ttk.Checkbutton(
            row2, text="Export GEXF (requires networkx)",
            variable=self.v_write_gexf,
        ).pack(side="left", padx=(4, 2))

        # ── Node Appearance ───────────────────────────────────────────
        sec = Section(inner, "Node Appearance")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner

        row = ttk.Frame(s)
        row.pack(fill="x", pady=1)
        ttk.Label(row, text="Colour by", width=LW, anchor="w").pack(
            side="left", padx=(4, 2))
        ttk.Combobox(
            row, textvariable=self.v_color_by,
            values=["arm_dist", "r", "none"],
            state="readonly", width=10,
        ).pack(side="left")

        SliderEntry(s, "Node size",
                    self.v_node_size, 0.1, 20.0, 0.1,
                    label_width=LW).pack(fill="x")
        ColorEntry(s, "Flat colour  (mode: none)",
                   self.v_node_color, label_width=LW).pack(fill="x")
        ColorEntry(s, "Gradient low  (near arm)",
                   self.v_grad_low, label_width=LW).pack(fill="x")
        ColorEntry(s, "Gradient high  (inter-arm)",
                   self.v_grad_high, label_width=LW).pack(fill="x")

        # ── Edge Appearance ───────────────────────────────────────────
        sec = Section(inner, "Edge Appearance")
        sec.pack(fill="x", padx=4, pady=3)
        s = sec.inner

        ttk.Checkbutton(
            s, text="Hide edges  (much faster for N > 1 000)",
            variable=self.v_no_edges,
        ).pack(anchor="w", padx=4, pady=2)
        ColorEntry(s, "Edge colour",
                   self.v_edge_color, label_width=LW).pack(fill="x")
        SliderEntry(s, "Edge width",
                    self.v_edge_width, 0.1, 5.0, 0.1,
                    label_width=LW).pack(fill="x")
        SliderEntry(s, "Edge opacity",
                    self.v_edge_alpha, 0.0, 1.0, 0.05,
                    label_width=LW).pack(fill="x")

    # ── Preview panel ─────────────────────────────────────────────────────

    def _build_preview_panel(self, parent: ttk.Frame) -> None:
        self._preview_frame = ttk.Frame(parent)
        self._preview_frame.pack(fill="both", expand=True)

        self._placeholder = ttk.Label(
            self._preview_frame,
            text=(
                "Press  Generate  to create galaxy data,\n"
                "then  Preview  to render it here."
            ),
            anchor="center",
            justify="center",
        )
        self._placeholder.pack(expand=True)

    # ── Action bar ────────────────────────────────────────────────────────

    def _build_action_bar(self) -> None:
        bar = ttk.Frame(self.root)
        bar.pack(side="bottom", fill="x", padx=6, pady=(0, 6))

        self.btn_generate = ttk.Button(
            bar, text="Generate", command=self._on_generate, width=12)
        self.btn_generate.pack(side="left", padx=(0, 4))

        self.btn_preview = ttk.Button(
            bar, text="Preview", command=self._on_preview, width=12)
        self.btn_preview.pack(side="left", padx=4)

        ttk.Separator(bar, orient="vertical").pack(
            side="left", fill="y", padx=8, pady=4)

        self.btn_png = ttk.Button(
            bar, text="Export PNG…", command=self._on_export_png, width=13)
        self.btn_png.pack(side="left", padx=4)

        self.btn_svg = ttk.Button(
            bar, text="Export SVG…", command=self._on_export_svg, width=13)
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
        """Build a GalaxyConfig from the current GUI variable values."""
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
        )

    def _build_plot_args(self) -> types.SimpleNamespace:
        """Build a plot-args namespace from the current GUI variable values.

        The spatial overlay fields (r_disk, r_core, n_arms, arm_b, r_arm_start)
        are taken directly from the generation variables so they always match
        the data — no manual synchronisation required.
        """
        return types.SimpleNamespace(
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
            # Spatial overlays — always mirrored from generation params
            r_disk              = self.v_r_disk.get(),
            r_core              = self.v_r_core.get(),
            n_arms              = self.v_n_arms.get(),
            arm_b               = self.v_arm_b.get(),
            r_arm_start         = 3.0,
        )

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

    # ── Generate action ───────────────────────────────────────────────────

    def _on_generate(self) -> None:
        if self._is_busy():
            return
        self._set_busy(True)
        self._status("Generating galaxy…")
        self._worker = threading.Thread(
            target=self._generate_worker, daemon=True)
        self._worker.start()

    def _generate_worker(self) -> None:
        try:
            cfg = self._build_config()
            os.makedirs(cfg.out_dir, exist_ok=True)
            gen = GalaxyGenerator(cfg)
            gen.run()

            # Persist params.json (mirrors what run_generate.py does)
            params_path = os.path.join(cfg.out_dir, "params.json")
            params = {f: getattr(cfg, f) for f in cfg.__dataclass_fields__}
            with open(params_path, "w") as fp:
                json.dump(params, fp, indent=2)

            n = cfg.n_nodes
            d = cfg.out_dir
            self.root.after(0, lambda: self._status(
                f"Done — {n:,} nodes written to '{d}'.  "
                f"Press Preview to render."
            ))
        except Exception as exc:
            msg = str(exc)
            self.root.after(0, lambda: (
                self._status(f"Generation failed: {msg}"),
                messagebox.showerror("Generation failed", msg),
            ))
        finally:
            self.root.after(0, lambda: self._set_busy(False))

    # ── Preview action ────────────────────────────────────────────────────

    def _on_preview(self) -> None:
        out_dir = self.v_out_dir.get()
        nodes_csv = os.path.join(out_dir, "nodes.csv")
        if not os.path.exists(nodes_csv):
            messagebox.showwarning(
                "No data",
                f"nodes.csv not found in '{out_dir}'.\n"
                "Run Generate first.",
            )
            return
        if self._is_busy():
            return
        self._set_busy(True)
        self._status("Rendering preview…")
        self._worker = threading.Thread(
            target=self._preview_worker, daemon=True)
        self._worker.start()

    def _preview_worker(self) -> None:
        try:
            args = self._build_plot_args()
            fig  = draw_galaxy(args)
            self.root.after(0, lambda: self._display_figure(fig))
            self.root.after(0, lambda: self._status("Preview ready."))
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
        # Remove placeholder on first preview
        self._placeholder.pack_forget()

        # Destroy previous canvas and toolbar
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

        self._canvas_widget  = canvas
        self._toolbar_frame  = toolbar_frame

    # ── Scroll-wheel zoom ─────────────────────────────────────────────────

    def _attach_scroll_zoom(self, canvas: FigureCanvasTkAgg) -> None:
        """Bind mouse-wheel events to zoom the first axes of *canvas*.

        Zooming is centred on the cursor position when it is inside the axes;
        otherwise it zooms around the axes centre.  Works alongside the
        NavigationToolbar (zoom-box / pan) without conflict.
        """
        FACTOR = 1.2  # zoom-in multiplier per scroll notch

        def _do_zoom(factor: float, ex: int, ey: int) -> None:
            if not canvas.figure.axes:
                return
            ax = canvas.figure.axes[0]
            # Convert tkinter pixel coords (origin = top-left of widget) to
            # matplotlib display coords (origin = bottom-left of figure).
            h = canvas.figure.bbox.height
            xd, yd = float(ex), h - float(ey)
            if ax.bbox.contains(xd, yd):
                inv = ax.transData.inverted()
                cx, cy = inv.transform((xd, yd))
            else:
                cx = sum(ax.get_xlim()) / 2.0
                cy = sum(ax.get_ylim()) / 2.0
            xl = ax.get_xlim()
            yl = ax.get_ylim()
            ax.set_xlim([cx - (cx - xl[0]) / factor,
                         cx + (xl[1] - cx) / factor])
            ax.set_ylim([cy - (cy - yl[0]) / factor,
                         cy + (yl[1] - cy) / factor])
            canvas.draw_idle()

        w = canvas.get_tk_widget()
        # Windows / macOS: event.delta is ±120 per notch
        w.bind("<MouseWheel>",
               lambda e: _do_zoom(FACTOR if e.delta > 0 else 1 / FACTOR,
                                   e.x, e.y))
        # Linux: Button-4 = scroll up, Button-5 = scroll down
        w.bind("<Button-4>", lambda e: _do_zoom(FACTOR,     e.x, e.y))
        w.bind("<Button-5>", lambda e: _do_zoom(1 / FACTOR, e.x, e.y))

    # ── Export actions ────────────────────────────────────────────────────

    def _on_export_png(self) -> None:
        self._export("png")

    def _on_export_svg(self) -> None:
        self._export("svg")

    def _export(self, fmt: str) -> None:
        out_dir = self.v_out_dir.get()
        if not os.path.exists(os.path.join(out_dir, "nodes.csv")):
            messagebox.showwarning("No data",
                                   "Run Generate first to produce nodes.csv.")
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
                save_kw: dict = dict(
                    bbox_inches="tight",
                    facecolor=fig.get_facecolor(),
                )
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
