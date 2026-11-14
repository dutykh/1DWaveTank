#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: four frames of the default experiment of 1DWaveTank (sloping beach).

Run wavetank_readme_data.m first (MATLAB, or GNU Octave), then, from this folder,
    python3 plot_readme_figure.py
which writes wavetank-sloping-beach.png. The style is the stamped copy manuscript_style.py beside it.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

TIMES = (6.0, 12.0, 18.0, 25.0)       # instants shown, among the stored frames [s]
WET = 1e-3                            # depth above which a cell is drawn as water [m]

data = loadmat(HERE / "wavetank_frames.mat")
t = data["t"].ravel()
x = data["xc"].ravel()
zb = data["zb"].ravel()
H = data["H"]
frames = [int(np.argmin(np.abs(t - s))) for s in TIMES]
dx = x[1] - x[0]
edges = (x[0] - dx / 2, x[-1] + dx / 2)

hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 265
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig, axes = plt.subplots(len(TIMES), 1, figsize=(hs.WIDTH, 4.6), sharex=True)
water = matplotlib.colors.to_rgba(hs.NAVY, 0.12)
for j, (axis, k) in enumerate(zip(axes, frames)):
    depth = H[k]
    surface = np.where(depth > WET, zb + depth, np.nan)
    axis.fill_between(x, zb, zb + depth, where=depth > WET, facecolor=water, edgecolor="none",
                      interpolate=False)
    axis.fill_between(x, -1.25, zb, facecolor=hs.SAND, edgecolor="none")
    axis.plot(x, zb, color=hs.INK, lw=0.8)
    sea = zb < 0.0
    axis.plot(x[sea], np.zeros(sea.sum()), color=hs.SLATE, lw=0.6, linestyle=(0, (1.6, 1.6)))
    axis.plot(x, surface, color=hs.NAVY, lw=1.1)
    axis.plot([edges[1], edges[1]], [zb[-1], 0.62], color=hs.INK, lw=1.6, solid_capstyle="butt")
    axis.set_xlim(edges[0], edges[1] + 0.25)
    axis.set_ylim(-1.25, 0.62)
    axis.set_yticks([-1.0, -0.5, 0.0, 0.5])
    hs.panel(axis, rf"({'abcdef'[j]}) $t = {t[k]:g}$ s")
axes[-1].set_xlabel(r"$x$ (m)")
fig.supylabel(r"elevation $z$ (m)", fontsize=plt.rcParams["axes.labelsize"])
axes[-1].set_xticks(np.arange(0, 41, 5))

handles = [Line2D([], [], color=hs.NAVY, lw=1.1),
           Line2D([], [], color=hs.SLATE, lw=0.6, linestyle=(0, (1.6, 1.6))),
           Patch(facecolor=hs.SAND, edgecolor=hs.INK, linewidth=0.8)]
labels = [r"free surface $\eta = z_b + H$", r"still-water level $z = 0$", r"bottom $z_b$"]
hs.legend_outside(fig, where="upper", ncol=3, handles=handles, labels=labels)

out = HERE / "wavetank-sloping-beach.png"
hs.save(fig, out)
# The background is opaque: without its alpha channel, and compressed a second time without
# loss, the file is smaller.
Image.open(out).convert("RGB").save(out, optimize=True)
print("written", out, "frames at t =", [float(t[k]) for k in frames])
