#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Second figure of the README: a hump leaving a flat tank through its two open ends.

Panel (a): the free surface left in the tank at t = 8 s, after both halves of the hump have left,
with the characteristic open end and with the extrapolation of the end cell (second-order scheme).
Panel (b): the largest elevation |eta - eta0| in the tank against time, for the second-order and
the first-order schemes with both open ends.
Run open_end_reflection_data.m first (MATLAB, or GNU Octave), then, from this folder,
    python3 plot_open_end_reflection.py
which writes open-end-reflection.png. The style is the stamped copy manuscript_style.py beside it.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

d = loadmat(HERE / "open_end_reflection.mat")
t = d["t"].ravel()
x = d["xc"].ravel()
t_snap = float(d["T_snap"].ravel()[0])

hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 265
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

plain = lambda i: {k: v for k, v in hs.series(i).items() if not k.startswith("marker")}
fig, (ax, bx) = plt.subplots(2, 1, figsize=(hs.WIDTH, 4.6), layout="constrained")

ax.plot(x, 1e3 * d["eta_snap_second_extrapolation"].ravel(), **plain(1), lw=1.1,
        label="extrapolation of the end cell")
ax.plot(x, 1e3 * d["eta_snap_second_characteristic"].ravel(), **plain(0), lw=1.3,
        label="characteristic open end")
ax.set_xlim(0, 20)
ax.set_xlabel(r"$x$ (m)")
ax.set_ylabel(r"$\eta - \eta_0$ (mm)")
hs.panel(ax, rf"(a) free surface left in the tank at $t = {t_snap:g}$ s, second order")
hs.legend_outside(ax, where="above", ncol=2)

every = 25
mask = t > 0
bx.semilogy(t[mask], d["etamax_second_extrapolation"].ravel()[mask], **hs.series(1), markevery=every,
            label="second order, extrapolation")
bx.semilogy(t[mask], d["etamax_second_characteristic"].ravel()[mask], **hs.series(0),
            markevery=(every // 2, every), label="second order, characteristic")
bx.semilogy(t[mask], d["etamax_first_extrapolation"].ravel()[mask], **hs.series(3), markevery=every,
            label="first order, extrapolation")
bx.semilogy(t[mask], d["etamax_first_characteristic"].ravel()[mask], **hs.series(2),
            markevery=(every // 2, every), label="first order, characteristic")
bx.set_xlim(0, float(t[-1]))
bx.set_ylim(1e-5, 0.2)
bx.set_xlabel(r"$t$ (s)")
bx.set_ylabel(r"$\max |\eta - \eta_0|$ (m)")
hs.grid(bx)
hs.panel(bx, "(b) largest elevation in the tank")
hs.legend_outside(bx, where="right", ncol=1)

out = HERE / "open-end-reflection.png"
hs.save(fig, out)
# The background is opaque: without its alpha channel, and compressed a second time without
# loss, the file is smaller.
Image.open(out).convert("RGB").save(out, optimize=True)
print("written", out)
