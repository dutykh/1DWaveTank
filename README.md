<!-- Author: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE) -->

# 1DWaveTank: a numerical wave tank for the nonlinear shallow water equations, by finite volumes

![Four frames of the default experiment of 1DWaveTank: sine waves made at the left end of the tank steepen into bores and run up a beach](docs/figures/wavetank-sloping-beach.png)

The default experiment of the code (`sloping_beach`), computed in MATLAB as `run_simulation.m` computes it, at
t = 6, 12, 18 and 25 s. The tank is 40 m long and 1 m deep at rest, with a beach of slope 0.1 over its last third and a
wall at its right end; at the left end the depth is driven as 1 + 0.2 sin(2πt/4) m. The first two sine waves have
steepened into bores over the flat bottom (a); the first one has run up the beach to 0.2 m above the still-water level
(b); the second bore climbs the slope while the water of the first one runs back down (c); the run ends with three bores
in the tank (d). The computation uses 1000 cells, the MUSCL reconstruction with the van Albada limiter, HLLE fluxes and
the integrator `ode113` of MATLAB with tolerances 10⁻⁵; it takes about 26 s.

The MATLAB code 1DWaveTank solves the nonlinear shallow water (Saint-Venant) equations in one space dimension by
finite volume methods, over a variable bottom that may emerge from the water, with walls, wave-making,
open or periodic ends, and an optional bottom friction. It is a laboratory to compare schemes: each ingredient
(numerical flux, reconstruction, slope limiter, time integrator, boundary condition, initial condition, bottom,
friction law) is a function in its own package folder, chosen by a function handle in the configuration. No article
accompanies the code, so the software itself is cited (see How to cite).

## The problem

With the total depth H(x, t), the discharge HU(x, t) (U is the velocity averaged over the depth), the bottom elevation
z_b(x) and the free surface η = z_b + H, all elevations being counted upwards from the still-water level, the code
solves

```text
∂H/∂t  + ∂(HU)/∂x               = 0,
∂(HU)/∂t + ∂(HU² + g H²/2)/∂x    = − g H ∂z_b/∂x + F,
```

where g = 9.81 m/s² and F is the bottom friction of the chosen law (none by default; Chézy, Manning and
Darcy–Weisbach laws are provided). A cell whose depth is below the dry tolerance 10⁻⁶ m counts as dry.

![Sketch of a wave tank with the notation of long-wave models](img/MyWaveTank.jpg)

The sketch above shows a tank in the notation of long-wave models (free surface y = η(x, t), flat bottom y = −d,
velocity ū averaged over the depth); the code uses a bottom z_b(x) of any shape, the total depth H and the discharge HU.

The tank is cut into N cells of width Δx, and the averages of H and HU over the cells obey a system of ordinary
differential equations: the differences of the numerical fluxes across the cell interfaces, given by an approximate
Riemann solver from the states on both sides, and a discrete bed slope. Three right-hand sides are provided. The
first-order one, `core.rhs_nsw_1st_order`, uses the hydrostatic reconstruction of Audusse, Bouchut, Bristeau, Klein and
Perthame (2004, cited below): before the Riemann solver is called, the depths on both sides of the interface between
cells i and i + 1 are lowered to the higher of the two bottoms, and the bed slope enters the momentum equation of cell i
through these lowered depths (Φ is the momentum component of the numerical flux, z[i] the bottom of cell i):

```text
H_L[i+1/2] = max(0, H[i]   + z[i]   − max(z[i], z[i+1]))       (left state at the interface i + 1/2)
H_R[i+1/2] = max(0, H[i+1] + z[i+1] − max(z[i], z[i+1]))       (right state at the interface i + 1/2)

d(HU)[i]/dt = − (Φ[i+1/2] − Φ[i−1/2])/Δx + g/(2Δx) (H_L[i+1/2]² − H_R[i−1/2]²)
```

so that a lake at rest is an exact steady state over any bottom, next to a wall as well (see below). The hybrid one,
`core.rhs_nsw_hybrid_order`, applies the same lowering to the interface values of a higher-order reconstruction and
keeps the first-order bed slope. The high-order one, `core.rhs_nsw_high_order` (that of the default experiment), passes
the reconstructed interface values to the Riemann solver and discretises the bed slope as − g H[i] (z[i+1/2] −
z[i−1/2])/Δx, the bottom of each interface being the mean of those of its two cells. The system is then integrated in
time either by a solver of MATLAB (`ode113`, `ode45`, ...) or by one of the explicit schemes of the folder `+time`, whose
step obeys the condition Δt = CFL Δx / max(|U| + √(gH)) or an estimate of the local error.

At a solid wall, the ghost cells beyond the end of the tank hold the mirror image of the cells next to it: the depth
and the bottom are reflected as even functions about the wall, the discharge as an odd one. For a wall closing the tank
after its last cell N,

```text
H[N+k] = H[N+1−k],     HU[N+k] = − HU[N+1−k],     z[N+k] = z[N+1−k],     k = 1, 2, ...
```

Data reflected in this way stay symmetric under the shallow water equations, so that the velocity vanishes on the wall
(the argument of LeVeque, Finite Volume Methods for Hyperbolic Problems, Cambridge University Press 2002, Section 7.3.3,
with the bottom reflected like the depth). At the interface on the wall, the hydrostatic reconstruction then meets the
same bottom on both sides and the same depth H[N], so that no water crosses the wall and a lake at rest stays at rest
next to it, over any bottom. The first-order and hybrid right-hand sides reflect the bottom in this way (function
`bc.wall_bathymetry`); the high-order one takes the bottom of its ghost cells from the bottom function beyond the wall
(see A check of the schemes).

## Quick start

In MATLAB, from the root folder of the repository:

```matlab
run_simulation
```

The script adds the repository to the path, builds the configuration of the experiment selected in
`+cfg/simulation_config.m` (by default `sloping_beach`), prints it, solves the equations with `core.solver`, animates the
free surface, the bottom and the velocity at every output time, and prints the CPU time (about 26 s for the default
experiment, with MATLAB R2025b on one core). The solution stays in the workspace: `results.t` holds the output times and
`results.H`, `results.HU` and `results.U` one row per output time and one column per cell, and `config` holds the
configuration. With `config.save_results = true` in the configuration, both are also saved to a file
`output/results_<date>_<time>.mat`.

Another experiment is chosen by editing the line

```matlab
experiment_setup = 'sloping_beach'; % CHANGE THIS TO SELECT SETUP
```

near the top of `+cfg/simulation_config.m`; the parameters of each experiment are set in its own `case` of that file.
The user guide `UserGuide.md` describes the configuration fields and every component that can be chosen.

## Experiments

The nine experiments are the cases of `+cfg/simulation_config.m`; each case completes the defaults of
`+cfg/default_config.m` and the common settings written at the top of the file, which select a tank of 20 m with 500
cells, the first-order scheme with PVM fluxes and `ode113` with tolerances 10⁻⁴, up to t = 10 s. All nine run to their
final time. The last column reports a run of each in GNU Octave 11.1 (with `ode45` at the same tolerances where an
experiment asks for `ode113`, and with a stand-in for the function `contains` of MATLAB, which the MUSCL reconstruction
calls with a limiter of two arguments; Octave has neither) and, for the closed tanks, the largest change of the volume of
water V(t) = Σ H[i] Δx during the run, relative to its initial value.

| Experiment | Set-up | Run |
|---|---|---|
| `sloping_beach` | the default experiment, shown above | to 25 s; the same solution as with version 1.1.1, to the last bit, since the statements it executes are unchanged |
| `flat_wave_gen` | sine waves of amplitude 0.1 m and period 2π s made at the left end of a 20 m tank of depth 0.5 m, wall at the right end; MUSCL with the van Leer limiter, PVM fluxes, the second-order strong-stability-preserving Runge–Kutta scheme (SSP2) at CFL = 0.95 | to 10 s; the same solution as with version 1.1.1, to the last bit |
| `flat_rest` | lake at rest 0.5 m deep over a flat bottom, walls at both ends; first order | to 10 s; the water stays at rest, the discharge exactly zero and the depth 0.5 m to 10⁻¹⁵; V constant to 2 × 10⁻¹⁶ |
| `gaussian_bump_rest` | lake at rest over a Gaussian bump of height 0.05 m, walls; hybrid scheme with MUSCL and the van Leer limiter, HLLC fluxes, `ode113` with tolerances 10⁻⁹ | to 10 s, but the tank holds no water: `ic.lake_at_rest` puts the free surface at z = 0, the level of the bottom on which the bump stands, so that every cell has the dry depth 10⁻⁶ m |
| `flat_gaussian` | hump 0.1 exp(−2 (x − 10)²) m of the free surface, at rest on a depth of 0.5 m, open ends; MUSCL with the van Leer limiter, PVM fluxes, SSP2 at CFL = 0.95 | to 10 s; the two halves of the hump leave through the open ends, and V falls from 10.125 to 9.996 m² |
| `flat_solitary` | solitary wave of amplitude 0.2 m on a depth of 0.5 m in the middle of the tank, walls; first order | to 10 s, the wave reaching the right wall at t = 3.3 s; V constant to 3 × 10⁻¹⁵ |
| `periodic_solitary` | the same wave in a periodic tank; first order | to 10 s; V constant to 2 × 10⁻¹⁵ |
| `dam_break` | breaking of a dam at x = 10 m between depths 0.8 m and 0.5 m, walls; first order, SSP2 at CFL = 0.95 | to 3 s; V constant to 4 × 10⁻¹⁵ |
| `dry_dam_break` | breaking of a dam at x = 10 m holding a depth of 0.5 m, dry bottom beyond it, walls; MUSCL in the Riemann invariants u ∓ 2√(gH) with the van Leer limiter, HLLE fluxes, forward Euler steps at CFL = 0.45 | to 1 s; V constant to 5 × 10⁻¹⁵, the depth never negative, the front at x = 14.1 m where the exact solution of Ritter puts it at 10 + 2√(g × 0.5) = 14.4 m |

The solitary wave is that of the Serre–Green–Naghdi equations (Carter and Cienfuegos, European Journal of Mechanics
B/Fluids 30, 2011, equation (9)); the shallow water equations, which have no dispersion, steepen it into a bore. Ritter's
solution of the dam break over a dry bottom is derived in Castro-Orgaz and Hager, Shallow Water Hydraulics, Springer
2019, Section 6.2: the front moves at the speed 2√(g h_L), h_L being the depth behind the dam.

Version 1.2.0 makes every experiment run. Four initial conditions, `ic.gaussian_bump`, `ic.solitary_wave`,
`ic.dam_break` and `ic.dry_dam_break`, took the cell centres as their first argument, while `core.solver` calls every
initial condition as `ic_handle(cfg)`; they now take `cfg` alone, as `ic.lake_at_rest` does, and read the cell centres
in `cfg.mesh.xc`. The experiments `flat_gaussian` and `periodic_solitary` read the field `config.param`, which the
common settings no longer define; they now set their still depth themselves (and `periodic_solitary` the amplitude of
its wave). The settings of `dry_dam_break` gave its MUSCL reconstruction neither a method nor an order, which
`cfg.validate_config` needs, and did not pass it to the high-order right-hand side, which would then have computed at
first order; both are set now, and its CFL number is 0.45 instead of 0.95, because forward Euler steps with a
second-order reconstruction keep the depth nonnegative under half the CFL bound of the first-order scheme (Audusse et
al. 2004, Section 3): at 0.95 the depth fell to −1.5 cm at the front. Four experiments, `flat_rest`, `flat_gaussian`,
`flat_solitary` and `periodic_solitary`, never set the ends of their mesh, which `cfg.validate_config` builds from
`mesh.domain`, and so were set in the 100 m tank of the defaults (where `flat_rest` ran) instead of the 20 m tank of
their names and initial conditions; they now copy the ends of the tank into `mesh.domain`. Finally, `flat_rest` sets its
flat bottom 0.5 m below the still-water level, so that its tank holds the 0.5 m of water of its name.

The solid wall of version 1.2.0 reflects the bottom as well (see The problem). The bottom of its ghost cells used to be
the bottom function evaluated beyond the wall, so that over a sloping bottom the two states of the hydrostatic
reconstruction at the wall differed in depth by about |dz/dx| Δx: a lake at rest was set in motion from the walls, and
water crossed them, at a rate that does not decrease with Δx. In a tank of 10 m closed by two walls over the bottom
z = −1 + 0.05 x, the first-order scheme gave the depth in the two cells next to the walls the time derivatives 0.078
and 0.055 m/s at rest, with 100, 200 or 400 cells alike; in 4000 steps of SSP2 at CFL = 0.5 (64 s, 100 cells), the
velocity reached 0.021 m/s and 2 % of the water left the tank. With the reflected bottom the velocity stays below
5 × 10⁻¹⁵ m/s and the volume does not change (A check of the schemes). Over a flat bottom the reflected and the former
ghost bottoms coincide, so that a wave reflected from a wall over a flat bottom is computed exactly as before, to the
last bit: `flat_wave_gen` is, and so is a hump of the free surface released between two walls over a flat bottom and
reflected back and forth for 20 s, by each of the three right-hand sides.

Version 1.1.1 corrected a misprint at line 84 of `+cfg/simulation_config.m`, where the common settings wrote the empty
reconstruction into the field `reconstructopenion` instead of `reconstruction`, the field that `cfg.validate_config`
reads. Through the misprint, the characteristic MUSCL reconstruction that `+cfg/default_config.m` defines as its default
reached every experiment that does not replace it whole: the configuration then announced a MUSCL reconstruction of
order 2 for experiments computed without one, and `gaussian_bump_rest` reconstructed in characteristic variables, which
its own settings do not ask for. The default experiment and `flat_wave_gen` gave the same solutions as before to the
last digit.

## A check of the schemes

A scheme is well balanced when a lake at rest (U = 0, and η = 0 wherever the bottom lies below the still-water level) is
an exact steady state of its discrete equations. The first-order scheme is well balanced, walls included. Three bottoms
were tried on 0 ⩽ x ⩽ 10 m, closed by walls at both ends and cut into 100 cells: the slope z = −1 + 0.05 x, the wavy
bottom z = −0.7 + 0.15 sin(2πx/7) + 0.02 x, and the beach z = −1 + 0.15 x, whose shoreline lies 3.3 m before the right
wall. In 4000 steps of SSP2 at CFL = 0.5 (64 to 74 s), the largest velocity stays below 6 × 10⁻¹⁵ m/s with each of the
fluxes HLL, HLLC, HLLE, Rusanov, Roe, FVCF, Lax–Friedrichs and central-upwind, the free surface within 2 × 10⁻¹⁵ m of
rest and the volume constant to 3 × 10⁻¹⁵. With PVM fluxes the velocity reaches 10⁻⁶ m/s: the PVM flux adds no
dissipation where the states on the two sides of an interface differ by less than 10⁻⁶ (its threshold is the dry
tolerance), so that rounding errors grow up to that size.

The hybrid and high-order schemes reconstruct the depth and the discharge, not the free surface H + z, and keep a lake
at rest only in part. Over the slope above, away from the walls, the time derivative of the depth at rest is
2.7 × 10⁻⁴ m/s in the hybrid scheme, halved with each halving of Δx (its reconstructed depths are lowered with the
bottoms of the cell centres), and below 3 × 10⁻⁹ m/s in the high-order one, while their momentum equations balance to
round-off. Next to a wall over a slope s, the limiter flattens the reconstructed depth in the last cell (as it does at
a kink of the bottom or at a shoreline), the two states at the inner interface of that cell differ by s Δx/2, and its
depth changes at a rate that does not decrease with Δx,

```text
dH/dt ≈ c s / 4,     c = √(g H),
```

0.039 m/s next to the left wall of the slope above. In 64 s the velocity reaches 0.024 m/s in the hybrid scheme and
0.011 m/s in the high-order one over that slope, and up to 0.8 m/s at the shoreline of the beach. The reflected bottom
stops the losses of water through the walls in the hybrid scheme (2 % of the water in 64 s over the slope before, none
now) and divides its imbalance in the momentum equation next to a wall by three. The high-order scheme loses no water
through its walls with either bottom; it keeps the bottom function beyond the wall, which leaves the default experiment
unchanged: its bores reach the right wall over the beach from about t = 12.4 s to 16.8 s and again after t = 23.6 s,
where the reflected bottom would change the bed slope of the last cell without making the lake at rest steady for this
scheme. Over the bottom of the default experiment, the hybrid and high-order right-hand sides at rest also reach 0.22
m²/s² in the momentum equation in the three cells at the foot of the beach, where the bottom has a kink. The
second-order scheme of Audusse et al. (2004, Section 3) reconstructs the free surface together with the depth for this
reason.

## Contents

| Folder or file | Contents |
|---|---|
| `run_simulation.m` | the main script: configuration, solution, animation and statistics |
| `+cfg/` | `default_config.m` (defaults), `simulation_config.m` (the experiments), `validate_config.m` (checks and completes the configuration) |
| `+core/` | `solver.m` (initial condition, time integration, output) and the right-hand sides `rhs_nsw_1st_order.m`, `rhs_nsw_hybrid_order.m` and `rhs_nsw_high_order.m`; the folder `multilayer/` holds a separate solver of the multilayer shallow water equations for atmospheric flows (layers of air with the densities of the standard atmosphere), with its own copies of these components |
| `+flux/` | sixteen numerical fluxes: AUSM+ (`AUSM.m`), AUSMDV, central-upwind (`CentralUpwind.m`), FORCE, FVCF, HLL, HLLC, HLLE, kinetic (`Kinetic.m`), Lax–Friedrichs (`LF.m`), Osher–Solomon, PVM, Roe, Rusanov, SLAU and Steger–Warming |
| `+reconstruct/` | reconstructions of the interface values: none (first order), MUSCL, ENO2, UNO2, PPM, WENO5, MP5 and THINC, named in `reconstruct_selector.m`; the folder `+limiters/` holds the slope limiters of MUSCL: minmod, van Leer, van Albada, superbee, MC, Koren, UMIST, OSPRE and Sweby |
| `+time/` | explicit integrators with a step set by the CFL condition (forward Euler, the strong-stability-preserving Runge–Kutta schemes SSP2 and SSP3, the classical fourth-order Runge–Kutta scheme, the two-step Adams–Bashforth scheme) or by an estimate of the local error (Bogacki–Shampine 3(2), Dormand–Prince 5(4)), and `integrate_matlab_ode.m` for the solvers of MATLAB |
| `+bc/` | boundary conditions: `wall`, `open` (zero gradient), `periodic`, and `generating` (a depth prescribed in time, with the outgoing Riemann invariant taken from the tank); `wall_bathymetry` gives the ghost cells of a wall the mirror image of the bottom |
| `+ic/` | initial conditions, each called as `ic(cfg)`: `lake_at_rest`, `gaussian_bump`, `solitary_wave`, `dam_break` and `dry_dam_break` |
| `+bathy/` | bottoms: `flat`, `gaussian_bump` and `sloping_beach` |
| `+friction/` | friction laws: none, Chézy, Manning, and Darcy–Weisbach with a constant factor or the Colebrook–White factor, named in `friction_selector.m` |
| `+utils/` | the CFL time step, the physical flux, Roe averages, eigenvectors, a uniform mesh and a text progress bar |
| `+vis/` | `plot_state.m`, one frame of the animation |
| `docs/figures/` | the two scripts of the figure at the top of this page |
| `img/` | the sketch of the tank |

The figure is made by `docs/figures/wavetank_readme_data.m` (MATLAB, or GNU Octave), which solves the default
experiment as `run_simulation.m` does and writes every output frame to `wavetank_frames.mat`, and by
`docs/figures/plot_readme_figure.py` (Python 3 with NumPy, SciPy and Matplotlib, and the style file beside it), which
draws it. Both are run from `docs/figures/`, in this order.

## Requirements

MATLAB, with no toolbox (checked with R2025b). GNU Octave cannot run the default experiment as it stands, since it has
no `ode113`; the MUSCL reconstruction with a limiter of two arguments also calls `contains`, which Octave lacks. The
script `docs/figures/wavetank_readme_data.m` runs in both: without `ode113` it integrates with `ode45` at the same
tolerances, about thirty times more slowly in GNU Octave 11.1, whose depths over the first 4 s agree with those of
MATLAB to 2 × 10⁻⁵ m.

## How to cite

No article accompanies this code, so please cite the software itself:

D. Dutykh, 1DWaveTank: a numerical wave tank for the nonlinear shallow water equations, by finite volumes, version
1.2.0 (2026).

https://github.com/dutykh/1DWaveTank

```bibtex
@software{Dutykh_1DWaveTank,
  author = {Dutykh, Denys},
  title = {{1DWaveTank}: a numerical wave tank for the nonlinear shallow water equations, by finite volumes},
  year = {2026},
  version = {1.2.0},
  url = {https://github.com/dutykh/1DWaveTank},
  license = {GPL-3.0-or-later}
}
```

The button "Cite this repository" in the right-hand column of this page gives the same reference, read from
`CITATION.cff`. The hydrostatic reconstruction of the first-order and hybrid schemes is that of the article

E. Audusse, F. Bouchut, M.-O. Bristeau, R. Klein and B. Perthame, A fast and stable well-balanced scheme with
hydrostatic reconstruction for shallow water flows, SIAM Journal on Scientific Computing 25, 2050–2065 (2004).

https://doi.org/10.1137/s1064827503431090

## Licence

GNU General Public License, version 3 or (at your option) any later version (see `LICENSE`). The functions `atmos.m` and
`pressurealt.m` of `+core/multilayer/` are by Sky Sartorius (MATLAB Central File Exchange, 2015), and the text progress
bar of `+utils/textprogressbar.m` is adapted from the one of Paul Proteus; they remain the work of their authors and are
not covered by the licence of this repository.

## Authors

Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE).

https://www.denys-dutykh.com/

Mehmet Ersoy (SEATECH and Institut de Mathématiques de Toulon, Université de Toulon) contributed the kinetic flux
`+flux/Kinetic.m`, written from his original Fortran code, and Francesco Carbone (Institute of Atmospheric Pollution
Research of the National Research Council of Italy, Rende) contributed the multilayer solver of `+core/multilayer/` and
pointed out several errors in the code. Questions, corrections and contributions are welcome through the issues and
pull requests of this repository.
