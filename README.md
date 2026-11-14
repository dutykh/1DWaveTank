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

so that, away from the ends of the tank, a lake at rest is an exact steady state over any bottom. The hybrid one, `core.rhs_nsw_hybrid_order`, applies
the same lowering to the interface values of a higher-order reconstruction and keeps the first-order bed slope. The
high-order one, `core.rhs_nsw_high_order` (that of the default experiment), passes the reconstructed interface values
to the Riemann solver and discretises the bed slope as − g H[i] (z[i+1/2] − z[i−1/2])/Δx, the bottom of each interface
being the mean of those of its two cells. The system is then integrated in time either by a solver of MATLAB
(`ode113`, `ode45`, ...) or by one of the explicit schemes of the folder `+time`, whose step obeys the condition
Δt = CFL Δx / max(|U| + √(gH)) or an estimate of the local error.

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
`+cfg/default_config.m` and the common settings written at the top of the file, which select the first-order scheme
with PVM fluxes and `ode113` with tolerances 10⁻⁴, up to t = 10 s. Their state as shipped was checked in MATLAB R2025b:

| Experiment | Set-up | As shipped |
|---|---|---|
| `sloping_beach` | the default experiment, shown above | runs |
| `flat_wave_gen` | sine waves of amplitude 0.1 m and period 2π s made at the left end of a 20 m tank of depth 0.5 m with 500 cells, wall at the right end; MUSCL with the van Leer limiter, PVM fluxes, the second-order strong-stability-preserving Runge–Kutta scheme at CFL = 0.95; 10 s | runs, in about 60 s |
| `flat_rest` | lake at rest over a flat bottom, walls at both ends; first order | runs, but the tank holds no water: `ic.lake_at_rest` puts the free surface at z = 0, the level of this bottom, so that every cell has the dry depth 10⁻⁶ m |
| `gaussian_bump_rest` | lake at rest over a Gaussian bump of height 0.05 m, walls; hybrid scheme with MUSCL and the van Leer limiter, HLLC fluxes, `ode113` with tolerances 10⁻⁹ | runs, but the tank holds no water, for the same reason (the bump stands on a bottom at z = 0) |
| `flat_gaussian`, `periodic_solitary` | a Gaussian hump of the free surface between open ends; a solitary wave in a periodic tank | stop while the configuration is built: they read the field `config.param`, which the common settings no longer define |
| `flat_solitary`, `dam_break` | a solitary wave between walls; the breaking of a dam between depths 0.8 m and 0.5 m | stop at the initial condition: `ic.solitary_wave` and `ic.dam_break` take the arguments `(xc, cfg)`, while `core.solver` calls every initial condition with `cfg` alone |
| `dry_dam_break` | the breaking of a dam of depth 0.5 m on a dry bottom | stops in `cfg.validate_config`, since it sets a reconstruction handle but neither a method nor an order (its initial condition has the same mismatch of arguments as above) |

Version 1.1.1 corrects a misprint at line 84 of `+cfg/simulation_config.m`, where the common settings wrote the empty
reconstruction into the field `reconstructopenion` instead of `reconstruction`, the field that `cfg.validate_config`
reads. Through the misprint, the characteristic MUSCL reconstruction that `+cfg/default_config.m` defines as its default
reached every experiment that does not replace it whole: the configuration then announced a MUSCL reconstruction of
order 2 for experiments computed without one, and `gaussian_bump_rest` reconstructed in characteristic variables, which
its own settings do not ask for. The default experiment and `flat_wave_gen` give the same solutions as before to the
last digit; `dry_dam_break`, which stopped at its initial condition, now stops earlier, in `cfg.validate_config`.

## A check of the schemes

A scheme is well balanced when a lake at rest (U = 0, and η = 0 wherever the bottom lies below the still-water level) is
an exact steady state of its discrete equations. At such a state over the bottom of the default experiment, with walls
at both ends, the right-hand side of the first-order scheme vanishes in every wet cell to round-off: its largest value,
3 × 10⁻¹⁴ m²/s², is in the momentum equation. Those of the hybrid and high-order schemes stay below 10⁻¹⁰ m²/s² on the
flat bottom and on the slope, but not in the three cells at the foot of the beach, where the bottom has a kink (0.22 m²/s² in
the momentum equation), nor in the last wet cell at the shoreline.

## Contents

| Folder or file | Contents |
|---|---|
| `run_simulation.m` | the main script: configuration, solution, animation and statistics |
| `+cfg/` | `default_config.m` (defaults), `simulation_config.m` (the experiments), `validate_config.m` (checks and completes the configuration) |
| `+core/` | `solver.m` (initial condition, time integration, output) and the right-hand sides `rhs_nsw_1st_order.m`, `rhs_nsw_hybrid_order.m` and `rhs_nsw_high_order.m`; the folder `multilayer/` holds a separate solver of the multilayer shallow water equations for atmospheric flows (layers of air with the densities of the standard atmosphere), with its own copies of these components |
| `+flux/` | sixteen numerical fluxes: AUSM+ (`AUSM.m`), AUSMDV, central-upwind (`CentralUpwind.m`), FORCE, FVCF, HLL, HLLC, HLLE, kinetic (`Kinetic.m`), Lax–Friedrichs (`LF.m`), Osher–Solomon, PVM, Roe, Rusanov, SLAU and Steger–Warming |
| `+reconstruct/` | reconstructions of the interface values: none (first order), MUSCL, ENO2, UNO2, PPM, WENO5, MP5 and THINC, named in `reconstruct_selector.m`; the folder `+limiters/` holds the slope limiters of MUSCL: minmod, van Leer, van Albada, superbee, MC, Koren, UMIST, OSPRE and Sweby |
| `+time/` | explicit integrators with a step set by the CFL condition (forward Euler, the strong-stability-preserving Runge–Kutta schemes SSP2 and SSP3, the classical fourth-order Runge–Kutta scheme, the two-step Adams–Bashforth scheme) or by an estimate of the local error (Bogacki–Shampine 3(2), Dormand–Prince 5(4)), and `integrate_matlab_ode.m` for the solvers of MATLAB |
| `+bc/` | boundary conditions: `wall`, `open` (zero gradient), `periodic`, and `generating` (a depth prescribed in time, with the outgoing Riemann invariant taken from the tank) |
| `+ic/` | initial conditions: `lake_at_rest`, `gaussian_bump`, `solitary_wave`, `dam_break` and `dry_dam_break` |
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
1.1.1 (2026).

https://github.com/dutykh/1DWaveTank

```bibtex
@software{Dutykh_1DWaveTank,
  author = {Dutykh, Denys},
  title = {{1DWaveTank}: a numerical wave tank for the nonlinear shallow water equations, by finite volumes},
  year = {2026},
  version = {1.1.1},
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
