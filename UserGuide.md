<!-- Author: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE) -->

# 1DWaveTank: user guide

This guide describes how a run of 1DWaveTank proceeds, the fields of its configuration and the calling conventions of
its components, so that a component can be chosen or added. The equations, the schemes, the experiments and their state
are described in the `README.md` of the repository.

## How a run proceeds

The script `run_simulation.m` adds the repository to the path and calls `cfg.simulation_config()`. That function loads
the defaults of `cfg.default_config()`, applies the common settings written at its top, then the `case` of the
experiment named in its variable `experiment_setup`, and finally passes the structure to `cfg.validate_config`, which
checks it and completes it (the mesh, the output times, the number of ghost cells, the reconstruction settings). The
script then calls `core.solver(config)`, which builds the initial state with `config.ic_handle(config)`, integrates the
semi-discrete equations given by `config.model` with the integrator `config.timeStepper` over the output times
`config.tspan`, and returns the structure `results` (fields `t`, `H`, `HU`, `U`, `cfg` and `total_steps`). The script
finally draws every output frame with `vis.plot_state` and, when `config.save_results` is true, saves `results` and
`config` to `config.outputPath`.

## The configuration

The fields below are those that the components read. A field that an experiment does not set keeps its value from
`+cfg/default_config.m` or from the common settings of `+cfg/simulation_config.m`.

| Field | Meaning |
|---|---|
| `mesh.domain.xmin`, `mesh.domain.xmax`, `mesh.N` | ends of the tank and number of cells, from which `cfg.validate_config` builds the cell width `mesh.dx` and the cell centres `mesh.xc` when they are not set |
| `domain.xmin`, `domain.xmax` | ends of the tank printed in the summary of the configuration; most experiments copy them into `mesh.domain`, from which the mesh is built (`sloping_beach` builds its cell centres from them directly, and `gaussian_bump_rest` sets `mesh.domain` itself), and `ic.solitary_wave` and `bathy.gaussian_bump` read them |
| `phys.g`, `phys.dry_tolerance` | gravity (9.81 m/s²) and the depth below which a cell counts as dry (10⁻⁶ m) |
| `model` | the right-hand side: `@core.rhs_nsw_1st_order`, `@core.rhs_nsw_hybrid_order` or `@core.rhs_nsw_high_order` |
| `numFlux` | the numerical flux, a handle to a function of `+flux` (the PVM flux also reads `pvm_degree`) |
| `reconstruct` | the reconstruction used by the hybrid and high-order right-hand sides: `handle` (a function of `+reconstruct`), `limiter` (a handle or a name for MUSCL), `characteristic` (true for MUSCL, WENO5 or THINC in characteristic variables), `ppm_mode` and `mp5_mode` (`'component'` or `'characteristic'`), `thinc_beta`. These right-hand sides reconstruct the depth, the free surface and the discharge component by component, and the velocity instead of the discharge where a reconstruction in characteristic variables is asked for, since characteristic variables are not flat in a lake at rest over a sloping bottom |
| `reconstruction` | the reconstruction settings that `cfg.validate_config` completes (`method`, `order`, number of ghost cells `ng`, `handle`); the experiments `sloping_beach`, `gaussian_bump_rest` and `dry_dam_break` copy it into `reconstruct` |
| `bc.left.handle`, `bc.right.handle` | the boundary conditions, handles to functions of `+bc` |
| `bc.left.param`, `bc.right.param` | their parameters: the amplitude `a` and the period `T` of `@bc.generating` |
| `param.H0` | the still depth about which `@bc.generating` oscillates, which should be the depth of the tank at its wave-making end at rest, since the ghost cells there take the bottom of the end cell |
| `bc.num_ghost_cells` | ghost cells on each side, used by the hybrid and high-order right-hand sides (the first-order one always uses one); WENO5 stops with fewer than 3, and the header of `muscl.m` asks for at least 2 |
| `ic_handle`, `ic_param` | the initial condition and its parameters (`a`, `lambda`, `x0` and `H0` for `ic.gaussian_bump`; `h0` and `a` for `ic.solitary_wave`) |
| `bathyHandle` | the bottom, a handle to a function of `+bathy`, with its parameters `bathy.params` (`sloping_beach`: `L`, `h0`, `slope`), `bathy_params.flat_elevation` (`flat`), or `bathy_bump_center`, `bathy_bump_height`, `bathy_bump_width`, `h0` and `L` (`gaussian_bump`) |
| `phys.friction_model` | the friction law, `@friction.no_friction` by default, with `phys.chezy_C`, `phys.manning_n`, `phys.darcy_f`, `phys.f_calculation` (`'constant'` or `'colebrook_white'`), `phys.ks`, `phys.kinematic_viscosity`, `phys.cw_iterations` and `phys.cw_tolerance` |
| `timeStepper` | the integrator, a handle to a function of `+time` |
| `time.matlab_solver`, `time.AbsTol`, `time.RelTol`, `time.ode_options`, `time.show_progress_bar` | the solver of MATLAB used by `@time.integrate_matlab_ode` (`'ode113'` by default), its tolerances (10⁻⁴ by default) and options |
| `time.cfl`, `time.adaptive_tol`, `time.num_progress_reports` | the CFL number of the integrators with a CFL step, the tolerance of the Bogacki–Shampine and Dormand–Prince integrators, the number of progress reports |
| `t0`, `tEnd`, `tspan` | start, end and output times of the run |
| `vis.plot_velocity`, `vis.show_legend` | whether the animation shows the velocity in a lower panel, and a legend |
| `save_results`, `outputPath` | whether `run_simulation.m` saves the results, and where (`./output/` by default) |

## The components

Every component is a function in a package folder and is selected by a handle in the configuration. A new component
is added by writing a function with the same arguments in the same folder and selecting its handle. In the conventions
below, N is the number of cells, ng the number of ghost cells, and a state holds the depth H in its first column (or
first half) and the discharge HU in its second.

| Folder | Calling convention | Functions provided |
|---|---|---|
| `+core` | `dwdt = rhs(t, w, cfg)`, with `w` and `dwdt` columns of length 2N, `[H; HU]` | `rhs_nsw_1st_order`, `rhs_nsw_hybrid_order`, `rhs_nsw_high_order`; and `[wL, wR, S] = hydrostatic_reconstruction_2nd_order(w_padded, z_padded, cfg)`, which the last two call: the (N + 1) × 2 states on both sides of each interface after the second-order hydrostatic reconstruction, and the bed slope of each cell |
| `+flux` | `F = flux(wL, wR, cfg)`, with the left and right states at the N + 1 interfaces as (N + 1) × 2 arrays, and the fluxes of mass and momentum in the same layout | `AUSM` (AUSM+), `AUSMDV`, `CentralUpwind`, `FORCE`, `FVCF`, `HLL`, `HLLC`, `HLLE`, `Kinetic`, `LF` (Lax–Friedrichs), `OsherSolomon`, `PVM`, `Roe`, `Rusanov`, `SLAU`, `StegerWarming` |
| `+reconstruct` | `[wL, wR] = method(w_padded, cfg)`, from the (N + 2ng) × 2 array of the cells and ghost cells to the (N + 1) × 2 states on both sides of each interface | `none`, `muscl`, `eno2`, `uno2`, `ppm`, `weno5`, `mp5`, `thinc`, named in `reconstruct_selector` |
| `+reconstruct/+limiters` | `slope = limiter(delta_minus, delta_plus)`, the limited slope from the two one-sided differences (`vanalbada` and `umist` also accept `cfg`, `sweby` a parameter) | `minmod`, `vanleer`, `vanalbada`, `superbee`, `mc`, `koren`, `umist`, `ospre`, `sweby`, named in `limiter_selector` |
| `+time` | `[sol_out, t_out, stats] = integrator(rhs, tspan, w0, cfg)`, with one row of `sol_out` per output time and `stats.nsteps` the number of steps | `integrate_euler_adaptive`, `integrate_ssp2_adaptive`, `integrate_ssp3_adaptive`, `integrate_rk4_adaptive`, `integrate_ab2_adaptive` (step from the CFL condition), `integrate_bogacki_shampine`, `integrate_dopri54` (step from an estimate of the local error), `integrate_matlab_ode` (the solvers of MATLAB) |
| `+bc` | `w_padded = bc(w_padded, t, side, cfg, ng)`, filling the ghost cells of the side `'left'` or `'right'` (a periodic tank is filled at once, with `'both'`) | `wall`, `open` (which gives its ghost cells the free surface and the discharge of the end cell over the bottom function `cfg.bathyHandle`), `periodic`, `generating`; and `z_padded = ghost_bathymetry(z_padded, cfg, ng)`, which gives the ghost cells the bottom that matches their state, for the three right-hand sides: the mirror image of the bottom at a wall (`wall_bathymetry`, with the same arguments), the bottom of the end cell at a wave-making end, the periodic image in a periodic tank, the bottom function at an open end |
| `+ic` | `w0 = ic(cfg)`, returning `[H; HU]`, an N × 2 array or the depth alone, at the cell centres `cfg.mesh.xc` | `lake_at_rest`, `gaussian_bump`, `solitary_wave`, `dam_break`, `dry_dam_break` |
| `+bathy` | `zb = bathy(cfg, x)`, the bottom elevation at the points `x` | `flat`, `gaussian_bump`, `sloping_beach` |
| `+friction` | `F = law(H, HU, g, cfg)`, the term added to the time derivative of the discharge in the wet cells | `no_friction`, `chezy`, `manning`, `darcy_weisbach` (with `colebrook_white` for its factor), named in `friction_selector` |
| `+vis` | `fig = plot_state(xc, H, zb, U, t, cfg, fig, x_limits, y_limits, u_limits)`, one frame of the animation | `plot_state` |

The integrators of `+time` call the right-hand side as `rhs(t, w, cfg)`; for `integrate_matlab_ode`, `core.solver`
wraps it as `rhs(t, w)`, the form that the solvers of MATLAB expect.

## The multilayer solver

The folder `+core/multilayer/` holds a separate code, contributed by Francesco Carbone: the function
`multilayer_well_balanced` (no argument) sets up and solves the multilayer shallow water equations for atmospheric
flows, with its own copies of the solver, of the high-order right-hand side, of a fourth-order Runge–Kutta integrator and
of a Gaussian bottom; the densities of its layers are those of the 1976 standard atmosphere, computed by `atmos.m` of Sky
Sartorius (whose `pressurealt.m` converts pressures into altitudes).
