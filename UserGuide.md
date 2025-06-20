# 1DWaveTank: User Guide

![1DWaveTank Simulation](img/MyWaveTank.jpg)

**Author:**
- Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)

**Contributors:**
- Prof. Mehmet ERSOY (SEATECH - École d'Ingénieurs de l'Université de Toulon, IMATH - Institut de Mathématiques de Toulon, France). The kinetic flux routine (+flux/Kinetic.m) is based on his original Fortran code and we gratefully acknowledge his scientific input and generosity in sharing this algorithm. Contact: http://ersoy.univ-tln.fr/
- Dr. Francesco Carbone (University of Calabria, Italy). We gratefully acknowledge his contribution of the multilayer nonlinear shallow water equations implementation for atmospheric modeling, which has been included in this distribution.

**Date:** May 15, 2025

## 1. Overview

**Latest Well-Balanced Functionality:**
- The code now features a robust well-balanced hydrostatic reconstruction for all high-order finite volume schemes (MUSCL, PPM, MP5, CWENO, THINC, WENO5). This ensures exact preservation of stationary "lake at rest" solutions over variable bathymetry, eliminating spurious velocity oscillations when properly configured.
- The explicit bed slope source term is discretized using the same cell-centered bathymetry as the hydrostatic reconstruction, guaranteeing consistency and improved balance.
- For stationary/well-balanced tests, it is recommended to use a high-order scheme (e.g., MUSCL+van Leer, PPM, MP5) and set MATLAB ODE solver tolerances (AbsTol, RelTol) to 1e-9 or tighter for machine-precision accuracy.
- The configuration file (`simulation_config.m`) allows easy switching between high-order schemes, limiters, and ODE tolerances for different experiments.
- Advanced reconstruction methods (MUSCL, PPM, MP5, CWENO, THINC) are now available and selectable in the config, supporting both component-wise and characteristic-based reconstruction. For MUSCL, use `cfg.reconstruct.characteristic`; for PPM, use `cfg.reconstruct.ppm_mode`; for MP5, use `cfg.reconstruct.mp5_mode` (default: 'characteristic', requires ng=3). The dam_break test case uses 5th-order characteristic MP5 with RK4.
- **Troubleshooting:** If you observe small spurious velocities in a stationary test, ensure you are using a well-balanced configuration and tight ODE tolerances as described above.

**Latest Well-Balanced High-Order Scheme:**
- The function [`+core/rhs_nsw_high_order.m`](./+core/rhs_nsw_high_order.m) implements a well-balanced high-order finite volume scheme for the 1D Nonlinear Shallow Water equations, based on hydrostatic reconstruction and a centered source term as described in Audusse et al. (2004) and related literature (see function header for details).
- Hydrostatic reconstruction at cell interfaces ensures exact preservation of stationary "lake at rest" states over arbitrary bathymetry.
- The explicit bed slope source term is discretized in a centered and consistent manner with the interface treatment, improving stability and well-balancing.
- Modular support for high-order reconstruction methods: MUSCL (with limiters), PPM (3rd order), MP5 (5th order), CWENO, THINC, WENO5, with both component-wise and characteristic-based options. Configure via `cfg.reconstruct` (see below).
- Robust handling of dry states, ghost cells, and boundary conditions.
- For stationary/well-balanced tests, use a high-order scheme (e.g., MUSCL+van Leer, PPM, MP5) and set MATLAB ODE solver tolerances (AbsTol/RelTol) to 1e-9 or tighter for machine-precision accuracy.
- The `lake_at_rest` initial condition is now truly well-balanced for arbitrary bathymetry.

**Performance Optimizations (NEW):**
- **MUSCL Reconstruction Acceleration:** The MUSCL reconstruction algorithm in `+reconstruct/muscl.m` has been optimized by vectorizing the inner loop, resulting in significant performance improvements for large-scale simulations.
- The van Albada limiter function now correctly receives the configuration (`cfg`) parameter, ensuring proper handling of epsilon values for near-zero gradients.

**New Test Cases (NEW):**
- **Sloping Beach Simulation:** A new bathymetry function `+bathy/sloping_beach.m` has been added that creates a flat bottom for the first 2/3 of the channel followed by a sloping beach. This enables simulation of wave runup and interaction with a beach.
- The visualization in `+vis/plot_state.m` has been enhanced to correctly display the free surface at y=0 and the bottom at y=-1 for the flat region, with proper positioning of boundary markers at the bottom of the tank.
- The sloping beach case can be configured in `+cfg/simulation_config.m` with the 'sloping_beach' experiment setup.
- See the function header in `rhs_nsw_high_order.m` for authorship and literature references.
- **Configuration:**
  - Select the reconstruction method and mode (component-wise or characteristic) via `cfg.reconstruct` fields:
    - For characteristic-based schemes, set `cfg.reconstruct.characteristic = true` (MUSCL, WENO5, THINC), `cfg.reconstruct.ppm_mode = 'characteristic'` (PPM), or `cfg.reconstruct.mp5_mode = 'characteristic'` (MP5).
  - Example configuration snippets are provided in the "Configuration Details" and "High-Order Reconstruction" sections below.
- **Troubleshooting:** If you observe small spurious velocities in a stationary test, ensure you are using a well-balanced configuration and tight ODE tolerances as described above. Consult the README and code documentation for further tips.

`1DWaveTank` is a MATLAB-based numerical laboratory for simulating 1D long wave phenomena using the Finite Volume method. It provides a modular framework for experimenting with different numerical schemes, boundary conditions, and initial states for the **Non-Linear Shallow Water (NSW) equations**. The design prioritizes modularity, readability, and ease of extension, allowing users to readily test and integrate new components, even if this means sacrificing raw computational speed compared to highly optimized, monolithic codes.

The NSW equations solved here typically conserve mass and momentum, represented by the state variables:
* `H`: Water depth [m]
* `HU`: Discharge (or momentum) [m²/s]

This document details the structure, components, and usage of the `1DWaveTank` code.

## 2. Project Structure

The code is organised into MATLAB packages for modularity:

- **`+cfg`**: Configuration files and functions.
- **`+core`**: Core solver logic (`solver.m`) and RHS functions (`rhs_nsw_*.m`).
- **`+flux`**: Numerical flux functions (e.g., `HLLC.m`, `Rusanov.m`).
- **`+reconstruct`**: High-order reconstruction methods (e.g., `muscl.m`, `weno5.m`, `ppm.m`, `mp5.m`), supporting both component-wise and characteristic-based approaches for MUSCL, PPM, and MP5. See below for configuration.
- **`+bc`**: Boundary condition implementations (e.g., `wall.m`, `open.m`, `generating.m`). The generating BC fills only the outermost ghost cell using Riemann invariants and uses constant extrapolation for additional ghost cells.
- **`+ic`**: Initial condition setups (e.g., `lake_at_rest.m` which sets a flat free surface over the given bathymetry).
- **`+time`**: Time integration schemes
- **`+vis`**: Visualisation tools
- **`+bathy`**: Bathymetry definition functions.
- **`+friction`**: Bottom friction models.
- **`run_simulation.m`**: Main script to execute a simulation
- **`README.md`**: Project overview and basic setup guide

## 3. Getting Started

1. **Configure:** Edit `+cfg/simulation_config.m` to select an `experiment_setup` (e.g., `'flat_rest'`, `'flat_gaussian'`, `'flat_wave_gen'`) and customise parameters like domain, mesh size (`N`), end time (`tEnd`), numerical flux (`numFlux`), time stepper (`timeStepper`), boundary conditions (`bc.left.handle`, `bc.right.handle`), initial condition (`ic_handle`), etc.
2. **Run:** Execute the main script:
   ```matlab
   run_simulation
   ```
3. **Observe:** The script displays simulation progress and animates the results. The final output is stored in the `results` structure.

## 4. Configuration (`+cfg`)

Configuration controls the simulation setup hierarchically: Defaults -> Experiment Setup -> User Overrides.

### 4.1. `default_config.m`

- Provides baseline parameters and default function handles.
- Called first by `simulation_config.m`.

### 4.2. `simulation_config.m`

- Main configuration file.
- Loads defaults and applies scenario-specific overrides.
- Sets mesh, time stepping, physical model, output path, etc.
- **Example Override:**
  ```matlab
  % Inside simulation_config.m, after loading defaults/setup:
  switch experiment_setup
      case 'flat_gaussian'
          % ... other settings ...
          config.tEnd = 10.0; % Override default end time
          config.numFlux = @flux.LaxFriedrichs; % Choose Lax-Friedrichs flux
          config.ic_param.a = 0.1; % Set IC parameter
          config.bathyHandle = @bathy.flat; % Set flat bathymetry
          config.icHandle = @ic.solitary_wave; % Set solitary wave IC
          % ...
  ```

### 4.3. `+bathy`

- Contains functions defining the bottom elevation `h(x)`.
- **`flat.m`**: Implements flat bathymetry:
  ```matlab
  function h = flat(x, cfg)
  ```

## 5. Core Components (`+core`)

### 5.1. `solver.m`

- Main simulation orchestrator.
- Sets IC, selects RHS, calls time integrator, computes `H`, `HU`, `U`.
- Returns results with metadata.

### 5.2. `rhs_nsw_1st_order.m`

- Computes RHS using 1st order FV scheme.
- Applies BCs, computes numerical fluxes, handles source terms.
- **Well-Balancing Note:** The bed slope source term is currently commented out. For simulations with non-flat bathymetry (`bathyHandle` is not `@bathy.flat`), implementing a **well-balanced** scheme (either by modifying the flux or the source term calculation) is crucial. This ensures that steady states (like a lake at rest over a varying bottom) are correctly maintained and prevents spurious flows generated by numerical imbalance between flux gradients and source terms.

### 5.3. `rhs_nsw_high_order.m`

- Computes RHS using a 2nd order FV scheme with MUSCL reconstruction.
- Requires a reconstruction method to be specified in `cfg.reconstruct`.
- Calls the selected reconstruction function (e.g., `@reconstruct.muscl` or `@reconstruct.muscl_characteristic`) to compute left and right states at cell interfaces.
- Applies BCs, computes numerical fluxes using reconstructed states, handles source terms.
- Suitable for use with higher-order time integrators (e.g., SSP2, SSP3, RK4).

### 5.4. Core Utilities (`+core/+utils`)

- `calculate_dt_cfl.m`, `get_bc_style.m`, `get_param.m`, `odetpbar.m`, `physical_flux.m`, `struct2str.m`, `textprogressbar.m`, `uniform.m`.

## 6. Numerical Fluxes (`+flux`)

Implemented schemes:

- `AUSM+.m`  (AUSM+ variant)
- `AUSMDV.m` (AUSM-Derivative Variant)
- `FORCE.m`
- `FVCF.m`
- `HLL.m`
- `HLLC.m`
- `LF.m`
- `LaxFriedrichs.m`
- `OsherSolomon.m`
- `Roe.m`
- `Rusanov.m`
- `StegerWarming.m`
- `HLLE.m`
- `SLAU.m`
- `CentralUpwind.m`: Central-Upwind (Kurganov-Noelle-Petrova) flux.
- `PVM.m`: Pressure Velocity Momentum flux (specific details TBD).
- `Kinetic.m`: Kinetic flux (contributed by Prof. Mehmet ERSOY, Université de Toulon, IMATH, based on original Fortran code).

All take `(wL, wR, cfg)` and return `[Flux_H, Flux_HU]`.

To use a specific flux, set the `config.numFlux` function handle in `+cfg/simulation_config.m`, for example:

```matlab
config.numFlux = @flux.HLLC; % Use the HLLC flux
```

## 7. High-Order Reconstruction (`+reconstruct`)

To achieve second-order spatial accuracy, the finite volume method requires reconstructing the solution within each cell to obtain distinct values at the left and right sides of cell interfaces. This package provides implementations for the MUSCL (Monotonic Upstream-centered Scheme for Conservation Laws) approach.

### Reconstruction Schemes (`+reconstruct`)

This package provides functions for higher-order spatial reconstruction of cell-averaged values to cell interfaces. This is essential for reducing numerical diffusion and capturing sharper gradients.

- **MUSCL** (`muscl.m`): Supports both component-wise and characteristic-based reconstruction (set `cfg.reconstruct.characteristic`).
- **PPM** (`ppm.m`): Supports both component-wise and characteristic-based reconstruction (set `cfg.reconstruct.ppm_mode`).
- **MP5** (`mp5.m`): Implements the 5th-order Monotonicity Preserving method (Suresh & Huynh, 1997), supports both component-wise and characteristic-based reconstruction (set `cfg.reconstruct.mp5_mode`, default 'characteristic'), requires `ng=3` ghost cells. Used in the dam_break test case with RK4.


- **`weno5.m`**: Implements the 5th-order Weighted Essentially Non-Oscillatory (WENO5) scheme. It provides high accuracy while controlling spurious oscillations near discontinuities. It can reconstruct variables component-wise or use a characteristic decomposition (activated by `cfg.reconstruct.characteristic = true`) for improved stability, especially for systems like the shallow water equations.
- **`muscl.m`**: Implements the 2nd-order Monotonic Upstream-centered Scheme for Conservation Laws (MUSCL). This is a widely used, robust scheme. It requires a slope limiter function (specified via `cfg.reconstruct.limiter`, e.g., `@reconstruct.limiters.minmod`) to prevent oscillations. Similar to WENO5, it can operate component-wise or use a characteristic decomposition (`cfg.reconstruct.characteristic = true`).
- **`ppm.m`**: Implements the 3rd-order Piecewise Parabolic Method (PPM). It offers high accuracy in smooth regions and robust shock-capturing. Supports both component-wise and characteristic-based reconstruction (controlled by `cfg.reconstruct.ppm_mode`). Requires `ng >= 2` ghost cells.
- **`mp5.m`**: Implements the 5th-order Monotonicity Preserving (MP5) scheme. Based on Suresh & Huynh (1997), it achieves 5th-order accuracy in smooth regions while strictly preserving monotonicity near discontinuities using specialized limiters. Supports both component-wise and characteristic-based reconstruction (controlled by `cfg.reconstruct.mp5_mode`, default is 'characteristic'). Requires `ng >= 3` ghost cells.
- **`cweno.m`**: Implements the Central Weighted Essentially Non-Oscillatory (CWENO) method. Uses a central optimal stencil and two sub-stencils for each interface, achieving high-order accuracy with fewer points than traditional WENO. It is especially accurate at smooth extrema and robust near discontinuities. Supports both component-wise and characteristic-based reconstruction (controlled by `cfg.reconstruct.cweno_mode`), with characteristic mode using Roe-averaged eigenvectors. Falls back to first-order at boundaries and handles dry states robustly.
- **`thinc.m`**: Implements the THINC (Tangent of Hyperbola for Interface Capturing) method. Supports both component-wise and characteristic-based reconstruction (controlled by `cfg.reconstruct.characteristic`). The characteristic mode uses Roe-averaged eigenvectors for improved shock resolution. Robust boundary handling is ensured by clamping stencil indices within valid bounds.
- **`cweno.m`**: Implements the Central Weighted Essentially Non-Oscillatory (CWENO) method (Levy, Puppo & Russo, 1999). Supports both component-wise and characteristic-based reconstruction (controlled by `cfg.reconstruct.cweno_mode`). CWENO uses a central optimal stencil and two sub-stencils for each interface, achieving high-order accuracy and improved performance at smooth extrema with fewer points than traditional WENO. Characteristic mode leverages Roe-averaged eigenvectors for robust shock capturing. Near boundaries, the method falls back to first-order. Water depths are enforced non-negative, and dry states are handled robustly.
- **`+limiters`**: A sub-package containing various slope limiter functions (`minmod.m`, `vanleer.m`, `mc.m`, `superbee.m`) for use with the MUSCL scheme. The choice of limiter affects the scheme's accuracy and dissipative properties.

### 7.2. Slope Limiters (`+reconstruct/+limiters`)

Slope limiters modify the calculated slopes within each cell to prevent oscillations while maintaining accuracy. Several limiters are available:

- `minmod.m`: Basic, robust, but relatively dissipative.
- `mc.m`: Monotonized Central limiter.
- `koren.m`: Third-order accurate limiter.
- `ospre.m`: A smoother limiter.
- `umist.m`: UMIST limiter.
- `sweby.m`: Sweby limiter (parameterizable between minmod and superbee).
- `superbee.m`: Superbee limiter (less dissipative than minmod).
- `vanleer.m`: Van Leer's limiter.
- `vanalbada.m`: Van Albada limiter.

### 7.3. Configuration

To use high-order reconstruction:

1.  Set the RHS function to use the high-order solver:
    ```matlab
    config.model = @core.rhs_nsw_high_order;
    ```
2.  Specify the reconstruction method and limiter in the `cfg.reconstruct` structure:
    ```matlab
    % Option 1: Component-wise MUSCL with OSPRE limiter
    config.reconstruct.method = 'muscl';
    config.reconstruct.handle = @reconstruct.muscl;
    config.reconstruct.order = 2;
    config.reconstruct.limiter = @reconstruct.limiters.ospre;

    % Option 2: Component-wise ENO2
    config.reconstruct.method = 'eno2';
    config.reconstruct.handle = @reconstruct.eno2;
    config.reconstruct.order = 2;
    config.model = @core.rhs_nsw_high_order;
    config.timeStepper = @time.integrate_ssp2_adaptive; % Need >= 2nd order time stepper

    % Option 3: Characteristic MUSCL with Koren limiter
    config.reconstruct.method = 'muscl_characteristic';
    config.reconstruct.handle = @reconstruct.muscl_characteristic;
    config.reconstruct.order = 2;
    config.reconstruct.limiter = @reconstruct.limiters.koren;

    % Option 4: UNO2 Reconstruction
    config.reconstruct.method = 'uno2';
    config.reconstruct.handle = @reconstruct.uno2;
    config.reconstruct.order = 2;
    % Limiter is not specified for UNO2 as it's built-in

    % Option 5: WENO5 Reconstruction
    config.reconstruct.method = 'weno5';
    config.reconstruct.handle = @reconstruct.weno5;
    config.reconstruct.order = 5;
    % Limiter is not specified for WENO5 as it's built-in

    % Option 6: MP5 Reconstruction
    config.reconstruct.method = 'mp5';
    config.reconstruct.handle = @reconstruct.mp5;
    config.reconstruct.order = 5;
    config.reconstruct.mp5_mode = 'characteristic';
    
    % Option 7: THINC Reconstruction (component-wise or characteristic)
    config.reconstruct.method = 'thinc';
    config.reconstruct.handle = @reconstruct.thinc;
    config.reconstruct.thinc_beta = 1.5; % Steepness parameter (default 1.5, typical 1.5-2.5)
    config.reconstruct.characteristic = true; % Set to true for characteristic-based THINC
    % Note: Stencil indices are automatically clamped for robust boundary handling.

    % Option 8: CWENO Reconstruction (component-wise or characteristic)
    config.reconstruct.method = 'cweno';
    config.reconstruct.handle = @reconstruct.cweno;
    config.reconstruct.cweno_mode = 'characteristic'; % or 'component'
    % Requires at least 2 ghost cells (ng >= 2)
    ```

## 8. Boundary Conditions (`+bc`)

Boundary conditions define how the simulation behaves at the edges of the domain (`xmin` and `xmax`). They are set using function handles in `cfg.bc.left.handle` and `cfg.bc.right.handle`.

*   **`@bc.wall`:** Simulates a solid, impermeable wall. This is the default condition. It sets the velocity (HU) to zero in the ghost cells and reflects the water height (H).
*   **`@bc.generating`:** Simulates a wave generation boundary. Currently configured for sinusoidal waves based on amplitude (`cfg.bc.left.param.a` or `cfg.bc.right.param.a`) and frequency (`cfg.bc.left.param.omega` or `cfg.bc.right.param.omega`).
*   **`@bc.periodic`:** Connects the left and right boundaries, treating the domain as if it wraps around. Values leaving one side enter the other. To use periodic conditions, *both* `cfg.bc.left.handle` and `cfg.bc.right.handle` must be set to `@bc.periodic`.

Example configuration for periodic BCs:
```matlab
cfg.bc.left.handle = @bc.periodic;
cfg.bc.right.handle = @bc.periodic;
```

## 9. Initial Conditions (`+ic`)

Defines the initial state (`H`, `HU`) at `t=0`.

- `gaussian_bump.m`: Gaussian perturbation on a lake at rest.
- `lake_at_rest.m`: Sets a flat free surface at elevation `z = cfg.h0` over the specified bathymetry (`cfg.bathyHandle`), with zero initial velocity. The initial depth `H` is calculated as `H(x, 0) = cfg.h0 - b(x)`, respecting `cfg.phys.dry_tolerance`. Requires the full `cfg` structure.
- `solitary_wave.m`: Solitary wave initial condition.
- `dam_break.m`: Piecewise constant initial water level. Uses `cfg.dam_break.h_L`, `cfg.dam_break.h_R`, `cfg.dam_break.x_dam`.

Signature:
```matlab
function w0 = ic_function(cfg)
```

## 10. Bathymetry (`+bathy`)

Bathymetry functions define the bottom elevation `h(x)`.

- `flat.m`: Flat bathymetry (constant depth `cfg.param.H0`).

Signature:
```matlab
function h = bathy_function(x, cfg)
```

## 11. Friction Models (`+friction`)

The 1DWaveTank code supports various friction models through the `+friction` package. By default, no friction is applied.

### Available Models

* **No Friction** (default): `config.phys.friction_model = @friction.no_friction`
* **Chézy**: `config.phys.friction_model = friction.friction_selector('chezy')`
  * Requires: `config.phys.chezy_C` (typical values: 30-90 m^(1/2)/s)
* **Manning**: `config.phys.friction_model = friction.friction_selector('manning')`
  * Requires: `config.phys.manning_n` (typical values: 0.01-0.05 s/m^(1/3))
  * Common values:
    * 0.01-0.02: Very smooth channels (concrete, metal)
    * 0.025-0.035: Clean, straight natural channels
    * 0.04-0.07: Natural channels with vegetation, stones
    * 0.08-0.15: Very rough channels with heavy vegetation
* **Darcy-Weisbach**: `config.phys.friction_model = friction.friction_selector('darcy_weisbach')`
  * With constant friction factor:
    * Set `config.phys.darcy_f` (typical values: 0.01-0.05)
  * With Colebrook-White formula (variable friction factor):
    * Set `config.phys.f_calculation = 'colebrook_white'`
    * Set `config.phys.ks` - Equivalent sand roughness [m]
    * Set `config.phys.kinematic_viscosity` - Kinematic viscosity [m²/s]
    * Optional: `config.phys.cw_iterations`, `config.phys.cw_tolerance`

### Using Friction in Simulations

To enable a friction model, add these lines to your configuration in `simulation_config.m`:

```matlab
% Option 1: Chézy model
config.phys.friction_model = friction.friction_selector('chezy');
config.phys.chezy_C = 50;  % Chézy coefficient [m^(1/2)/s]

% Option 2: Manning model
config.phys.friction_model = friction.friction_selector('manning');
config.phys.manning_n = 0.03;  % Manning coefficient [s/m^(1/3)]

% Option 3: Darcy-Weisbach with constant friction factor
config.phys.friction_model = friction.friction_selector('darcy_weisbach');
config.phys.darcy_f = 0.02;  % Constant Darcy friction factor

% Option 4: Darcy-Weisbach with Colebrook-White formula
config.phys.friction_model = friction.friction_selector('darcy_weisbach');
config.phys.f_calculation = 'colebrook_white';
config.phys.ks = 0.001;  % Equivalent sand roughness [m]
config.phys.kinematic_viscosity = 1e-6;  % Kinematic viscosity [m²/s]
```

### Adding New Friction Models

To add a new friction model:

1. Create a new file in `+friction/` (e.g., `manning.m`)
2. Implement the friction model with the standard interface: `friction_term = model_name(H, HU, g, cfg)`
3. Add the model to `friction_selector.m`
4. Update documentation

## 12. Time Integration (`+time`)

### Custom Adaptive Steppers

- `integrate_ab2_adaptive.m`: 2nd-order Adams-Bashforth.
- `integrate_bogacki_shampine.m`: 3(2) embedded Runge-Kutta (Bogacki-Shampine).
- `integrate_dopri54_adaptive.m`: 5(4) embedded Runge-Kutta (Dormand-Prince).
- `integrate_euler_adaptive.m`: Forward Euler.
- `integrate_rk4_adaptive.m`: 4th-order Runge-Kutta.
- `integrate_ssp2_adaptive.m`: 2nd-order SSP Runge-Kutta.
- `integrate_ssp3_adaptive.m`: 3rd-order SSP Runge-Kutta.

All follow signature:
```matlab
function [sol_out, t_out, stats] = integrate_...(rhs_func, tspan, w0, cfg)
```

All main simulation scripts, solver, time integration, configuration, CFL utility, and visualization files now include comprehensive comments and improved documentation for clarity.
### MATLAB Solver Wrapper

- `integrate_matlab_ode.m`
- Uses `cfg.time.matlab_solver`, `AbsTol`, `RelTol`, `ode_options`
- Uses internal error control for time stepping (not CFL)

## 13. Visualisation (`+vis`)

### `plot_state.m`

- Plots water surface and bathymetry.
- Optionally plots velocity.
- Handles figure updating and axis management.
- Used by `run_simulation.m`