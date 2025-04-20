# 1DWaveTank: User Guide

**Author:** Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)  
**Date:** April 20, 2025

## 1. Overview

`1DWaveTank` is a MATLAB-based numerical laboratory for simulating 1D long wave phenomena using the Finite Volume method. It provides a modular framework for experimenting with different numerical schemes, boundary conditions, and initial states for the Non-Linear Shallow Water (NSW) equations.

This document details the structure, components, and usage of the `1DWaveTank` code.

## 2. Project Structure

The code is organised into MATLAB packages for modularity:

- **`+cfg`**: Configuration files and functions
- **`+core`**: Core solver logic and RHS (Right-Hand Side) functions
- **`+flux`**: Numerical flux implementations for the finite volume method
- **`+bc`**: Boundary condition implementations
- **`+ic`**: Initial condition setups
- **`+time`**: Time integration schemes
- **`+vis`**: Visualisation tools
- **`run_simulation.m`**: Main script to execute a simulation
- **`README.md`**: Project overview and basic setup guide

## 3. Getting Started

1. **Configure:** Edit `+cfg/simulation_config.m` to select an `experiment_setup` (e.g., `'flat_rest'`, `'flat_gaussian'`, `'flat_wave_gen'`) and customise parameters.
2. **Run:** Execute the main script:
   ```matlab
   run_simulation
   ```
3. **Observe:** The script displays simulation progress and animates the results. The final output is stored in the `results` structure.

## 4. Configuration (`+cfg`)

### 4.1. `default_config.m`

- Provides baseline parameters and default function handles.
- Called first by `simulation_config.m`.

### 4.2. `simulation_config.m`

- Main configuration file.
- Loads defaults and applies scenario-specific overrides.
- Sets mesh, time stepping, physical model, output path, etc.

### 4.3. `+cfg/+bathy`

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

### 5.3. Core Utilities (`+core/+utils`)

- `calculate_dt_cfl.m`, `get_bc_style.m`, `get_param.m`, `odetpbar.m`, `physical_flux.m`, `struct2str.m`, `textprogressbar.m`, `uniform.m`.

## 6. Numerical Fluxes (`+flux`)

Implemented schemes:

- `FORCE.m`
- `FVCF.m`
- `HLL.m`
- `HLLC.m`
- `LF.m`
- `OsherSolomon.m`
- `Roe.m`
- `Rusanov.m`
- `StegerWarming.m`

All take `(wL, wR, cfg)` and return `[Flux_H, Flux_HU]`.

## 7. Boundary Conditions (`+bc`)

Boundary condition functions:

- `generating.m`
- `open.m`
- `placeholder_bc_1st_order.m`
- `wall.m`

Signature:
```matlab
function w_padded = bc_function(w_padded, t, side, cfg, num_ghost_cells)
```

## 8. Initial Conditions (`+ic`)

Initial condition functions:

- `gaussian_bump.m`: Gaussian water surface bump, zero velocity.
- `lake_at_rest.m`: Constant water depth, zero velocity.

Signature:
```matlab
function w0 = ic_function(xc, cfg)
```

## 9. Time Integration (`+time`)

### Custom Adaptive Steppers

- `integrate_ab2_adaptive.m`
- `integrate_euler_adaptive.m`
- `integrate_rk4_adaptive.m`
- `integrate_ssp2_adaptive.m`
- `integrate_ssp3_adaptive.m`

All follow signature:
```matlab
function [sol_out, t_out, stats] = integrate_...(rhs_func, tspan, w0, cfg)
```

### MATLAB Solver Wrapper

- `integrate_matlab_ode.m`
- Uses `cfg.time.matlab_solver`, `AbsTol`, `RelTol`, `ode_options`

## 10. Visualisation (`+vis`)

### `plot_state.m`

- Plots water surface and bathymetry.
- Optionally plots velocity.
- Handles figure updating and axis management.
- Used by `run_simulation.m` for animation.

## 11. Main Script (`run_simulation.m`)

- Clears environment, adds paths.
- Loads config and runs the simulation.
- Calls visualisation tools.
- Prints CPU time and global statistics.