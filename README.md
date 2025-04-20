# 1DWaveTank: A Finite Volume Numerical Wave Tank

## Overview

`1DWaveTank` is a MATLAB-based numerical laboratory designed for simulating and analyzing long wave phenomena in one dimension. It serves as a flexible framework for implementing, testing, and comparing various mathematical models (both dispersive and non-dispersive) and numerical schemes, primarily focusing on the finite volume method.

The core philosophy is to provide a modular and extensible platform where researchers and students can easily experiment with different physical setups, numerical algorithms, and boundary conditions related to shallow water wave dynamics.

## Author

This project was initiated and is maintained by:

*   **Dr. Denys Dutykh**
    *   *Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE*

## Project Structure

The codebase is organized using MATLAB packages (directories starting with `+`) to promote modularity and clarity:

*   **`+cfg`**: Configuration files ([`simulation_config.m`](./+cfg/simulation_config.m), [`default_config.m`](./+cfg/default_config.m)). Defines simulation parameters, physical setup, numerical choices, and run control.
*   **`+core`**: Core solver components ([`solver.m`](./+core/solver.m), [`rhs_*.m`](./+core/), utils). Contains the main time-stepping logic and the functions defining the right-hand side (RHS) of the governing equations.
*   **[`+flux`](./+flux/)**: Numerical flux functions (e.g., `FVCF.m`, `OsherSolomon.m`, `StegerWarming.m`, `FORCE.m`). Implements different finite volume flux calculators.
*   **[`+bc`](./+bc/)**: Boundary condition implementations (e.g., `wall.m`, `generating.m`). Defines how the boundaries of the computational domain are handled.
*   **[`+ic`](./+ic/)**: Initial condition setups (e.g., `lake_at_rest.m`). Defines the initial state of the system (water elevation, velocity).
*   **[`+time`](./+time/)**: Time integration schemes (e.g., `integrate_euler_adaptive.m`). Contains different methods for advancing the solution in time.
*   **[`+vis`](./+vis/)**: Visualization tools (`plot_state.m`). Functions for plotting the simulation results.
*   **`+test`**: (Optional/Future) Unit tests and validation cases.
*   **[`run_simulation.m`](./run_simulation.m)**: The main script to configure, run, and visualize a simulation.

## Features

*   Non-linear Shallow Water (NSW) equations solver
*   1st Order Finite Volume Method framework
*   **Numerical Fluxes:** Modular functions available in [`+flux/`](./+flux/) including:
    *   FVCF, HLL, HLLC, Rusanov, Roe, Osher-Solomon, Steger-Warming, FORCE
*   **Time Integration:** Adaptive time stepping based on a CFL condition available in [`+time/`](./+time/) for:
    *   Forward Euler (`integrate_euler_adaptive.m`)
    *   SSP(2,2) (`integrate_ssp2_adaptive.m`)
    *   SSP(3,3) (`integrate_ssp3_adaptive.m`)
    *   Explicit RK4 (`integrate_rk4_adaptive.m`)
    *   Adams-Bashforth 2nd Order (AB2) (`integrate_ab2_adaptive.m`)
    *   Standardized output format (`[sol_out, t_out, stats]`) for custom adaptive steppers.
*   **MATLAB ODE Solvers:** Wrapper (`integrate_matlab_ode.m`) for standard solvers (e.g., `ode45`, `ode113`, `ode23`) with:
    *   Configurable solver choice (`cfg.time.matlab_solver`)
    *   Optional custom `odeset` options (`cfg.time.ode_options`)
*   **Boundary Conditions:** Implementations in [`+bc/`](./+bc/) including:
    *   Solid Wall (`wall.m`)
    *   Wave Generating (`generating.m`)
*   **Initial Conditions:** Setups in [`+ic/`](./+ic/) including:
    *   Lake at Rest (`lake_at_rest.m`)
    *   Gaussian Bump (`gaussian_bump.m`)
*   Configurable domain, mesh, and simulation parameters
*   Modular, extensible configuration system
*   Visualization tools for water surface, velocity, and bathymetry
*   Output control: results saving can be toggled via configuration
*   MATLAB package-based project structure for clarity and extensibility
*   Example experiment setups for rapid testing

## Getting Started

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/dutykh/1DWaveTank 1DWaveTank
    cd 1DWaveTank
    ```
2.  **Configure Simulation:**
    *   Open [`+cfg/simulation_config.m`](./+cfg/simulation_config.m).
    *   **Select Experiment Setup:** Choose a pre-defined setup by uncommenting the corresponding `cfg = cfg.experiment_setups.[setup_name](cfg);` line (e.g., `flat_wave_gen`). These setups define specific combinations of initial conditions, boundary conditions, and physical parameters.
    *   **Customize Parameters:** Modify parameters directly within `simulation_config.m` *after* loading the default and experiment setups. Key areas include:
        *   **Domain & Mesh:** `cfg.domain.xmin`, `cfg.domain.xmax`, `cfg.mesh.N`
        *   **Time:** `cfg.time.T`, `cfg.time.CFL`, `cfg.vis.dt_plot`
        *   **Physics:** `cfg.phys.g`, `cfg.phys.Cf` (friction)
        *   **Numerics:**
            *   `cfg.time.integrator`: Select the time integration method (e.g., `@time.integrate_rk4_adaptive`, `@time.integrate_matlab_ode`).
            *   `cfg.time.matlab_solver`: If using `@time.integrate_matlab_ode`, specify the solver (e.g., `'ode45'`).
            *   `cfg.numFlux`: Choose the numerical flux function (e.g., `@flux.FORCE`, `@flux.OsherSolomon`).
        *   **Initial/Boundary Conditions:** Handles are typically set within the experiment setups, but can be overridden (e.g., `cfg.icHandle`, `cfg.bc.left.handle`, `cfg.bc.right.handle`). Specific parameters for BCs/ICs are often nested (e.g., `cfg.bc.left.param.a` for wave amplitude).
    *   See the **Configuration Details** section below for more information.
3.  **Run Simulation:**
    *   Execute the main script from the MATLAB command window:
        ```matlab
        run_simulation
        ```
4.  **Observe Results:** The script will output simulation progress to the console and display an animation of the wave tank evolution (if `cfg.vis.do_vis = true`). Solution data (`sol_out`, `t_out`) and statistics (`stats`) are returned by `core.solver`.

## Configuration Details (`+cfg/simulation_config.m`)

Configuration is handled hierarchically:

1.  **Defaults:** [`+cfg/default_config.m`](./+cfg/default_config.m) provides baseline values for all parameters.
2.  **Experiment Setup:** A function within `+cfg/experiment_setups/` (e.g., [`flat_wave_gen.m`](./+cfg/experiment_setups/flat_wave_gen.m)) is called from `simulation_config.m`. This function modifies the default `cfg` structure to define a specific scenario (IC, BCs, specific parameters).
3.  **User Overrides:** Direct assignments within `simulation_config.m` *after* calling the experiment setup function allow fine-tuning or overriding specific parameters for the current run.

Key `cfg` fields to customize:

*   `cfg.mesh`: Domain (`xmin`, `xmax`) and discretization (`N`, `dx`).
*   `cfg.phys`: Physical constants (`g`, `Cf`).
*   `cfg.param`: Model-specific parameters (e.g., `H0` for still water depth).
*   `cfg.time`: Time integration settings:
    *   `T`: Final simulation time.
    *   `CFL`: Courant-Friedrichs-Lewy number for adaptive steppers.
    *   `integrator`: Function handle (`@time...`) for the chosen time stepper.
    *   `matlab_solver`: String name (e.g., `'ode45'`) if using the MATLAB ODE wrapper.
    *   `ode_options`: Optional `odeset` structure for MATLAB solvers.
*   `cfg.numFlux`: Function handle (`@flux...`) for the numerical flux calculation.
*   `cfg.icHandle`: Function handle (`@ic...`) for the initial condition.
*   `cfg.bc`: Structure containing boundary condition settings:
    *   `cfg.bc.left.handle`, `cfg.bc.right.handle`: Function handles (`@bc...`).
    *   `cfg.bc.left.param`, `cfg.bc.right.param`: Structures holding parameters specific to the chosen BC functions (e.g., amplitude `a` and period `T` for `@bc.generating`).
*   `cfg.vis`: Visualization settings:
    *   `do_vis`: Boolean to enable/disable real-time plotting.
    *   `dt_plot`: Time interval between plot updates.
*   `cfg.output`: Output settings (e.g., saving results).

## Extending the Solver

The package structure makes adding new components straightforward:

1.  **New Numerical Flux (`+flux`)**:
    *   Create a new `.m` file in the [`+flux/`](./+flux/) directory (e.g., `my_new_flux.m`).
    *   Implement your flux function with the signature: `F = my_new_flux(wL, wR, cfg)`
        *   `wL`, `wR`: State vectors [H; HU] to the left and right of the interface.
        *   `cfg`: Configuration structure (can be used for parameters like `g`).
        *   `F`: Returned numerical flux vector [Flux_H; Flux_HU].
    *   Select your flux in `simulation_config.m`: `cfg.numFlux = @flux.my_new_flux;`

2.  **New Time Integrator (`+time`)**:
    *   Create a new `.m` file in the [`+time/`](./+time/) directory (e.g., `integrate_my_scheme.m`).
    *   Implement your integrator. Aim for the standard output signature for consistency: `[sol_out, t_out, stats] = integrate_my_scheme(rhs_func, tspan, w0, cfg)`
        *   `rhs_func`: Function handle to the RHS evaluation function (provided by `core.solver`).
        *   `tspan`: Vector of requested output times.
        *   `w0`: Initial condition vector (flattened).
        *   `cfg`: Configuration structure.
        *   `sol_out`: Solution matrix (rows are time points, columns are state variables).
        *   `t_out`: Vector of actual output times.
        *   `stats`: Structure with solver statistics (e.g., `nsteps`, `nfevals`).
    *   Select your integrator in `simulation_config.m`: `cfg.time.integrator = @time.integrate_my_scheme;`

3.  **New Boundary Condition (`+bc`)**:
    *   Create a new `.m` file in the [`+bc/`](./+bc/) directory (e.g., `my_boundary.m`).
    *   Implement your BC function with the signature: `w_padded = my_boundary(w_padded, t, side, cfg, num_ghost_cells)`
        *   `w_padded`: State array including ghost cells.
        *   `t`: Current time.
        *   `side`: String ('left' or 'right').
        *   `cfg`: Configuration structure (can hold BC-specific parameters in `cfg.bc.(side).param`).
        *   `num_ghost_cells`: Number of ghost cells to fill.
        *   The function should modify the appropriate ghost cell rows in `w_padded` and return the modified array.
    *   Select your BC in `simulation_config.m` (or an experiment setup): `cfg.bc.left.handle = @bc.my_boundary;` (and potentially set parameters in `cfg.bc.left.param`).

4.  **New Initial Condition (`+ic`)**:
    *   Create a new `.m` file in the [`+ic/`](./+ic/) directory (e.g., `my_initial_state.m`).
    *   Implement your IC function with the signature: `w0 = my_initial_state(cfg)`
        *   `cfg`: Configuration structure (contains mesh info like `cfg.mesh.xc`).
        *   `w0`: Returned initial state vector [H; HU] (N x 2 array or flattened 2N x 1).
    *   Select your IC in `simulation_config.m` (or an experiment setup): `cfg.icHandle = @ic.my_initial_state;`

## Contributing

Contributions are highly welcome! This project aims to be a collaborative environment for exploring finite volume methods for wave modeling.

*   If you add a new feature or fix a bug, please document your changes and, if possible, provide a validation or unit test.
*   New tests and beta testers are especially welcome—if you find issues or have suggestions, please open an issue or submit a pull request.
*   Please follow the existing MATLAB package structure for any new code.

## License

This project is licensed under the **GNU General Public License v3.0**.

*   Copyright (C) 2025 Dr. Denys Dutykh
*   See the `LICENSE` file for the full license text.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.

*   Main author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
*   Please cite appropriately if you use this code for research or teaching.