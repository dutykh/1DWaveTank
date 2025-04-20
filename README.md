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

*   **`+cfg`**: Configuration files (`simulation_config.m`, `default_config.m`). Defines simulation parameters, physical setup, numerical choices, and run control.
*   **`+core`**: Core solver components (`solver.m`, `rhs_*.m`, utils). Contains the main time-stepping logic and the functions defining the right-hand side (RHS) of the governing equations.
*   **`+flux`**: Numerical flux functions (e.g., `FVCF.m`). Implements different finite volume flux calculators.
*   **`+bc`**: Boundary condition implementations (e.g., `wall.m`, `generating.m`). Defines how the boundaries of the computational domain are handled.
*   **`+ic`**: Initial condition setups (e.g., `lake_at_rest.m`). Defines the initial state of the system (water elevation, velocity).
*   **`+time`**: Time integration schemes (e.g., `integrate_euler_adaptive.m`). Contains different methods for advancing the solution in time.
*   **`+vis`**: Visualization tools (`plot_state.m`). Functions for plotting the simulation results.
*   **`+test`**: (Optional/Future) Unit tests and validation cases.
*   **`run_simulation.m`**: The main script to configure, run, and visualize a simulation.

## Features

The current implementation provides:

- Non-linear Shallow Water (NSW) equations solver
- 1st Order Finite Volume Method
- Modular numerical flux functions (FVCF, HLL, HLLC, LF, Rusanov, Roe implemented)
- Adaptive Forward Euler time stepping based on a CFL condition
- Adaptive Strong Stability Preserving Runge-Kutta (SSP) schemes:
  - SSP(2,2)
  - SSP(3,3)
- Wrapper for standard MATLAB ODE solvers (e.g., `ode45`, `ode113`, `ode23`) with:
  - Configurable solver choice (`cfg.time.matlab_solver`)
  - Optional custom `odeset` options (`cfg.time.ode_options`)
  - Optional text progress bar (`cfg.time.show_progress_bar`)
- Configurable domain, mesh, and simulation parameters
- Solid wall and wave-generating boundary conditions
- Modular initial condition interface (lake at rest, Gaussian bump, etc.)
- Modular, extensible configuration system
- Visualization tools for water surface, velocity, and bathymetry
- Output control: results saving can be toggled via configuration
- MATLAB package-based project structure for clarity and extensibility
- Example experiment setups for rapid testing
- (Optional) Unit test and validation framework structure

## Getting Started

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/dutykh/1DWaveTank 1DWaveTank
    cd 1DWaveTank
    ```
2.  **Configure Simulation:**
    *   Open `+cfg/simulation_config.m`.
    *   Select the desired `experiment_setup` (e.g., `'flat_wave_gen'`).
    *   Adjust parameters within the selected setup or the common settings (domain size, mesh resolution, wave parameters, end time, CFL number, etc.).
3.  **Run Simulation:**
    *   Execute the main script from the MATLAB command window:
        ```matlab
        run_simulation
        ```
4.  **Observe Results:** The script will output simulation progress to the console and display an animation of the wave tank evolution.

## Contributing

Contributions are highly welcome! This project aims to be a collaborative environment for exploring finite volume methods for wave modeling.

- If you add a new feature or fix a bug, please document your changes and, if possible, provide a validation or unit test.
- New tests and beta testers are especially welcome—if you find issues or have suggestions, please open an issue or submit a pull request.
- Please follow the existing MATLAB package structure for any new code.

## License

This project is licensed under the **GNU General Public License v3.0**.

- Copyright (C) 2024 Dr. Denys Dutykh
- See the `LICENSE` file for the full license text.

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

- Main author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
- Please cite appropriately if you use this code for research or teaching.
