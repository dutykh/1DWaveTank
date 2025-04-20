# 1DWaveTank: A Finite Volume Numerical Wave Tank

## Overview

`1DWaveTank` is a MATLAB-based numerical laboratory designed for simulating and analyzing long wave phenomena in one dimension. It serves as a flexible framework for implementing, testing, and comparing various mathematical models (both dispersive and non-dispersive) and numerical schemes, primarily focusing on the finite volume method.

The core philosophy is to provide a modular and extensible platform where researchers and students can easily experiment with different physical setups, numerical algorithms, and boundary conditions related to shallow water wave dynamics.

## Author

This project was initiated and is maintained by:

*   **Dr. Denys Dutykh**
    *   *Khalifa University of Science and Technology, Abu Dhabi, UAE*

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
- Modular numerical flux functions (Flux Vector Central Flux (FVCF) implemented)
- Adaptive Forward Euler time stepping based on a CFL condition
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

This project is released under the MIT License. See the LICENSE file for details.

- Main author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
- Please cite appropriately if you use this code for research or teaching.
