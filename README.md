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

## Current Implementation (Example)

The repository currently includes an example setup featuring:

*   **Model**: Non-linear Shallow Water (NSW) equations.
*   **Scheme**: 1st Order Finite Volume Method.
*   **Numerical Flux**: Flux Vector Central Flux (FVCF) scheme.
*   **Time Stepping**: Adaptive Forward Euler method based on a CFL condition.
*   **Boundary Conditions**: Solid wall, wave generation.
*   **Initial Conditions**: Lake at rest.

## Getting Started

1.  **Clone the Repository:**
    ```bash
    git clone <repository_url> 1DWaveTank
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

Contributions are highly welcome! This project aims to be a collaborative environment for exploring finite volume methods for wave modeling. If you have implemented or are interested in implementing:

*   Different numerical flux functions (e.g., Roe, HLL, HLLC, FORCE)
*   Higher-order finite volume reconstructions (e.g., MUSCL-Hancock, WENO)
*   Alternative time-stepping schemes (e.g., Runge-Kutta methods)
*   Different shallow water models (e.g., Boussinesq-type, Serre-Green-Naghdi)
*   New boundary or initial conditions

Please feel free to fork the repository, add your contributions following the existing package structure, and submit a pull request. We encourage clear documentation and, ideally, validation test cases for new numerical methods.

## Future Directions

*   Implementation of higher-order spatial and temporal schemes.
*   Inclusion of dispersive Boussinesq-type models.
*   Addition of more sophisticated boundary conditions (e.g., absorbing, moving).
*   Development of a comprehensive test suite for validation.
*   Optimization for performance.

---
