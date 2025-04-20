function config = default_config()
% DEFAULT_CONFIG Provides a default set of configuration parameters.
%
%   This function defines baseline parameters for the simulation, which
%   can be overridden by specific experiment setups in simulation_config.m.
%
%   Outputs:
%     config - Structure containing default configuration parameters.

    fprintf('Loading default configuration...\n');

    % --- Physical Parameters ---
    % Define physical constants and parameters
    config.phys.g = 9.81; % [m/s^2] Acceleration due to gravity
    config.phys.Cf = 0;   % [optional] Bottom friction coefficient (e.g., Manning's n or Chezy C). Default 0 = no friction.

    % --- Domain and Mesh ---
    % Define computational domain and mesh parameters
    config.mesh.domain = [0, 100]; % [m] Computational domain [x_min, x_max]
    config.mesh.N = 200;        % Number of finite volume cells
    % Define bathymetry function h(x) (depth positive downwards from z=0)
    config.mesh.h_fun = @(x) zeros(size(x)); % Default: flat bottom at z=0

    % --- Time Integration ---
    % Define time integration parameters
    config.time.t_span = [0, 30.0]; % [s] Simulation time interval [t_start, t_end]
    config.time.cfl = 0.5;          % CFL number for adaptive time stepping
    config.time.dt_plot = 0.1;      % [s] Time interval for plotting/saving output
    % Default time integrator handle (can be overridden)
    config.time.integrator = @time.integrate_euler_adaptive;

    % --- Numerics ---
    % Define numerical method parameters
    % Default RHS function handle for the governing equations
    config.numerics.rhs_handle = @core.rhs_nsw_1st_order;
    % Default numerical flux function handle
    config.numerics.flux_handle = @flux.FVCF;

    % --- Boundary Conditions ---
    % Define boundary condition parameters
    % Default handles for left and right boundary conditions
    config.bc.left_handle = @bc.wall;
    config.bc.right_handle = @bc.wall;

    % --- Initial Conditions ---
    % Define initial condition parameters
    % Default initial condition handle
    config.prob.ic_handle = @ic.lake_at_rest;

    % --- Visualization ---
    % Define visualization parameters
    config.vis.plot_vars = {'H', 'U'}; % Variables to plot ('H', 'U', 'HU')
    config.vis.ylim_margin = 0.1;      % Margin factor for y-axis limits in plots

    fprintf('Default configuration loaded.\n');
end