%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +cfg/default_config.m
%
% Purpose:
%   Provides a default set of configuration parameters for the 1DWaveTank code.
%   This function returns a structure containing baseline values for all
%   simulation parameters (physical, mesh, numerics, BCs, ICs, visualization).
%   These defaults can be overridden by scenario-specific or user settings.
%
% Syntax:
%   config = default_config()
%
% Inputs:
%   (none)
%
% Outputs:
%   config - [struct] Complete configuration structure with default values for
%            all major simulation settings and function handles.
%
% Dependencies:
%   None directly, but function handles (e.g., for fluxes, BCs, ICs) must
%   correspond to available files in the codebase.
%
% References:
%   - See 1DWaveTank UserGuide.md for structure and usage.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function config = default_config()

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print status for user clarity (optional, can be commented)  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Loading default configuration...\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Physical Parameters ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set gravitational acceleration (g) and bottom friction (Cf)
    config.phys.g = 9.81;   % [m/s^2] Acceleration due to gravity
    config.phys.Cf = 0;     % [optional] Bottom friction coefficient (e.g., Manning's n or Chezy C). Default 0 = no friction.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Domain and Mesh ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define computational domain and mesh resolution
    config.mesh.domain = [0, 100];   % [m] Computational domain [x_min, x_max]
    config.mesh.N = 200;             % [integer] Number of finite volume cells
    % Default bathymetry: flat bottom (can be overridden)
    config.mesh.h_fun = @(x) zeros(size(x)); % [m] Bathymetry function (depth positive downward from z=0)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Time Integration ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set time span, CFL number, plotting interval, and default integrator
    config.time.t_span = [0, 30.0];  % [s] Simulation time interval [t_start, t_end]
    config.time.dt_init = 0.01;      % [s] Initial time step guess for adaptive methods
    config.time.cfl = 0.5;           % [unitless] CFL number for adaptive time stepping
    config.time.dt_plot = 0.1;       % [s] Time interval for plotting/saving output
    config.time.adaptive_tol = 1e-4; % Tolerance for adaptive RK methods
    % Default time integrator handle (can be overridden by user/scenario)
    config.time.integrator = @time.integrate_euler_adaptive;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Numerics ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set handles for the numerical RHS and flux functions
    config.numerics.rhs_handle = @core.rhs_nsw_1st_order; % [function handle] RHS for 1st order FV scheme
    config.numerics.flux_handle = @flux.FVCF;             % [function handle] Default numerical flux (FVCF)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Boundary Conditions ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set default boundary condition function handles
    % Valid BC options include: @bc.wall, @bc.generating, @bc.periodic
    config.bc.left_handle = @bc.wall;   % [function handle] Left boundary (solid wall)
    config.bc.right_handle = @bc.wall;  % [function handle] Right boundary (solid wall)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Initial Conditions ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set default initial condition function handle and parameters
    config.ic_handle = @ic.lake_at_rest; % [function handle] IC: still water everywhere
    config.ic_param = struct();          % [struct] Parameters for IC (empty by default)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Visualization ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set plotting options for simulation output
    config.vis.plot_vars = {'H', 'U'}; % [cell array] Variables to plot ('H', 'U', 'HU')
    config.vis.ylim_margin = 0.1;      % [unitless] Margin factor for y-axis limits in plots

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Output Control ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    config.save_results = false;        % [logical] Do not save results by default

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print status for user clarity (optional, can be commented)  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Default configuration loaded.\n');

end