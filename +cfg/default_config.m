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
    % Set gravitational acceleration (g) and friction model
    config.phys.g = 9.81;   % [m/s^2] Acceleration due to gravity
    config.phys.Cf = 0;     % [optional] Bottom friction coefficient (legacy, kept for compatibility)
    % IMPORTANT: dry_tolerance is used for well-balanced hydrostatic reconstruction and safe division.
    config.phys.dry_tolerance = 1e-6; % [m] Water depth below which a cell is considered dry
    
    % Default friction model: no friction
    config.phys.friction_model = @friction.no_friction;
    
    % Additional parameters for advanced friction models
    config.phys.darcy_f = 0.02;              % Default Darcy friction factor
    config.phys.ks = 0.001;                  % Default roughness height [m]
    config.phys.kinematic_viscosity = 1e-6;  % Default kinematic viscosity [m²/s]
    config.phys.cw_iterations = 20;          % Default max iterations for Colebrook-White
    config.phys.cw_tolerance = 1e-6;         % Default tolerance for Colebrook-White

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Reconstruction Configuration ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Default is characteristic reconstruction
    config.reconstruction.method = 'muscl';
    config.reconstruction.order = 2;
    config.reconstruction.limiter = 'minmod';
    config.reconstruction.handle = @reconstruct.muscl;
    config.reconstruction.characteristic = true; % Default: characteristic-wise for all methods
    config.reconstruction.mp5_mode = 'characteristic';
    config.reconstruction.ppm_mode = 'characteristic';
    config.reconstruction.weno_mode = 'characteristic';
    config.reconstruction.uno2_mode = 'characteristic';
    config.reconstruction.limiter_handle = reconstruct.limiters.limiter_selector(config.reconstruction.limiter);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Domain and Mesh ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define computational domain and mesh resolution
    config.mesh.domain.xmin = 0;     % [m] Domain start
    config.mesh.domain.xmax = 100;   % [m] Domain end
    config.mesh.N = 200;             % [integer] Number of finite volume cells
    % DO NOT set config.mesh.dx here! It will be computed automatically in validation.
    % Default bathymetry: flat bottom (can be overridden)
    % Now, bathy.flat returns a constant bottom elevation (default 0.0, z=0 datum)
    if ~isfield(config, 'bathy_params'), config.bathy_params = struct(); end
    config.bathy_params.flat_elevation = 0.0; % Default flat bottom elevation (z=0 datum)
    config.mesh.h_fun = @(x) zeros(size(x)); % [m] (Legacy, not used if bathyHandle is set)
    % Default bathymetry handle (used by well-balanced schemes):
    config.bathyHandle = @bathy.flat; % Must accept (cfg, x) and return bathymetry at x

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
    if ~isfield(config, 'numerics')
        config.numerics = struct();
    end
    config.numerics.epsilon = 1e-10; % Default numerical epsilon for stability/zero division
    if isfield(config, 'reconstruct') && isfield(config.reconstruct, 'order') && config.reconstruct.order > 1
        config.numerics.rhs_handle = @core.rhs_nsw_high_order;  % High-order RHS
    else
        config.numerics.rhs_handle = @core.rhs_nsw_1st_order;  % 1st order RHS (default)
    end
    config.numerics.flux_handle = @flux.FVCF;             % [function handle] Default numerical flux (FVCF)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Boundary Conditions ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set default boundary condition function handles
    % Valid BC options include: @bc.wall, @bc.generating, @bc.periodic
    config.bc.left_handle = @bc.wall;   % [function handle] Left boundary (solid wall)
    config.bc.right_handle = @bc.wall;  % [function handle] Right boundary (solid wall)

    % Number of ghost cells depends on reconstruction order
    if isfield(config, 'reconstruct') && isfield(config.reconstruct, 'order')
        config.bc.num_ghost_cells = max(1, config.reconstruct.order); 
    else
        config.bc.num_ghost_cells = 1;  % Default: 1 ghost cell (1st order)
    end

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