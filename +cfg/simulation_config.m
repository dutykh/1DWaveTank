%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +cfg/simulation_config.m
%
% Purpose:
%   Sets up the specific configuration for a simulation run in the 1DWaveTank code.
%   This function allows the user to select and customize an experiment setup
%   (initial/boundary conditions, numerics, domain, etc.) by overriding the
%   defaults from default_config.m. The result is a complete config structure
%   ready for use by the solver.
%
% Syntax:
%   config = simulation_config()
%
% Inputs:
%   (none)
%
% Outputs:
%   config - [struct] Complete configuration structure for the selected simulation run.
%
% Dependencies:
%   Calls cfg.default_config() for baseline values. Requires all referenced
%   function handles (fluxes, BCs, ICs, etc.) to exist in the codebase.
%
% References:
%   - See 1DWaveTank UserGuide.md for structure and usage.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function config = simulation_config()

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print status for user clarity (optional, can be commented)  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('--- Setting up Simulation Configuration ---\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Load Default Configuration ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Start with the baseline parameters defined in default_config.m
    config = cfg.default_config();
    fprintf('Default config loaded. Overriding for specific experiment...\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Experiment Selection ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Choose a predefined setup or define a custom one below
    % Available setups: 'flat_rest', 'flat_gaussian', 'flat_wave_gen'
    % To change the simulation run, modify the 'experiment_setup' variable below.
    experiment_setup = 'flat_wave_gen'; % CHANGE THIS TO SELECT SETUP
    config.experiment_setup = experiment_setup; % Store the chosen setup name in config

    fprintf('Selected experiment setup: %s\n', experiment_setup);

    % ======================================================================
    % --- Common Settings (can be overridden by specific setups below) ---
    % ======================================================================
    % These are the default values for domain, mesh, numerics, etc.
    % They can be overridden in the switch block for each experiment.

    % --- Domain and Mesh ---
    config.domain.xmin = 0.0;              % [m] Left endpoint of domain
    config.domain.xmax = 20.0;             % [m] Right endpoint of domain
    config.mesh.N      = 500;              % [integer] Number of spatial cells
    config.param.H0    = 0.50;             % [m] Undisturbed water depth

    % --- Model and Numerics ---
    config.model = @core.rhs_nsw_1st_order;        % [function handle] RHS function (1st order FV)
    config.numFlux = @flux.PVM;                    % [function handle] Numerical flux
    config.reconstructopenion = [];                % [empty/struct] No reconstruction (1st order)
    config.timeStepper = @time.integrate_euler_adaptive; % [function handle] Time integration wrapper
    % config.timeStepper = @time.integrate_matlab_ode; % Alternative: MATLAB ODE
    % config.time.matlab_solver = 'ode45';           % MATLAB ODE solver
    % config.time.ode_options = odeset();            % MATLAB ODE options
    % config.time.AbsTol = 1e-4;                     % Absolute tolerance for MATLAB ODE
    % config.time.RelTol = 1e-4;                     % Relative tolerance for MATLAB ODE
    % config.time.show_progress_bar = true;          % Show progress bar for MATLAB ODE
    config.time.num_progress_reports = 10;         % [integer] Number of progress updates
    config.time.cfl = 0.95;                        % [unitless] CFL number (not used by MATLAB ODE)

    % --- Run Control ---
    config.t0 = 0.0;                               % [s] Simulation start time
    config.tEnd = 10.0;                             % [s] Simulation end time
    config.vis.dt_plot = 0.1;                      % [s] Output interval for visualization/saving
    config.vis.plot_velocity = true;               % [logical] Plot velocity in a subpanel
    config.vis.show_legend = false;                % [logical] Show legend in wave tank plot
    config.tspan = config.t0:config.vis.dt_plot:config.tEnd; % [vector] Output time points
    config.save_results = false;                   % [logical] Save results by default
    config.output_path = 'results/';               % [char] Output directory

    % ======================================================================
    % --- Specific Experiment Setups ---
    % ======================================================================
    % Use a switch block to override settings for each experiment.

    switch experiment_setup
        case 'flat_rest'
            % Simple case: flat bottom, lake at rest, wall boundaries.
            % Useful for testing stability.
            config.caseName = 'flat_rest_L20m_H0.5m_N400';

            % Bathymetry
            config.bathyHandle = @bathy.flat;

            % Initial Condition
            config.ic_handle = @ic.lake_at_rest;

            % Boundary Conditions
            config.bc.left.handle = @bc.wall;
            config.bc.right.handle = @bc.wall;

            % (Other settings can be added/overridden here)

        case 'flat_gaussian'
            % Gaussian bump in water surface, open boundaries.
            config.caseName = 'flat_gaussian_L20m_H0.5m_N500';

            config.bathyHandle = @bathy.flat;
            config.ic_handle = @ic.gaussian_bump;

            % Open boundaries
            config.bc.left.handle = @bc.open;
            config.bc.right.handle = @bc.open;

            % Parameters for initial condition (Gaussian bump)
            config.ic_param.a = 0.1;    % [m] Amplitude
            config.ic_param.lambda = 2; % [m] Width
            config.ic_param.x0 = 10.0;  % [m] Center location
            config.ic_param.H0 = config.param.H0; % [m] Reference depth

        case 'flat_wave_gen'
            % Sine wave generated at left boundary, wall at right.
            config.caseName = 'flat_wave_gen_L20m_H0.5m_N500';

            config.bathyHandle = @bathy.flat;
            config.ic_handle = @ic.lake_at_rest;

            % Generating BC at left, wall at right
            config.bc.left.handle = @bc.generating;
            config.bc.right.handle = @bc.wall;
            config.bc.left.param.a = 0.1;    % [m] Amplitude
            config.bc.left.param.T = 2*pi;   % [s] Period
            config.bc.right.param = struct(); % No params needed for wall

            % (Other settings can be added/overridden here)

        case 'flat_solitary'
            % Solitary wave on flat bottom, wall boundaries.
            config.caseName = 'flat_solitary_L20m_H0.5m_A0.2m_N500';
            config.L = 20.0;                               % [m] Tank length (matches default)
            config.H0 = 0.5;                               % [m] Still water depth (matches default)
            config.Nx = 500;                               % [-] Number of grid points (matches default)
            config.param.a = 0.2;                          % [m] Solitary wave amplitude
            config.bathyHandle = @bathy.flat;          % Bathymetry function handle
            config.ic_handle = @ic.solitary_wave;          % Initial condition handle
            % Define boundary conditions explicitly (using defaults)
            config.bc.left.handle = @bc.wall;
            config.bc.right.handle = @bc.wall;

        case 'periodic_solitary'
            config.bc.left.handle = @bc.periodic;
            config.bc.right.handle = @bc.periodic;
            config.bathyHandle = @bathy.flat;
            config.ic_handle = @ic.solitary_wave;
            % Use default IC parameters (h0=0.5, a=0.2) unless overridden
            % Get potentially overridden h0 from default_config or command line
            if isfield(config.param, 'h0'); h0 = config.param.h0; else h0 = 0.5; end
            % Get potentially overridden a from default_config or command line
            if isfield(config.param, 'a'); a = config.param.a; else a = 0.2; end
            config.caseName = sprintf('periodic_solitary_L%.0fm_H%.1fm_A%.1fm_N%d', ...
                                    config.domain.xmax - config.domain.xmin, h0, a, config.mesh.N);
            fprintf('--- Configuration for periodic solitary wave ---\n');
            fprintf('      Bathymetry: Flat\n');
            fprintf('      Initial Condition: Solitary Wave (h0=%.2f, a=%.2f)\n', h0, a);
            fprintf('      Boundary Conditions: Periodic\n');

        otherwise
            error('Unknown experiment setup: %s', experiment_setup);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Mesh Generation ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate mesh spacing and cell centers based on domain and N
    config.mesh.xmin = config.domain.xmin;
    config.mesh.xmax = config.domain.xmax;
    config.mesh.N = config.mesh.N; % Ensure consistency
    config.mesh.dx = (config.mesh.xmax - config.mesh.xmin) / config.mesh.N;
    config.mesh.x = linspace(config.mesh.xmin, config.mesh.xmax, config.mesh.N+1); % Cell edges
    config.mesh.xc = 0.5*(config.mesh.x(1:end-1) + config.mesh.x(2:end));         % Cell centers

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%a%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Print Summary (optional) ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('--- Simulation Configuration Summary ---\n');
    fprintf('  Case: %s\n', config.caseName);
    fprintf('  Domain: [%.2f, %.2f] m, N = %d\n', config.domain.xmin, config.domain.xmax, config.mesh.N);
    fprintf('  Bathymetry: %s\n', func2str(config.bathyHandle));
    fprintf('  Initial Condition: %s\n', func2str(config.ic_handle));
    fprintf('  BCs: left = %s, right = %s\n', func2str(config.bc.left.handle), func2str(config.bc.right.handle));
    fprintf('  Flux: %s\n', func2str(config.numFlux));
    fprintf('  Time Span: [%.2f, %.2f] s\n', config.t0, config.tEnd); % Removed CFL as it's not directly used
    fprintf('-------------------------------------------\n');

end