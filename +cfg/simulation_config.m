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
% Author: Dr. Denys Dutykh
% Date:   April 26, 2025
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

    % Ensure h0 and L are always present for bathymetry and visualization compatibility
    % (These may be overwritten in the scenario block, but this guarantees presence)
    if isfield(config, 'param') && isfield(config.param, 'H0')
        config.h0 = config.param.H0;
    elseif isfield(config, 'H0')
        config.h0 = config.H0;
    end
    if isfield(config, 'domain') && isfield(config.domain, 'xmax') && isfield(config.domain, 'xmin')
        config.L = config.domain.xmax - config.domain.xmin;
    elseif isfield(config, 'L')
        config.L = config.L;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Experiment Selection ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Choose a predefined setup or define a custom one below
    % Available setups: 'flat_rest', 'flat_gaussian', 'flat_wave_gen', 'flat_solitary', 'periodic_solitary', 'dam_break'
    % To change the simulation run, modify the 'experiment_setup' variable below.
    experiment_setup = 'gaussian_bump_rest'; % CHANGE THIS TO SELECT SETUP
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
    % config.param.H0    = 0.50;             % [m] Undisturbed water depth (legacy, not used for flat bottom)

    % --- Model and Numerics ---
    config.model = @core.rhs_nsw_1st_order;        % [function handle] RHS function (1st order FV)
    config.numFlux = @flux.PVM;                    % [function handle] Numerical flux
    config.reconstructopenion = [];                % [empty/struct] No reconstruction (1st order)
    % config.timeStepper = @time.integrate_euler_adaptive; % [function handle] Time integration wrapper
    config.timeStepper = @time.integrate_matlab_ode; % Alternative: MATLAB ODE
    config.time.matlab_solver = 'ode113';           % MATLAB ODE solver
    config.time.ode_options = odeset();            % MATLAB ODE options
    config.time.AbsTol = 1e-4;                     % Absolute tolerance for MATLAB ODE
    config.time.RelTol = 1e-4;                     % Relative tolerance for MATLAB ODE
    config.time.show_progress_bar = true;          % Show progress bar for MATLAB ODE
    config.time.num_progress_reports = 10;         % [integer] Number of progress updates
    config.time.cfl = 0.95;                        % [unitless] CFL number (not used by MATLAB ODE)

    % --- Run Control ---
    config.t0 = 0.0;                               % [s] Simulation start time
    config.tEnd = 10.0;                            % [s] Simulation end time
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
            
            % Example of how to use friction models
            % Comment/uncomment to enable/disable friction
            % Option 1: Chézy model
            % config.phys.friction_model = friction.friction_selector('chezy');
            % config.phys.chezy_C = 50;  % Chézy coefficient [m^(1/2)/s]
            % Option 2: Manning model
            % config.phys.friction_model = friction.friction_selector('manning');
            % config.phys.manning_n = 0.03;  % Manning coefficient [s/m^(1/3)]
            % Option 3: Darcy-Weisbach with constant f
            % config.phys.friction_model = friction.friction_selector('darcy_weisbach');
            % config.phys.darcy_f = 0.02;  % Constant Darcy friction factor
            % Option 4: Darcy-Weisbach with Colebrook-White
            % config.phys.friction_model = friction.friction_selector('darcy_weisbach');
            % config.phys.f_calculation = 'colebrook_white';
            % config.phys.ks = 0.001;  % Equivalent sand roughness [m]
            % config.phys.kinematic_viscosity = 1e-6;  % Kinematic viscosity [m²/s]

            % Bathymetry (flat bottom, z_b = 0 by default)
            config.bathyHandle = @bathy.flat;
            % Optionally set bottom elevation (default is 0.0)
            config.bathy_params.flat_elevation = 0.0; % [m] Bottom elevation (z=0 datum)

            % Initial Condition (lake at rest)
            config.h0 = 0.5; % [m] Still Water Level (SWL) elevation z_s
            config.ic_handle = @ic.lake_at_rest;

            % Boundary Conditions
            config.bc.left.handle = @bc.wall;
            config.bc.right.handle = @bc.wall;

            % (Other settings can be added/overridden here)

        case 'flat_gaussian'
            % Gaussian bump in water surface, open boundaries.
            config.caseName = 'flat_gaussian_L20m_H0.5m_N500';

            % Example of how to use friction models
            % Comment/uncomment to enable/disable friction
            % Option 1: Chézy model
            % config.phys.friction_model = friction.friction_selector('chezy');
            % config.phys.chezy_C = 50;  % Chézy coefficient [m^(1/2)/s]
            % Option 2: Manning model
            % config.phys.friction_model = friction.friction_selector('manning');
            % config.phys.manning_n = 0.03;  % Manning coefficient [s/m^(1/3)]
            % Option 3: Darcy-Weisbach with constant f
            % config.phys.friction_model = friction.friction_selector('darcy_weisbach');
            % config.phys.darcy_f = 0.02;  % Constant Darcy friction factor
            % Option 4: Darcy-Weisbach with Colebrook-White
            % config.phys.friction_model = friction.friction_selector('darcy_weisbach');
            % config.phys.f_calculation = 'colebrook_white';
            % config.phys.ks = 0.001;  % Equivalent sand roughness [m]
            % config.phys.kinematic_viscosity = 1e-6;  % Kinematic viscosity [m²/s]

            config.bathyHandle = @bathy.flat;
            config.bathy_params.flat_elevation = 0.0; % [m] Bottom elevation (z=0 datum)
            config.ic_handle = @ic.gaussian_bump;

            % Open boundaries
            config.bc.left.handle = @bc.open;
            config.bc.right.handle = @bc.open;

            % Parameters for initial condition (Gaussian bump)
            config.ic_param.a = 0.1;    % [m] Amplitude
            config.ic_param.lambda = 2; % [m] Width
            config.ic_param.x0 = 10.0;  % [m] Center location
            config.ic_param.H0 = config.param.H0; % [m] Reference depth

            % Set reconstruction to MUSCL (2nd order)
            config.reconstruct.method = 'muscl';
            config.reconstruct.handle = @reconstruct.muscl;
            config.reconstruct.order = 2;
            config.reconstruct.limiter = @reconstruct.limiters.vanleer; % or minmod, superbee, etc.
            config.reconstruct.theta = 1/3; % Third-order accuracy in smooth regions

            % Update RHS to high-order version
            config.model = @core.rhs_nsw_high_order;

            % Higher-order time integrator for matching temporal accuracy
            config.timeStepper = @time.integrate_ssp2_adaptive; % or rk4, etc.

        case 'flat_wave_gen'
            % Sine wave generated at left boundary, wall at right.
            config.caseName = 'flat_wave_gen_L20m_H0.5m_N500_MUSCL'; % Updated name
            
            % Example of how to use friction models
            % Using Darcy-Weisbach with Colebrook-White formula
            % config.phys.friction_model = friction.friction_selector('darcy_weisbach');
            % config.phys.f_calculation = 'colebrook_white';
            % config.phys.ks = 0.002;  % Equivalent sand roughness [m]
            % config.phys.kinematic_viscosity = 1e-6;  % Kinematic viscosity [m²/s]
            
            % Other friction model options (commented out)
            config.phys.friction_model = friction.friction_selector('no_friction');
            % config.phys.manning_n = 0.03;  % Manning coefficient [s/m^(1/3)]

            config.bathyHandle = @bathy.flat;
            config.bathy_params.flat_elevation = 0.0; % [m] Bottom elevation (z=0 datum)
            config.h0 = 0.5; % [m] Still Water Level (SWL) elevation z_s
            config.ic_handle = @ic.lake_at_rest;

            % Generating BC at left, wall at right
            config.bc.left.handle = @bc.generating;
            config.bc.right.handle = @bc.wall;
            config.bc.left.param.a = 0.1;    % [m] Amplitude
            config.bc.left.param.T = 2*pi;   % [s] Period
            config.bc.right.param = struct(); % No params needed for wall

            % --- High-Order Configuration ---
            % Set reconstruction to MUSCL (2nd order)
            config.reconstruct.method = 'muscl';
            config.reconstruct.handle = @reconstruct.muscl;
            config.reconstruct.order = 2;
            config.reconstruct.limiter = @reconstruct.limiters.vanleer; % Using van Leer limiter
            config.reconstruct.theta = 1/3; % Third-order accuracy in smooth regions
            
            % Update RHS to high-order version
            config.numerics.rhs_handle = @core.rhs_nsw_high_order;
            
            % Higher-order time integrator for matching temporal accuracy
            config.timeStepper = @time.integrate_ssp2_adaptive; % Using SSP2 adaptive

        case 'flat_solitary'
            % Solitary wave on flat bottom, wall boundaries.
            config.caseName = 'flat_solitary_L20m_H0.5m_A0.2m_N500';
            
            % Example of how to use friction models
            % Comment/uncomment to enable/disable friction
            % Option 1: Chézy model
            % config.phys.friction_model = friction.friction_selector('chezy');
            % config.phys.chezy_C = 50;  % Chézy coefficient [m^(1/2)/s]
            % Option 2: Manning model
            % config.phys.friction_model = friction.friction_selector('manning');
            % config.phys.manning_n = 0.03;  % Manning coefficient [s/m^(1/3)]
            % Option 3: Darcy-Weisbach with constant f
            % config.phys.friction_model = friction.friction_selector('darcy_weisbach');
            % config.phys.darcy_f = 0.02;  % Constant Darcy friction factor
            % Option 4: Darcy-Weisbach with Colebrook-White
            % config.phys.friction_model = friction.friction_selector('darcy_weisbach');
            % config.phys.f_calculation = 'colebrook_white';
            % config.phys.ks = 0.001;  % Equivalent sand roughness [m]
            % config.phys.kinematic_viscosity = 1e-6;  % Kinematic viscosity [m²/s]
            config.L = 20.0;                               % [m] Tank length (matches default)
            config.H0 = 0.5;                               % [m] Still water depth (matches default)
            config.Nx = 500;                               % [-] Number of grid points (matches default)
            config.param.a = 0.2;                          % [m] Solitary wave amplitude
            config.bathyHandle = @bathy.flat;          % Bathymetry function handle (returns bottom elevation)
            config.bathy_params.flat_elevation = 0.0;  % [m] Bottom elevation (z=0 datum)
            config.ic_handle = @ic.solitary_wave;          % Initial condition handle
            % Define boundary conditions explicitly (using defaults)
            config.bc.left.handle = @bc.wall;
            config.bc.right.handle = @bc.wall;

        case 'periodic_solitary'
            config.bc.left.handle = @bc.periodic;
            config.bc.right.handle = @bc.periodic;
            
            % Example of how to use friction models
            % Comment/uncomment to enable/disable friction
            % Option 1: Chézy model
            % config.phys.friction_model = friction.friction_selector('chezy');
            % config.phys.chezy_C = 50;  % Chézy coefficient [m^(1/2)/s]
            % Option 2: Manning model
            % config.phys.friction_model = friction.friction_selector('manning');
            % config.phys.manning_n = 0.03;  % Manning coefficient [s/m^(1/3)]
            % Option 3: Darcy-Weisbach with constant f
            % config.phys.friction_model = friction.friction_selector('darcy_weisbach');
            % config.phys.darcy_f = 0.02;  % Constant Darcy friction factor
            % Option 4: Darcy-Weisbach with Colebrook-White
            % config.phys.friction_model = friction.friction_selector('darcy_weisbach');
            % config.phys.f_calculation = 'colebrook_white';
            % config.phys.ks = 0.001;  % Equivalent sand roughness [m]
            % config.phys.kinematic_viscosity = 1e-6;  % Kinematic viscosity [m²/s]
            config.bathyHandle = @bathy.flat;
            config.bathy_params.flat_elevation = 0.0; % [m] Bottom elevation (z=0 datum)
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

        case 'dam_break'
            % Dam break problem on a flat bottom, wall boundaries.
            % Important test case for shock capturing.
            config.caseName = 'dam_break_L20m_H0.8-0.5m_N500_MP5Char_RK4'; % Updated name
            
            % --- Mesh Override ---
            config.mesh.N = 500; % Ensure correct dx calculation during validation
            % Force removal of dx to prevent validation warning
            if isfield(config.mesh, 'dx'); config.mesh = rmfield(config.mesh, 'dx'); end
            
            % --- Final Time Override ---
            config.tEnd = 3.0;           % [s] Shorter time for typical dam break evolution
            
            % --- Domain and Mesh ---
            config.domain.xmin = 0.0;    % [m] Left endpoint of domain
            config.domain.xmax = 20.0;   % [m] Right endpoint of domain
            config.mesh.domain.xmin = config.domain.xmin;
            config.mesh.domain.xmax = config.domain.xmax;
            
            % --- Bathymetry (flat bottom) ---
            config.bathyHandle = @bathy.flat;
            config.param.H0 = 0.0;       % Reference depth for bathymetry (not used)
            
            % --- Dam Break Initial Condition ---
            config.ic_handle = @ic.dam_break;
            config.dam_break.h_L = 0.8;  % [m] Water depth left of dam
            config.dam_break.h_R = 0.5;  % [m] Water depth right of dam
            config.dam_break.x_dam = 10.0; % [m] Dam position at domain center
            
            % --- Boundary Conditions (walls on both sides) ---
            config.bc.left.handle = @bc.wall;
            config.bc.right.handle = @bc.wall;
            
            % --- Numerical Model ---
            % config.model = @core.rhs_nsw_1st_order; % Use 1st order
            config.model = @core.rhs_nsw_high_order; % USE HIGH-ORDER RHS for reconstruction
            config.numFlux = @flux.SLAU;             % SLAU flux (robust for shocks)
            
            % --- Reconstruction Settings --- 
            % Use MUSCL with van Albada limiter
            config.reconstruction = struct();
            config.reconstruction.method = 'muscl';
            config.reconstruction.limiter = 'vanalbada';
            config.bc.num_ghost_cells = 2; % MUSCL requires at least 2 ghost cells
            config.reconstruction.handle = reconstruct.reconstruct_selector(config.reconstruction.method);
            config.reconstruct = config.reconstruction; % Ensure compatibility with core solver
            
            % --- Time Integration --- 
            config.timeStepper = @time.integrate_ssp2_adaptive; 
            % config.timeStepper = @time.integrate_rk4_adaptive; % Use RK4 for 4th order
            config.cfl_target = 0.95;                   % CFL for RK4 (can be higher than SSP)
            
            % --- Visualization ---
            config.vis.dt_plot = 0.1;                 % [s] Output interval
            config.vis.plot_velocity = true;
            config.tspan = config.t0:config.vis.dt_plot:config.tEnd; % Use tEnd, not T

            fprintf('--- Configuration for dam break problem ---\n');
            fprintf('      Bathymetry: Flat\n');
            fprintf('      Initial Condition: Dam break (h_L=%.2f, h_R=%.2f)\n', config.dam_break.h_L, config.dam_break.h_R);
            fprintf('      Boundary Conditions: Wall on both sides\n');
            fprintf('      Time Integration: Adaptive Euler with CFL=%.2f\n', config.time.cfl);

        case 'dry_dam_break'
            % Dry dam break test case: water only on left side, dry on right
            config.caseName = 'dry_dam_break_L20m_H0.5m_N500';
            
            % --- Domain setup ---
            config.domain.xmin = 0.0;
            config.domain.xmax = 20.0;
            config.mesh.domain.xmin = config.domain.xmin;
            config.mesh.domain.xmax = config.domain.xmax;
            config.mesh.N = 500;
            % config.param.H0 = 0.5; % [m] (Legacy, not used for flat bottom convention)
            
            % --- Bathymetry ---
            config.bathyHandle = @bathy.flat;
            
            % --- Initial Condition ---
            config.ic_handle = @ic.dry_dam_break;
            
            % --- Boundary Conditions ---
            config.bc.left.handle = @bc.wall;
            config.bc.right.handle = @bc.wall;
            
            % --- Numerical scheme ---
            config.model = @core.rhs_nsw_high_order; % High-order FV method
            config.numFlux = @flux.HLLE;               % FVCF flux for shock capturing
            config.reconstruction.handle = @reconstruct.muscl; % MUSCL reconstruction
            config.reconstruction.limiter = 'vanleer';         % van Leer limiter
            config.reconstruction.characteristic = true;      % Component-wise reconstruction
            config.bc.num_ghost_cells = 2; % MUSCL requires at least 2 ghost cells
            
            % --- Time integration ---
            config.timeStepper = @time.integrate_euler_adaptive;  % Adaptive Euler
            config.time.cfl = 0.95;                  % CFL number

            % config.time.matlab_solver = 'ode113';          % MATLAB ODE solver
            % config.time.ode_options = odeset();            % MATLAB ODE options
            % config.time.AbsTol = 1e-4;                     % Absolute tolerance for MATLAB ODE
            % config.time.RelTol = 1e-4;                     % Relative tolerance for MATLAB ODE
            % config.time.show_progress_bar = true;          % Show progress bar for MATLAB ODE
            
            % --- Parameters for dry dam break ---
            config.dry_dam_break.h_L = 0.5;          % Water depth on left side [m]
            config.dry_dam_break.x_dam = 10.0;       % Dam location at x=10m
            
            % --- Run Control ---
            config.t0 = 0.0;                         % Start time [s]
            config.tEnd = 1.0;                       % End time [s]
            config.vis.dt_plot = 0.1;                % Output interval [s]
            config.tspan = config.t0:config.vis.dt_plot:config.tEnd;
            
        case 'gaussian_bump_rest'
            % Underwater Gaussian bump, lake at rest, wall boundaries (Scenario 1)
            config.caseName = 'gaussian_bump_rest_L20m_H0.5m_N500';

            % Domain (Set within mesh structure)
            config.mesh.domain.xmin = 0.0;
            config.mesh.domain.xmax = 20.0;
            config.mesh.N      = 500;
            config.param.H0    = 0.50;

            % Guarantee h0 and L for bathy.gaussian_bump (for all code paths)
            config.h0 = config.param.H0;
            config.L  = config.mesh.domain.xmax - config.mesh.domain.xmin; % Use mesh.domain

            % Bathymetry: Gaussian bump
            config.bathyHandle = @bathy.gaussian_bump;
            config.bathy_bump_center = (config.mesh.domain.xmin + config.mesh.domain.xmax)/2;   % Use mesh.domain
            config.bathy_bump_height = 0.15 * config.param.H0;                        % Bump height
            config.bathy_bump_width  = (config.mesh.domain.xmax - config.mesh.domain.xmin)/15;  % Use mesh.domain

            % Initial Condition: lake at rest
            config.ic_handle = @ic.lake_at_rest;

            % Boundary Conditions: walls
            config.bc.left.handle = @bc.wall;
            config.bc.right.handle = @bc.wall;

            % Numerical scheme: 1st order FV (no reconstruction)
            config.model = @core.rhs_nsw_1st_order;   % 1st order FV method
            config.numFlux = @flux.HLLE;              % Riemann solver (HLLE)
            config.reconstruction = [];               % No reconstruction
            % config.bc.num_ghost_cells = 1;          % 1st order usually needs 1 ghost cell (set if required)

            % Time integration
            config.timeStepper = @time.integrate_euler_adaptive;
            config.time.cfl = 0.5; % More conservative for stability over bump

            % Run control
            config.t0 = 0.0;
            config.tEnd = 10.0;
            config.vis.dt_plot = 0.1;
            config.tspan = config.t0:config.vis.dt_plot:config.tEnd;

        otherwise
            error('Unknown experiment setup: %s', experiment_setup);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Mesh Generation ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The following mesh property assignments are now commented out to prevent
    % inconsistent dx warnings. Mesh properties will be set by validate_config.
    % config.mesh.xmin = config.domain.xmin;
    % config.mesh.xmax = config.domain.xmax;
    % config.mesh.N = config.mesh.N; % Ensure consistency
    % config.mesh.dx = (config.mesh.xmax - config.mesh.xmin) / config.mesh.N;
    % config.mesh.x = linspace(config.mesh.xmin, config.mesh.xmax, config.mesh.N+1); % Cell edges
    % config.mesh.xc = 0.5*(config.mesh.x(1:end-1) + config.mesh.x(2:end));         % Cell centers

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