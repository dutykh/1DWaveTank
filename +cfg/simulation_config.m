function config = simulation_config()
% SIMULATION_CONFIG Sets up the specific configuration for a simulation run.
%
%   This function defines different experimental setups (cases) and returns
%   the complete configuration structure 'config' for the selected experiment.
%
%   Workflow:
%   1. Loads the default configuration from cfg.default_config.
%   2. Selects an experiment setup based on the 'experiment_setup' variable.
%   3. Uses a switch statement to apply specific overrides for the chosen experiment
%      (e.g., different initial conditions, boundary conditions, domain, time).
%   4. Calculates mesh properties (dx, x, xc) based on the final domain and N.
%   5. Returns the fully configured 'config' structure.
%
%   To change the simulation run, modify the 'experiment_setup' variable below.
%
%   Outputs:
%     config - The complete configuration structure for the selected simulation run.

% config.simulation_config: Central configuration file for 1DWaveTank simulations.
%   config = config.simulation_config() returns a structure 'config' containing all
%   parameters and function handles for a specific simulation setup.
%   Modify this file to set up different test cases.

fprintf('--- Setting up Simulation Configuration ---\n');

% --- Load Default Configuration ---
% Start with the baseline parameters defined in default_config.m
config = cfg.default_config();
fprintf('Default config loaded. Overriding for specific experiment...\n');

% --- Experiment Selection ---
% Choose a predefined setup or define a custom one below
% Available setups: 'flat_rest', 'flat_gaussian', 'flat_wave_gen'
% To change the simulation run, modify the 'experiment_setup' variable below.
experiment_setup = 'flat_gaussian'; % CHANGE THIS TO SELECT SETUP
config.experiment_setup = experiment_setup; % Store the chosen setup name in config

fprintf('Selected experiment setup: %s\n', experiment_setup);

% ======================================================================
% --- Common Settings (can be overridden by specific setups below) ---
% ======================================================================

% --- Domain and Mesh ---
config.domain.xmin = 0.0;
config.domain.xmax = 20.0; % Default domain length [m]
config.mesh.N      = 500;  % Default number of cells
config.param.H0    = 0.50; % Default undisturbed water depth [m]

% --- Model and Numerics ---
config.model = @core.rhs_nsw_1st_order;        % RHS function (1st order FV)
config.numFlux = @flux.FVCF;                   % Numerical flux
config.reconstructopenion = [];                    % No reconstruction (1st order)
config.timeStepper = @time.integrate_euler_adaptive; % Time integration
config.time.CFL = 0.95;                        % CFL number

% --- Run Control ---
config.t0 = 0.0;
config.tEnd = 10.0;            % Default end time [s]
config.vis.dt_plot = 0.1; % Output interval for visualization and saving [s]
config.vis.plot_velocity = true; % Set to true to plot velocity in a subpanel
config.vis.show_legend = false; % Set to true to show legend in wave tank plot
config.tspan = config.t0:config.vis.dt_plot:config.tEnd;
config.save_results = false; % Save results by default
config.output_path = 'results/'; % Default output directory

% ======================================================================
% --- Specific Experiment Setups ---
% ======================================================================

switch experiment_setup
    case 'flat_rest'
        % Simple case: flat bottom, lake at rest, wall boundaries.
        % Useful for testing stability.
        config.caseName = 'flat_rest_L20m_H0.5m_N400';

        % Bathymetry
        config.bathyHandle = @cfg.bathy.flat;

        % Initial Condition
        config.ic_handle = @ic.lake_at_rest;

        % Boundary Conditions
        config.bc.left.handle = @bc.wall;
        config.bc.right.handle = @bc.wall;

    case 'flat_gaussian'
        % Flat bottom, Gaussian bump IC, wall boundaries.
        config.caseName = 'flat_gaussian_L20m_H0.5m_N400';
        config.tEnd = 10.0; % Longer time for bump propagation
        config.tspan = linspace(config.t0, config.tEnd, 101);

        % Bathymetry
        config.bathyHandle = @cfg.bathy.flat;

        % Initial Condition
        config.ic_handle = @ic.gaussian_bump;
        config.ic_param.a = 0.1;       % Amplitude
        config.ic_param.lambda = 0.5;  % Decay rate
        config.ic_param.x0 = config.domain.xmax / 2; % Center position

        % Boundary Conditions
        config.bc.left.handle = @bc.wall;
        % config.bc.right.handle = @bc.wall;
        % Alternatively, use open boundaries:
        % config.bc.left.handle = @bc.open;
        config.bc.right.handle = @bc.open;

    case 'flat_wave_gen'
        % Flat bottom, wave generation at left, wall at right.
        config.caseName = 'flat_wave_gen_L20m_H0.5m_N400';

        % Bathymetry
        config.bathyHandle = @cfg.bathy.flat;

        % Initial Condition
        config.ic_handle = @ic.lake_at_rest;
        config.ic_param = struct();

        % Boundary Conditions
        % Left: Generating
        config.bc.left.handle = @bc.generating;
        config.bc.left.param.a = 0.05; % Wave amplitude (m)
        config.bc.left.param.T = 2.0;  % Wave period (s)

        % Right: Wall
        config.bc.right.handle = @bc.wall;

    otherwise
        error('Unknown experiment_setup: %s', experiment_setup);
    end

    % ======================================================================
    % --- Final Setup Steps ---
    % ======================================================================

    % Generate mesh based on final domain and N
    [config.mesh.xc, config.mesh.dx, config.mesh.x_edge] = core.utils.uniform(config.domain, config.mesh.N);

    % Define output path based on case name
    config.outputPath = fullfile('./results/', config.caseName);
    if ~isfolder(config.outputPath)
        mkdir(config.outputPath);
        fprintf('Created results directory: %s\n', config.outputPath);
    end

    fprintf('Configuration loaded for case: %s\n', config.caseName);
    fprintf('  Domain: [%.2f, %.2f] m, Cells: %d\n', config.domain.xmin, config.domain.xmax, config.mesh.N);
    fprintf('  Bathymetry: %s\n', func2str(config.bathyHandle));
    fprintf('  Initial Condition: %s\n', func2str(config.ic_handle));
    if isfield(config, 'ic_param') && ~isempty(fieldnames(config.ic_param))
        fprintf('    IC Params: %s\n', core.utils.struct2str(config.ic_param));
    end
    fprintf('  BC Left: %s', func2str(config.bc.left.handle));
    if isfield(config.bc.left, 'param') && ~isempty(fieldnames(config.bc.left.param))
        fprintf(' (Params: %s)', core.utils.struct2str(config.bc.left.param)); % Requires a helper struct2str or similar
    end
    fprintf('\n');
    fprintf('  BC Right: %s', func2str(config.bc.right.handle));
     if isfield(config.bc.right, 'param') && ~isempty(fieldnames(config.bc.right.param))
         fprintf(' (Params: %s)', core.utils.struct2str(config.bc.right.param)); % Requires a helper struct2str or similar
     end
    fprintf('\n');
    fprintf('  Time Span: [%.2f, %.2f] s, CFL: %.2f\n', config.t0, config.tEnd, config.time.CFL);

end