function config = simulation_config()
% config.simulation_config: Central configuration file for 1DWaveTank simulations.
%   config = config.simulation_config() returns a structure 'config' containing all
%   parameters and function handles for a specific simulation setup.
%   Modify this file to set up different test cases.

config = cfg.default_config(); % Start from safe defaults

% --- Experiment Selection ---
% Choose a predefined setup or define a custom one below
% Available setups: 'flat_rest', 'flat_gaussian', 'flat_wave_gen'
experiment_setup = 'flat_wave_gen'; % CHANGE THIS TO SELECT SETUP

fprintf('Setting up simulation configuration: %s\n', experiment_setup);

% ======================================================================
% --- Common Settings (can be overridden by specific setups below) ---
% ======================================================================

% --- Domain and Mesh ---
    config.domain.xmin = 0.0;
    config.domain.xmax = 20.0; % Default domain length [m]
    config.mesh.N      = 400;  % Default number of cells
    config.param.H0    = 0.50; % Default undisturbed water depth [m]

    % --- Model and Numerics ---
    config.model = @core.rhs_nsw_1st_order;        % RHS function (1st order FV)
    config.numFlux = @flux.FVCF;                   % Numerical flux
    config.reconstruction = [];                    % No reconstruction (1st order)
    config.timeStepper = @time.integrate_euler_adaptive; % Time integration
    config.time.CFL = 0.8;                         % CFL number

    % --- Run Control ---
    config.t0 = 0.0;
    config.tEnd = 15.0;            % Default end time [s]
    num_output_points = 151;    % Default number of output snapshots
    config.tspan = linspace(config.t0, config.tEnd, num_output_points);

    % ======================================================================
    % --- Specific Experiment Setups ---
    % ======================================================================

    switch experiment_setup
        case 'flat_rest'
            % --- Flat bottom, lake at rest, wall boundaries ---
            config.caseName = 'flat_rest_L20m_H0.5m_N400';

            % Bathymetry
            config.bathyHandle = @cfg.bathy.flat;

            % Initial Condition
            config.ic = @ic.lake_at_rest;

            % Boundary Conditions
            config.bc.left.handle = @bc.wall;
            config.bc.right.handle = @bc.wall;

        case 'flat_gaussian'
             % --- Flat bottom, Gaussian bump IC, wall boundaries ---
            config.caseName = 'flat_gaussian_L20m_H0.5m_N400';
            config.tEnd = 20.0; % Longer time for bump propagation
            config.tspan = linspace(config.t0, config.tEnd, 201);

            % Bathymetry
            config.bathyHandle = @config.bathy.flat;

            % Initial Condition
            config.ic = @ic.gaussian_bump;
            % Parameters for Gaussian bump (can be omitted to use defaults)
            config.ic.param.a = 0.1;       % Amplitude
            config.ic.param.lambda = 0.5;  % Decay rate
            config.ic.param.x0 = config.domain.xmax / 4; % Center position

            % Boundary Conditions
            config.bc.left.handle = @bc.wall;
            config.bc.right.handle = @bc.wall;
            % Alternatively, use open boundaries:
            % config.bc.left.handle = @bc.open;
            % config.bc.right.handle = @bc.open;

        case 'flat_wave_gen'
            % --- Flat bottom, wave generation at left, wall at right ---
            config.caseName = 'flat_wave_gen_L20m_H0.5m_N400';

            % Bathymetry
            config.bathyHandle = @cfg.bathy.flat;

            % Initial Condition
            config.ic = @ic.lake_at_rest;

            % Boundary Conditions
            % Left: Generating
            config.bc.left.handle = @bc.generating;
            config.bc.left.param.a = 0.05; % Wave amplitude (m)
            config.bc.left.param.T = 2.0;  % Wave period (s)

            % Right: Wall
            config.bc.right.handle = @bc.wall;
            % Alternatively, use open boundary:
            % config.bc.right.handle = @bc.open;

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
    fprintf('  Initial Condition: %s\n', func2str(config.ic));
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