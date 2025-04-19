function cfg = flat_bottom_config()

    %FLAT_BOTTOM_CONFIG Configuration for a 5m flat-bottom tank.
    %   cfg = FLAT_BOTTOM_CONFIG() returns a fully-initialised configuration
    %   structure for a 5-meter long tank with constant depth H0 = 0.5m,
    %   discretized with N = 500 cells.
    
    cfg = cfg.default_config();    % Start from safe defaults
    
    % --- Override Defaults for this Specific Case ---
    
    % Domain and spatial mesh
    cfg.domain.xmin = 0.0;         % Left wall [m]
    cfg.domain.xmax = 5.0;         % Right boundary [m]
    cfg.mesh.N      = 500;         % Number of control volumes
    [cfg.mesh.xc, cfg.mesh.dx, cfg.mesh.x_edge] = mesh.uniform(cfg.domain, cfg.mesh.N); % Generate mesh
    
    % Bathymetry and initial condition
    cfg.param.H0    = 0.50;        % Undisturbed water depth [m]
    cfg.bathyHandle = @cfg.bathy.flat; % Use the flat bathymetry function
    cfg.ic          = @ic.lake_at_rest; % Use lake at rest IC
    
    % Governing model and numerical options (Using placeholders for now)
    cfg.model = []; % Placeholder: Set wrapper in run_simulation.m
    cfg.numFlux = @flux.placeholder_flux;   % Placeholder numerical flux
    cfg.reconstruction = [];       % No reconstruction (1st order method)
    cfg.timeStepper = @time.placeholder_stepper; % Placeholder time stepper
    
    % Boundary conditions (Using placeholders for now)
    cfg.bcL = @bc.placeholder_bc;  % Placeholder BC at x = 0
    cfg.bcR = @bc.placeholder_bc;  % Placeholder BC at x = L
    
    % Run-control parameters
    cfg.tEnd = 10.0;               % Simulated time span [s]
    cfg.outputEvery = 100;         % Snapshot frequency (relative to integrator steps)
    
    % House-keeping
    cfg.caseName = 'flat_bottom_L5m_H0.5m_N500'; % Descriptive name
    cfg.outputPath = fullfile('./results/', cfg.caseName); % Output path
    
end