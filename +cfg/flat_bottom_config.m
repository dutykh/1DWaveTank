function cfg = flat_bottom_config()

    %FLAT_BOTTOM_CONFIG Configuration for a 5m flat-bottom tank (1st Order Euler).
    %   cfg = FLAT_BOTTOM_CONFIG() returns a fully-initialised configuration
    %   structure for a 5-meter long tank with constant depth H0 = 0.5m,
    %   discretized with N = 500 cells. Uses 1st order FVCF flux and
    %   adaptive explicit Euler time stepping.
    
    cfg = default_config(); % Start from safe defaults
    
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
    
    % Governing model and numerical options
    % The RHS function handle itself (will be called by the time stepper)
    cfg.model = @core.rhs_nsw_1st_order;
    
    % Select the numerical flux and time stepper
    cfg.numFlux = @flux.FVCF;              % Use the FVCF numerical flux
    cfg.reconstruction = [];               % No reconstruction (1st order method)
    cfg.timeStepper = @time.integrate_euler_adaptive; % Use adaptive Euler
    cfg.time.CFL = 0.9;                    % Set the CFL number for the adaptive stepper
    
    % Boundary conditions (Using placeholders for now - NEED TO BE IMPLEMENTED)
    % Replace these with actual boundary condition functions (e.g., @bc.wall)
    % IMPORTANT: The placeholder BC needs to be updated to handle ghost cells
    cfg.bcL = @bc.placeholder_bc_1st_order;  % Placeholder BC function for left
    cfg.bcR = @bc.placeholder_bc_1st_order;  % Placeholder BC function for right
    
    % Run-control parameters
    cfg.tEnd = 10.0;               % Simulated time span [s]
    
    % Define output times explicitly using tspan
    num_output_points = 101; % Example: 101 points including t0 and tEnd
    cfg.tspan = linspace(cfg.t0, cfg.tEnd, num_output_points);
    
    % House-keeping
    cfg.caseName = 'flat_bottom_L5m_H0.5m_N500_1stOrderEuler'; % Descriptive name
    cfg.outputPath = fullfile('./results/', cfg.caseName); % Output path
    
    % Create the output directory if it doesn't exist
    if ~isfolder(cfg.outputPath)
        mkdir(cfg.outputPath);
        fprintf('Created results directory: %s\n', cfg.outputPath);
    end

end