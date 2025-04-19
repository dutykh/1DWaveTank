function out = flat_bottom_config()

    %FLAT_BOTTOM_CONFIG Configuration for a 5m flat-bottom tank (1st Order Euler).
    %   cfg = FLAT_BOTTOM_CONFIG() returns a fully-initialised configuration
    %   structure for a 5-meter long tank with constant depth H0 = 0.5m,
    %   discretized with N = 500 cells. Uses 1st order FVCF flux and
    %   adaptive explicit Euler time stepping.
    
    out = cfg.default_config(); % Start from safe defaults
    
    % --- Override Defaults for this Specific Case ---
    
    % Domain and spatial mesh
    out.domain.xmin = 0.0;         % Left wall [m]
    out.domain.xmax = 5.0;         % Right boundary [m]
    out.mesh.N      = 500;         % Number of control volumes
    [out.mesh.xc, out.mesh.dx, out.mesh.x_edge] = core.utils.uniform(out.domain, out.mesh.N); % Generate mesh
    
    % Bathymetry and initial condition
    out.param.H0    = 0.50;        % Undisturbed water depth [m]
    out.bathyHandle = @cfg.bathy.flat; % Use the flat bathymetry function
    out.ic          = @ic.lake_at_rest; % Use lake at rest IC
    
    % Governing model and numerical options
    % The RHS function handle itself (will be called by the time stepper)
    out.model = @core.rhs_nsw_1st_order;
    
    % Select the numerical flux and time stepper
    out.numFlux = @flux.FVCF;              % Use the FVCF numerical flux
    out.reconstruction = [];               % No reconstruction (1st order method)
    out.timeStepper = @time.integrate_euler_adaptive; % Use adaptive Euler
    out.time.CFL = 0.9;                    % Set the CFL number for the adaptive stepper
    
    % Boundary conditions (Using placeholders for now - NEED TO BE IMPLEMENTED)
    % Replace these with actual boundary condition functions (e.g., @bc.wall)
    % IMPORTANT: The placeholder BC needs to be updated to handle ghost cells
    out.bcL = @bc.placeholder_bc_1st_order;  % Placeholder BC function for left
    out.bcR = @bc.placeholder_bc_1st_order;  % Placeholder BC function for right
    
    % Run-control parameters
    out.tEnd = 10.0;               % Simulated time span [s]
    
    % Define output times explicitly using tspan
    num_output_points = 101; % Example: 101 points including t0 and tEnd
    out.tspan = linspace(out.t0, out.tEnd, num_output_points);
    
    % House-keeping
    out.caseName = 'flat_bottom_L5m_H0.5m_N500_1stOrderEuler'; % Descriptive name
    out.outputPath = fullfile('./results/', out.caseName); % Output path
    
    % Create the output directory if it doesn't exist
    if ~isfolder(out.outputPath)
        mkdir(out.outputPath);
        fprintf('Created results directory: %s\n', out.outputPath);
    end

end