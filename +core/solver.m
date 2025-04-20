function results = solver(cfg)

    %SOLVER Main simulation driver for the 1DWaveTank.
    %   results = SOLVER(cfg) runs the simulation based on the settings
    %   provided in the configuration structure 'cfg'.
    %
    %   Inputs:
    %       cfg - Structure containing all configuration parameters and function
    %             handles for the simulation (e.g., initial condition, RHS model,
    %             time stepper, numerical flux, boundary conditions, mesh, etc.).
    %
    %   Outputs:
    %       results - Structure containing the simulation results, including:
    %                 .t: Vector of time points where the solution is saved.
    %                 .H: Array (time x space) of water depth H.
    %                 .HU: Array (time x space) of discharge HU.
    %                 .U: Array (time x space) of velocity U (calculated).
    %                 .cfg: The configuration structure used for this run.
    
    fprintf('--- Starting Core Solver ---\n');
    
    % --- Input Validation (Basic) ---
    required_fields = {'ic', 'tspan', 'model', 'timeStepper', 'mesh'};
    for i = 1:length(required_fields)
        if ~isfield(cfg, required_fields{i})
            error('Configuration structure `cfg` is missing required field: %s', required_fields{i});
        end
    end
    if ~isfield(cfg.mesh, 'N')
        error('Configuration structure `cfg.mesh` is missing required field: N');
    end
    
    % --- Setup ---
    fprintf('Setting up initial condition...\n');
    w0 = cfg.ic(cfg); % Get initial state vector [H0; HU0]
    
    % Get time span for output
    tspan = cfg.tspan;
    if isempty(tspan) || length(tspan) < 2
        error('cfg.tspan must contain at least start and end times.');
    end
    fprintf('Time span: [%.2f s, %.2f s]\n', tspan(1), tspan(end));
    
    % Get function handles from config
    rhs_handle = cfg.model; % Handle to the RHS function (e.g., @core.rhs_nsw_1st_order)
    time_stepper = cfg.timeStepper; % Handle to the time integration function
    
    if isempty(rhs_handle) || ~isa(rhs_handle, 'function_handle')
        error('cfg.model must be a valid function handle to the RHS function.');
    end
    if isempty(time_stepper) || ~isa(time_stepper, 'function_handle')
        error('cfg.timeStepper must be a valid function handle to the time integration function.');
    end
    
    fprintf('Using time stepper: %s\n', func2str(time_stepper));
    fprintf('Using RHS model: %s\n', func2str(rhs_handle));
    
    % --- Run Time Integration ---
    fprintf('Calling time integrator...\n');
    % The time stepper function (e.g., integrate_euler_adaptive) performs the actual loop
    [t_out, sol_out] = time_stepper(rhs_handle, tspan, w0, cfg);
    fprintf('Time integration completed.\n');
    
    % --- Process and Store Results ---
    fprintf('Processing results...\n');
    N = cfg.mesh.N; % Number of spatial cells
    M_out = length(t_out); % Number of output time steps
    
    % Reshape the flat solution vector returned by the time stepper
    % sol_out is expected to be M_out x (2*N)
    fprintf('[DEBUG] core.solver: Checking sol_out dimensions. Expected cols = %d, Received size = [%d, %d]\n', 2*N, size(sol_out, 1), size(sol_out, 2)); % DEBUG
    if size(sol_out, 2) ~= 2*N
        error('Time stepper returned solution array with unexpected dimensions.');
    end
    
    results = struct();
    results.t = t_out(:); % Ensure time is a column vector
    
    % Extract H and HU (M_out x N)
    results.H = sol_out(:, 1:N);
    results.HU = sol_out(:, N+1:2*N);
    
    % Calculate velocity U (handle division by zero)
    results.U = zeros(M_out, N);
    wet_mask = results.H > 1e-10; % Avoid division by zero in dry cells
    results.U(wet_mask) = results.HU(wet_mask) ./ results.H(wet_mask);
    
    % Store the configuration used for this run
    results.cfg = cfg;
    
    fprintf('--- Core Solver Finished ---\n');
 
end