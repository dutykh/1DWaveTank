function results = solver(cfg)
% SOLVER Main simulation driver for the 1DWaveTank.
%
%   This function orchestrates the simulation process:
%   1. Sets up the initial conditions based on the configuration.
%   2. Selects the appropriate Right-Hand Side (RHS) function for the chosen model.
%   3. Calls the specified time integration scheme.
%   4. Processes and stores the results (time vector, H, HU).
%
%   Inputs:
%     cfg - Configuration structure containing all simulation parameters:
%           cfg.mesh: Mesh details (N, x, xc, dx).
%           cfg.time: Time integration parameters (t_span, integrator handle, dt_plot, cfl).
%           cfg.phys: Physical parameters (g).
%           cfg.prob: Problem-specific parameters (initial conditions handle, etc.).
%           cfg.numerics: Numerical scheme details (rhs handle, flux handle, etc.).
%           cfg.bc: Boundary condition handles.
%
%   Outputs:
%     results - Structure containing the simulation output:
%               results.t: Column vector of output time points.
%               results.H: Matrix of water depth H at cell centers (M_out x N).
%               results.HU: Matrix of discharge HU at cell centers (M_out x N).
%               results.xc: Vector of cell center coordinates.
%               results.cfg: The configuration structure used for the run.
%               results.total_steps: Total number of time steps taken.
%               results.dt_history: Vector of dt values used at each step.

    fprintf('--- Starting Core Solver ---\n');
    
    % --- Input Validation (Basic) ---
    required_fields = {'tspan', 'model', 'timeStepper', 'mesh', 'ic_handle', 'ic_param'};
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
    w_init = cfg.ic_handle(cfg.mesh.xc, cfg.ic_param);
    N = cfg.mesh.N;
    if isvector(w_init) && length(w_init) == 2*N
        % Already in [H; HU] format
        w0 = w_init(:);
    elseif isvector(w_init) && length(w_init) == N
        % Only H is provided, assume HU = 0
        w0 = [w_init(:); zeros(N,1)];
    elseif size(w_init,2) >= 2 && size(w_init,1) == N
        % Use only the first two columns (H and HU provided, ignore extras)
        w0 = [w_init(:,1); w_init(:,2)];
    else
        error('Initial condition function must return N x 1, N x 2 (or more), or 2*N x 1 array.');
    end
    
    % Get time span for output
    tspan = cfg.tspan;
    if isempty(tspan) || length(tspan) < 2
        error('cfg.tspan must contain at least start and end times.');
    end
    
    % Get function handles from config
    rhs_handle = cfg.model; % Handle to the RHS function (e.g., @core.rhs_nsw_1st_order)
    time_stepper = cfg.timeStepper; % Handle to the time integration function
    
    if isempty(rhs_handle) || ~isa(rhs_handle, 'function_handle')
        error('cfg.model must be a valid function handle to the RHS function.');
    end
    if isempty(time_stepper) || ~isa(time_stepper, 'function_handle')
        error('cfg.timeStepper must be a valid function handle to the time integration function.');
    end
    
    % --- Run Time Integration ---
    fprintf('Calling time integrator...\n');
    % The time stepper function (e.g., integrate_euler_adaptive) performs the actual loop
    [t_out, sol_out, total_steps, dt_history] = time_stepper(rhs_handle, tspan, w0, cfg);
    fprintf('Time integration completed.\n');
    
    % --- Process and Store Results ---
    fprintf('Processing results...\n');
    N = cfg.mesh.N; % Number of spatial cells
    M_out = length(t_out); % Number of output time steps
    
    % Reshape the flat solution vector returned by the time stepper
    % sol_out is expected to be M_out x (2*N)
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
    
    % Store total steps and dt history
    results.total_steps = total_steps;
    results.dt_history = dt_history;
    
    fprintf('--- Core Solver Finished ---\n');
 
end