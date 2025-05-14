%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +core/solver.m
%
% Purpose:
%   Main simulation driver for the 1DWaveTank code. Orchestrates the
%   simulation process: sets up initial conditions, selects the appropriate
%   RHS function and time integrator, runs the simulation, and processes output.
%
% Syntax:
%   results = solver(cfg)
%
% Inputs:
%   cfg - [struct] Configuration structure containing all simulation parameters:
%         cfg.mesh: Mesh details (N, x, xc, dx)
%         cfg.time: Time integration parameters (t_span, integrator handle, dt_plot, cfl)
%         cfg.phys: Physical parameters (g)
%         cfg.prob: Problem-specific parameters (initial conditions handle, etc.)
%         cfg.numerics: Numerical scheme details (rhs handle, flux handle, etc.)
%         cfg.bc: Boundary condition handles
%
% Outputs:
%   results - [struct] Structure containing the simulation output:
%             results.t:    [M_out x 1] Column vector of output time points
%             results.H:    [M_out x N] Matrix of water depth H at cell centers
%             results.HU:   [M_out x N] Matrix of discharge HU at cell centers
%             results.U:    [M_out x N] Matrix of velocity U at cell centers
%             results.xc:   [N x 1] Vector of cell center coordinates
%             results.cfg:  [struct] The configuration structure used for the run
%             results.total_steps: [integer] Total number of time steps taken
%
% Dependencies:
%   Expects properly configured cfg structure and function handles for IC, RHS, integrator, etc.
%
% References:
%   - See 1DWaveTank UserGuide.md for structure and usage.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = solver(cfg)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Validation (Basic)                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('--- Starting Core Solver ---\n');
    required_fields = {'tspan', 'model', 'timeStepper', 'mesh', 'ic_handle', 'ic_param'};
    for i = 1:length(required_fields)
        if ~isfield(cfg, required_fields{i})
            error('Configuration structure `cfg` is missing required field: %s', required_fields{i});
        end
    end
    if ~isfield(cfg.mesh, 'N')
        error('Configuration structure `cfg.mesh` is missing required field: N');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial Condition Setup                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ensure mesh center coordinates are a row vector (1xN)
    if isfield(cfg.mesh, 'xc')
        if iscolumn(cfg.mesh.xc)
            cfg.mesh.xc = cfg.mesh.xc.'; % Transpose to row vector
            fprintf('Transposed cfg.mesh.xc to row vector for compatibility.\n');
        end
    end
    fprintf('Setting up initial condition...\n');
    % The initial condition function should return the state vector [H; HU]
    w_init = cfg.ic_handle(cfg); % Pass full config to IC handle (refactored for lake_at_rest)
    
    % --- Reshape Initial Condition --- 
    N = cfg.mesh.N;
    if isvector(w_init) && length(w_init) == 2*N
        % Already in [H; HU] format (flattened)
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time Span Preparation                                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tspan = cfg.tspan;
    if isempty(tspan) || length(tspan) < 2
        error('cfg.tspan must contain at least start and end times.');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare Function Handles                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine the function handle for the time stepper
    time_stepper = cfg.timeStepper;
    
    % Handle for the RHS function (differs based on solver type)
    if isequal(time_stepper, @time.integrate_matlab_ode)
        % MATLAB solvers need f(t,w), so wrap cfg.model to include cfg
        rhs_handle = @(t, w) cfg.model(t, w, cfg);
    else
        % Our custom solvers will receive cfg separately and call f(t,w,cfg)
        rhs_handle = cfg.model; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time Integration                                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic; % Start timer for integration
    % Call the selected time stepper with the appropriate RHS handle
    [sol_out, t_out, stats] = time_stepper(rhs_handle, tspan, w0, cfg);
    integration_time = toc; % Stop timer
    total_steps = stats.nsteps;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Post-processing & Output                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Time integration completed in %.2f seconds.\n', integration_time);
    fprintf('Total time steps taken: %d\n', total_steps);
    
    % Check output dimensions (sol_out should have time points as rows)
    if length(t_out) ~= size(sol_out, 1)
        error('Dimension mismatch: length(t_out)=%d vs size(sol_out,1)=%d.', length(t_out), size(sol_out, 1));
    end
    
    % Prepare results structure
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
    
    % Store total steps and dt history (if available)
    results.total_steps = total_steps;
    % results.dt_history = stats.dt_history; % dt_history is no longer reliably available
    
    fprintf('--- Core Solver Finished ---\n');
 
end