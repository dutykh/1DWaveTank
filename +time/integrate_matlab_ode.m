function [t_out, sol_out, k, dt_history] = integrate_matlab_ode(rhs_func, tspan, w0, cfg)

% INTEGRATE_MATLAB_ODE Solves ODEs using a built-in MATLAB ODE solver.
%   Acts as a wrapper for standard MATLAB ODE solvers like ode45, ode113, etc.
%
%   [T_OUT, SOL_OUT, K, DT_HISTORY] = INTEGRATE_MATLAB_ODE(RHS_FUNC, TSPAN, W0, CFG)
%   integrates the system defined by RHS_FUNC over TSPAN with initial
%   condition W0.
%
%   Inputs:
%   RHS_FUNC   - Function handle for the right-hand side, expected to have the
%                signature rhs_func(t, w, cfg).
%   TSPAN      - Vector specifying the time points for output [t0, t1, ..., tf].
%   W0         - Initial condition vector.
%   CFG        - Configuration structure. Should contain CFG.TIME.MATLAB_SOLVER
%                (e.g., 'ode45', 'ode113') and optionally CFG.TIME.ODE_OPTIONS
%                (created using odeset).
%
%   Outputs:
%   T_OUT      - Row vector of time points corresponding to the solution points.
%   SOL_OUT    - Matrix of solution vectors at times in T_OUT. Each row
%                corresponds to a time point (transposed from solver output).
%   K          - Returns NaN (step count not directly comparable).
%   DT_HISTORY - Returns NaN (dt history not directly available).

    % --- Configuration ---
    % Get solver name, default to ode113 if not specified
    if isfield(cfg, 'time') && isfield(cfg.time, 'matlab_solver')
        solver_name = cfg.time.matlab_solver;
    else
        solver_name = 'ode113'; % Default solver
        warning('cfg.time.matlab_solver not specified, defaulting to ''ode113''.');
    end
    solver_handle = str2func(solver_name); % Get function handle from name

    % Get ODE options, use defaults if not specified
    if isfield(cfg, 'time') && isfield(cfg.time, 'ode_options')
        options = cfg.time.ode_options;
        fprintf('Using custom ODE options specified in cfg.time.ode_options.\n');
    else
        options = odeset(); % Use default options
        fprintf('Using default ODE options for %s.\n', solver_name);
    end

    % Check if progress bar should be shown (default: true)
    show_progress = true; % Default value
    if isfield(cfg, 'time') && isfield(cfg.time, 'show_progress_bar')
        show_progress = cfg.time.show_progress_bar;
    end

    % Add odetpbar to OutputFcn if requested
    if show_progress
        existing_outputfcn = odeset(options).OutputFcn;
        if isempty(existing_outputfcn)
            options = odeset(options, 'OutputFcn', @core.utils.odetpbar);
        elseif isa(existing_outputfcn, 'function_handle')
            % If single handle, create cell array
            options = odeset(options, 'OutputFcn', {existing_outputfcn, @core.utils.odetpbar});
        elseif iscell(existing_outputfcn)
            % If already cell array, append
            options = odeset(options, 'OutputFcn', [existing_outputfcn, {@core.utils.odetpbar}]);
        end
        fprintf('Adding progress bar (@core.utils.odetpbar) to OutputFcn.\n');
    else
        fprintf('Progress bar disabled (cfg.time.show_progress_bar = false).\n');
    end

    % Adapt the RHS function handle to the f(t, y) signature expected by solvers
    % Capture cfg in the anonymous function context
    ode_rhs = @(t, w) feval(rhs_func, t, w, cfg);

    % --- Call Solver ---
    fprintf('Starting integration with MATLAB solver: %s from t=%.3f to t=%.3f\n', ...
            solver_name, tspan(1), tspan(end));

    [t_solver, sol_solver] = solver_handle(ode_rhs, tspan, w0, options);

    % --- Format Output ---
    % MATLAB solvers return time as a column vector and solutions where rows
    % correspond to time points (M_out x 2*N). core.solver expects this format.
    t_out = t_solver'; % Transpose time to be a row vector as expected by our setup
    sol_out = sol_solver; % NO transpose needed for solution matrix

    % Step count and dt history are not directly comparable/available
    k = NaN;
    dt_history = NaN;

    fprintf('Integration with %s finished at t = %.3f s. Output generated at %d time points.\n', ...
            solver_name, t_out(end), length(t_out));

end % Function end
