function [sol_out, t_out, stats] = integrate_matlab_ode(rhs_func, tspan, w0, cfg)

    % INTEGRATE_MATLAB_ODE Solves ODEs using a built-in MATLAB ODE solver.
    %   Acts as a wrapper for standard MATLAB ODE solvers like ode45, ode113, ode23, etc.
    %
    %   [SOL_OUT, T_OUT, STATS] = INTEGRATE_MATLAB_ODE(RHS_FUNC, TSPAN, W0, CFG)
    %   integrates the system defined by RHS_FUNC over TSPAN with initial
    %   condition W0.
    %
    %   Inputs:
    %   RHS_FUNC   - Function handle for the right-hand side, expected to have the
    %                signature rhs_func(t, w, cfg).
    %   TSPAN      - Vector specifying the time points for output [t0, t1, ..., tf].
    %   W0         - Initial condition vector.
    %   CFG        - Configuration structure. Should contain:
    %                - CFG.TIME.MATLAB_SOLVER (e.g., 'ode45', 'ode113', 'ode23')
    %                - CFG.TIME.AbsTol (Absolute tolerance, default 1e-4)
    %                - CFG.TIME.RelTol (Relative tolerance, default 1e-4)
    %                - CFG.TIME.ODE_OPTIONS (Optional, structure from odeset)
    %
    %   Outputs:
    %   SOL_OUT    - Matrix of solution vectors at times in T_OUT. Each row
    %                corresponds to a time point (transposed from solver output).
    %   T_OUT      - Row vector of time points corresponding to the solution points.
    %   STATS      - Basic statistics structure (contains nsteps).

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
    base_options = odeset(); % Start with default options
    if isfield(cfg, 'time') && isfield(cfg.time, 'ode_options')
        base_options = cfg.time.ode_options; % Use user-provided base options if available
        fprintf('Using base ODE options specified in cfg.time.ode_options.\n');
    end

    % Set Tolerances (use defaults if not specified)
    abs_tol = 1e-4; % Default absolute tolerance
    if isfield(cfg, 'time') && isfield(cfg.time, 'AbsTol')
        abs_tol = cfg.time.AbsTol;
        fprintf('Using specified AbsTol = %g.\n', abs_tol);
    else
        fprintf('Using default AbsTol = %.4e.\n', abs_tol);
    end

    rel_tol = 1e-4; % Default relative tolerance
    if isfield(cfg, 'time') && isfield(cfg.time, 'RelTol')
        rel_tol = cfg.time.RelTol;
        fprintf('Using specified RelTol = %g.\n', rel_tol);
    else
        fprintf('Using default RelTol = %.4e.\n', rel_tol);
    end

    % Combine base options with specific tolerances
    options = odeset(base_options, 'AbsTol', abs_tol, 'RelTol', rel_tol);

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
    % rhs_func already has the correct signature from core.solver
    ode_rhs = rhs_func;

    % --- Call Solver ---
    fprintf('Starting integration with MATLAB solver: %s from t=%.3f to t=%.3f\n', ...
            solver_name, tspan(1), tspan(end));

    [t_solver, sol_solver] = solver_handle(ode_rhs, tspan, w0, options);

    % --- Format Output ---
    % MATLAB solvers return time as a column vector and solutions where rows
    % correspond to time points (M_out x 2*N). core.solver expects this format.
    t_out = t_solver'; % Transpose time to be a row vector as expected by our setup
    sol_out = sol_solver; % NO transpose needed for solution matrix

    % Initialize stats structure
    stats = struct('nsteps', NaN, 'nfevals', NaN);

    % Step count and dt history are not directly comparable/available
    stats.nsteps = length(t_out) - 1; % Approximate number of steps

    fprintf('Integration with %s finished at t = %.3f s. Output generated at %d time points.\n', ...
            solver_name, t_out(end), length(t_out));

end % Function end