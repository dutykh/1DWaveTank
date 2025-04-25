%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +time/integrate_matlab_ode.m
%
% Purpose:
%   Provides a standardized wrapper function to integrate a system of Ordinary
%   Differential Equations (ODEs) using one of MATLAB's built-in, high-quality
%   ODE solvers (e.g., ode45, ode113, ode15s). This approach leverages
%   MATLAB's sophisticated adaptive time-stepping algorithms and robust error
%   control mechanisms, which can be beneficial for complex or stiff problems,
%   while still fitting into the 1DWaveTank's modular structure.
%
% Syntax:
%   [sol_out, t_out, stats] = integrate_matlab_ode(rhs_func, tspan, w0, cfg)
%
% Inputs:
%   rhs_func - [function handle] Handle to the function defining the RHS of
%                the ODE system: dw/dt = f(t, w). IMPORTANT: The `core.solver`
%                passes a wrapped function handle here. The original RHS functions
%                (like `+core/rhs_nsw_1st_order.m`) expect the signature f(t, w, cfg),
%                but the wrapper adapts it to the f(t, w) signature required by
%                MATLAB's standard ODE solvers.
%   tspan    - [vector, double] Time points [t0, t1, ..., tf] at which the
%                solution output is requested. The solver integrates from
%                tspan(1) to tspan(end), potentially taking many internal steps
%                but only returning the solution at these specified times.
%   w0       - [vector, double] Initial state vector (column vector) at time tspan(1).
%   cfg      - [struct] Configuration structure. Relevant fields examined here:
%                cfg.time.matlab_solver: [char] Name of the MATLAB ODE solver
%                                       to use (e.g., 'ode113', 'ode45'). Default: 'ode113'.
%                                       'ode113' is often a good default for non-stiff
%                                       problems requiring moderate accuracy.
%                cfg.time.AbsTol: [double] Absolute tolerance for the solver's
%                                error control. Controls the threshold below which
%                                the solution component value is considered zero.
%                                Default: 1e-4.
%                cfg.time.RelTol: [double] Relative tolerance for the solver's
%                                error control. Controls the number of correct
%                                digits relative to the magnitude of the solution.
%                                Default: 1e-4.
%                cfg.time.ode_options: [struct] Optional base `odeset` structure
%                                     created using `odeset`. Tolerances and OutputFcn
%                                     specified here will override the base settings.
%                                     Default: odeset().
%                cfg.time.show_progress_bar: [logical] If true, adds a text
%                                          progress bar to the command window via
%                                          the OutputFcn mechanism. Default: true.
%                (Other cfg fields are passed implicitly to the wrapped rhs_func).
%
% Outputs:
%   sol_out  - [M x length(w0), double] Solution matrix. Each row `sol_out(i,:)`
%                is the state vector corresponding to the time point `t_out(i)`.
%   t_out    - [1 x M, double] Row vector of time points where the solution is output.
%                These are exactly the times specified in the input `tspan` vector.
%   stats    - [struct] Basic statistics structure:
%                stats.nsteps: Approximated as the number of output intervals.
%                              This does *not* reflect the number of internal steps
%                              taken by the adaptive solver.
%                stats.nfevals: Number of RHS evaluations. This information
%                               is not directly available from standard MATLAB odeXX
%                               solver outputs without using a more complex OutputFcn
%                               or examining the `sol.stats` structure (if using the
%                               `sol = odeXX(...)` syntax). Set to NaN for simplicity.
%
% Dependencies:
%   - Requires the specified MATLAB ODE solver function (e.g., ode113) to be available.
%   - Requires +utils/odetpbar.m if cfg.time.show_progress_bar is true.
%
% References:
%   - MATLAB documentation for `ode45`, `ode113`, `odeset`, 'OutputFcn'.
%   - Shampine, L. F., & Reichelt, M. W. (1997). The MATLAB ODE Suite.
%     SIAM Journal on Scientific Computing, 18(1), 1-22.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol_out, t_out, stats] = integrate_matlab_ode(rhs_func, tspan, w0, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Configuration and Solver Selection                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('--- Using MATLAB ODE Solver Wrapper ---\n');

    % --- Select Solver ---
    % Determine which built-in MATLAB ODE solver to use based on configuration.
    % Provide a sensible default if the user hasn't specified one.
    if isfield(cfg, 'time') && isfield(cfg.time, 'matlab_solver') && ~isempty(cfg.time.matlab_solver)
        solver_name = cfg.time.matlab_solver;
    else
        solver_name = 'ode113'; % Default: variable-order Adams-Bashforth-Moulton solver, good for moderate accuracy.
        warning('integrate_matlab_ode:DefaultSolver', 'cfg.time.matlab_solver not specified, defaulting to ''%s''.', solver_name);
    end
    % Convert the solver name string into a function handle for calling.
    try
        solver_handle = str2func(solver_name);
    catch ME
        error('integrate_matlab_ode:InvalidSolver', 'Could not find MATLAB ODE solver: %s. Check name and availability.', solver_name);
    end
    fprintf('  Solver: %s\n', solver_name);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup ODE Solver Options (odeset)                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The `odeset` function creates/modifies an options structure that controls
    % various aspects of the ODE solver's behavior.

    % --- Base Options ---
    % Allow the user to provide a base set of options (e.g., for events, mass matrix).
    base_options = odeset(); % Start with MATLAB's default set.
    if isfield(cfg, 'time') && isfield(cfg.time, 'ode_options') && ~isempty(cfg.time.ode_options)
        if isstruct(cfg.time.ode_options)
            base_options = cfg.time.ode_options;
            fprintf('  Using base ODE options provided in cfg.time.ode_options.\n');
        else
            warning('integrate_matlab_ode:InvalidOptions', 'cfg.time.ode_options is not a valid structure. Using default odeset.');
        end
    end

    % --- Tolerances ---
    % Set absolute and relative tolerances for the solver's internal error control.
    % Lower tolerances generally lead to higher accuracy but more computational effort.
    abs_tol = 1e-4; % Default absolute tolerance
    if isfield(cfg, 'time') && isfield(cfg.time, 'AbsTol') && isnumeric(cfg.time.AbsTol)
        abs_tol = cfg.time.AbsTol;
    end
    rel_tol = 1e-4; % Default relative tolerance
    if isfield(cfg, 'time') && isfield(cfg.time, 'RelTol') && isnumeric(cfg.time.RelTol)
        rel_tol = cfg.time.RelTol;
    end
    fprintf('  Tolerances: AbsTol=%.2e, RelTol=%.2e\n', abs_tol, rel_tol);
    % Update the options structure with the specified/default tolerances.
    options = odeset(base_options, 'AbsTol', abs_tol, 'RelTol', rel_tol);

    % --- Progress Bar via OutputFcn ---
    % The 'OutputFcn' option allows calling a function after each successful internal step.
    % We use this to display a text progress bar.
    show_progress = true; % Default to showing the progress bar
    if isfield(cfg, 'time') && isfield(cfg.time, 'show_progress_bar') && islogical(cfg.time.show_progress_bar)
        show_progress = cfg.time.show_progress_bar;
    end

    if show_progress
        % Need to handle potentially pre-existing OutputFcn(s) defined in base_options.
        existing_outputfcn = odeset(options).OutputFcn; % Get current OutputFcn(s)
        if isempty(existing_outputfcn)
            % No existing OutputFcn, just add ours.
            options = odeset(options, 'OutputFcn', @utils.odetpbar);
        elseif isa(existing_outputfcn, 'function_handle')
            % Existing OutputFcn is a single handle, create a cell array including both.
            options = odeset(options, 'OutputFcn', {existing_outputfcn, @utils.odetpbar});
        elseif iscell(existing_outputfcn)
            % Existing OutputFcn is already a cell array, append ours.
            options = odeset(options, 'OutputFcn', [existing_outputfcn(:)', {@utils.odetpbar}]); % Ensure existing is row cell & append
        else
             warning('integrate_matlab_ode:InvalidOutputFcn', 'Existing OutputFcn has an unsupported type. Progress bar might not be added correctly.');
        end
        fprintf('  Progress Bar: Enabled (@utils.odetpbar)\n');
    else
        fprintf('  Progress Bar: Disabled.\n');
    end

    % --- RHS Function Handle ---
    % The rhs_func passed to this wrapper *already* incorporates the 'cfg'
    % structure via an anonymous function created in core.solver. Thus, it
    % correctly matches the f(t, w) signature required by MATLAB ODE solvers.
    ode_rhs = rhs_func;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Call the Selected MATLAB ODE Solver                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('  Starting integration...\n');
    tic; % Start timer

    % Ensure tspan is a column vector if needed by some solvers, although typically
    % row or column works. Row is requested by documentation.
    tspan_vec = tspan(:)';
    % Ensure w0 is a column vector as required by MATLAB ODE solvers.
    w0_col = w0(:);

    % Call the selected solver (e.g., ode113, ode45).
    % The solver integrates from tspan(1) to tspan(end) using adaptive steps,
    % but only returns the solution values at the specific times listed in tspan_vec.
    [T_sol, Y_sol] = solver_handle(ode_rhs, tspan_vec, w0_col, options);

    integration_time = toc; % Stop timer
    fprintf('  Integration finished (%.2f seconds).\n', integration_time);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Format Output to Match Custom Integrator Style             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MATLAB ODE solvers return:
    %   T_sol: Column vector of time points (matching input tspan, potentially rearranged).
    %   Y_sol: Solution matrix where each ROW Y_sol(i,:) corresponds to time T_sol(i).
    % We need:
    %   t_out: Row vector of output times matching the *original input* tspan order.
    %   sol_out: Solution matrix where row sol_out(i,:) corresponds to t_out(i).

    % Ensure the output times match the requested tspan exactly, handling potential
    % minor floating point differences or reordering by the solver.
    if length(T_sol) == length(tspan_vec) && max(abs(sort(T_sol(:)) - sort(tspan_vec(:)))) < 1e-9 * (tspan_vec(end)-tspan_vec(1))
        % If times match reasonably well, use the requested tspan for t_out
        t_out = tspan_vec;
        sol_out = Y_sol;
        % Reorder Y_sol if T_sol wasn't monotonic (rare but possible)
        [~, sort_idx_T] = sort(T_sol(:));
        [~, orig_idx_tspan] = sort(tspan_vec(:));
        if ~isequal(sort_idx_T, orig_idx_tspan)
             [~, remap_idx] = ismember(tspan_vec, T_sol);
             sol_out = Y_sol(remap_idx,:);
        end
    else
        warning('integrate_matlab_ode:TimeMismatch', 'Output times from ODE solver do not precisely match input tspan. Returning solver output times.');
        t_out = T_sol(:)';    % Return solver times as row vector
        sol_out = Y_sol;      % Return corresponding solutions
    end

    % --- Statistics --- (Limited information available for MATLAB solvers here)
    % stats.nsteps: We approximate this as the number of output intervals.
    % This is NOT the number of internal steps taken by the adaptive solver.
    stats.nsteps = length(t_out) - 1;
    % stats.nfevals: Number of RHS evaluations. This is available in the 'stats'
    % output if the solver is called as `sol = odeXX(...)`, but not directly
    % from the `[T,Y] = odeXX(...)` syntax without a custom OutputFcn.
    stats.nfevals = NaN;

    fprintf('--- MATLAB ODE Solver Wrapper Finished ---\n');

end % Function end