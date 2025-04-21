%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +time/integrate_ssp2_adaptive.m
%
% Purpose:
%   Solves a system of ODEs dw/dt = rhs_func(t, w, cfg) using the explicit
%   2nd-order Strong Stability Preserving Runge-Kutta method (SSP2, also
%   known as Heun's method or improved Euler) with adaptive time stepping.
%   SSP methods are designed to preserve the stability properties (e.g., TVD)
%   of the spatial discretization when coupled with Forward Euler time stepping.
%   The time step dt is determined dynamically at each step using a CFL
%   condition (see core.utils.calculate_dt_cfl). Solution is stored at
%   user-specified output times, adjusting step size to hit these times exactly.
%
% Syntax:
%   [sol_out, t_out, stats] = integrate_ssp2_adaptive(rhs_func, tspan, w0, cfg)
%
% Inputs:
%   rhs_func - [function handle] RHS of the ODE system. Signature:
%                f = rhs_func(t, w, cfg)
%   tspan    - [vector, double] Time points [t0, t1, ..., tf] at which the
%                solution output is requested. Must be monotonically increasing.
%   w0       - [vector, double] Initial state vector (column vector) at time t0.
%   cfg      - [struct] Configuration structure. Must contain:
%                cfg.phys.g, cfg.time.cfl, cfg.mesh.N, cfg.mesh.dx,
%                cfg.time.num_progress_reports (for progress bar).
%
% Outputs:
%   sol_out  - [M x length(w0), double] Solution matrix. Each row `sol_out(i,:)`
%                is the state vector corresponding to the time point `t_out(i)`.
%   t_out    - [1 x M, double] Row vector of output times.
%   stats    - [struct] Statistics:
%                stats.nsteps:   Total number of internal SSP2 time steps taken.
%                stats.nfevals:  Total number of RHS evaluations (2 * nsteps for SSP2).
%
% Dependencies:
%   - core.utils.calculate_dt_cfl.m (for adaptive time step)
%   - Progress bar utility (optional)
%
% References:
%   - Gottlieb, S., Shu, C.-W., & Tadmor, E. (2001). Strong Stability-Preserving
%     High-Order Time Discretization Methods. SIAM Review, 43(1), 89-112.
%   - LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems.
%     Cambridge University Press.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol_out, t_out, stats] = integrate_ssp2_adaptive(rhs_func, tspan, w0, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Validation and Setup                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Validate input arguments
    if nargin < 4
        error('integrate_ssp2_adaptive:NotEnoughInputs', 'Not enough input arguments.');
    end
    if ~isfield(cfg, 'time') || ~isfield(cfg.time, 'cfl') || isempty(cfg.time.cfl) || cfg.time.cfl <= 0
        error('integrate_ssp2_adaptive:MissingCFL', 'CFL number must be specified and positive in cfg.time.cfl');
    end
    if ~isvector(tspan) || ~issorted(tspan) || tspan(1) < 0 || length(tspan) < 2
        error('integrate_ssp2_adaptive:InvalidTSPAN', 'TSPAN must be a monotonically increasing vector with at least two elements, starting from t0 >= 0.');
    end

    %% Extract initial and final times
    t0 = tspan(1);         % [s] Initial time
    tf = tspan(end);       % [s] Final time
    t_out_req = tspan(:)'; % Ensure requested output times is a row vector

    %% Initialize state vector and time
    w = w0(:); % Ensure w0 is a column vector for internal calculations
    t = t0;    % [s] Current simulation time
    k = 0;     % Step counter
    nfevals = 0; % RHS evaluation counter

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preallocate Output Arrays                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Preallocate output arrays
    num_outputs = length(t_out_req);
    % Store solution as rows for direct compatibility with solver output format
    sol_out = zeros(num_outputs, length(w));
    t_out = zeros(1, num_outputs);
    % Estimate max steps for dt_history (can be resized if needed)
    max_diff_val = max(max(diff(t_out_req), 1e-6)); % Get the single maximum value
    estimated_steps = ceil(10 * (tf - t0) / max_diff_val) + 100;
    dt_history = zeros(1, estimated_steps);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Store Initial Condition                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Store initial condition
    output_idx = 1;
    if abs(t - t_out_req(output_idx)) < 1e-12 % Check if t0 is the first output time
        sol_out(output_idx,:) = w'; % Store initial state (as row)
        t_out(output_idx) = t;
        output_idx = output_idx + 1;
    end
    if num_outputs >= output_idx
         t_next_plot = t_out_req(output_idx);
    else
         t_next_plot = tf + 1; % No more plotting needed
    end

    %% Print simulation details
    fprintf('Starting adaptive SSP(2,2) integration from t=%.3f to t=%.3f\n', t0, tf);
    fprintf('Output requested at %d time points.\n', num_outputs);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Progress Reporting Setup                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set up progress reporting
    last_report_time = t0;
    num_reports = 10; % Default number of reports
    if isfield(cfg, 'time') && isfield(cfg.time, 'num_progress_reports') && cfg.time.num_progress_reports > 0
        num_reports = cfg.time.num_progress_reports;
    end
    report_interval = (tf - t0) / num_reports; % Report progress roughly num_reports times

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main Time Stepping Loop                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Main time stepping loop
    max_internal_steps = 1e7; % Safety break
    while t < tf
        if k >= max_internal_steps
             warning('integrate_ssp2_adaptive:MaxStepsExceeded', 'Maximum internal steps (%d) exceeded. Aborting.', max_internal_steps);
             break;
        end

        % --- Calculate Adaptive Time Step Based on CFL ---
        %% Calculate adaptive time step based on CFL condition
        % The CFL condition is used to determine the maximum time step size
        % that maintains stability. This is done by calculating the CFL number
        % based on the current state and configuration.
        dt = core.utils.calculate_dt_cfl(w, cfg);

        % --- Adjust dt to Hit Output Times Exactly ---
        %% Adjust dt to hit output times exactly
        % To ensure that the solution is stored at the requested output times,
        % the time step size is adjusted to hit these times exactly.
        dt_to_tf = tf - t;
        dt_to_plot = t_next_plot - t;
        % Choose smallest of CFL dt, time to tf, time to next plot
        dt = min([dt, dt_to_tf, dt_to_plot]);

        % --- Safety Checks for dt ---
        %% Safety checks for dt
        % To prevent excessively small time steps, a minimum time step size
        % is enforced. If the time step size is smaller than this minimum,
        % the simulation is aborted.
        if dt <= 1e-12 % Prevent excessively small steps
            if abs(t-tf) < 1e-9
                 fprintf('Reached final time tf=%.4f\n', tf);
                 break; % Exit loop if effectively at the end time
            else
                warning('integrate_ssp2_adaptive:SmallDt', 'Time step dt=%.3e is too small at t=%.3f. Aborting integration.', dt, t);
                break;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Strong Stability Preserving RK(2,2) Stages              %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Strong stability preserving RK(2,2) stages
        % This method combines Forward Euler steps to achieve 2nd order
        % accuracy while maintaining stability properties under a CFL limit.
        % Formula:
        %   w^(1) = w^n + dt * F(t^n, w^n)
        %   w^(n+1) = 0.5*w^n + 0.5*( w^(1) + dt*F(t^n+dt, w^(1)) )

        % --- Stage 1: Forward Euler Step ---
        %% Stage 1: Forward Euler step
        % The first stage of the RK(2,2) method is a Forward Euler step.
        F_n = rhs_func(t, w, cfg);       % F(t^n, w^n)
        w1 = w + dt * F_n;               % w^(1)
        nfevals = nfevals + 1;

        % --- Stage 2: Combine Stages ---
        %% Stage 2: Combine stages
        % The second stage of the RK(2,2) method combines the results of the
        % first stage with the current state to produce the final state.
        F_1 = rhs_func(t + dt, w1, cfg); % F(t^n+dt, w^(1))
        w_new = 0.5*w + 0.5*(w1 + dt*F_1);
        nfevals = nfevals + 1;

        % --- Advance Time ---
        %% Advance time
        t_new = t + dt;
        k = k + 1;

        % --- Store dt History (Resize if needed) ---
        %% Store dt history (resize if needed)
        if k > length(dt_history)
             dt_history = [dt_history, zeros(1, estimated_steps)];
        end
        dt_history(k) = dt;

        % --- Progress Reporting ---
        %% Progress reporting
        if report_interval > 0 && t_new >= last_report_time + report_interval
            fprintf('  t = %.3f s (%.1f%%), dt = %.3e s\n', t_new, 100*(t_new-t0)/(tf-t0), dt);
            last_report_time = t_new;
        end

        % --- Output Handling: Store Solution at Requested Times ---
        %% Output handling: Store solution at requested times
        % Check if the new time step landed exactly on a requested output time.
        if abs(t_new - t_next_plot) < 1e-12
            sol_out(output_idx,:) = w_new'; % Store as row vector
            t_out(output_idx) = t_new;
            output_idx = output_idx + 1;
            if output_idx <= num_outputs
                 t_next_plot = t_out_req(output_idx);
            else
                 t_next_plot = tf + 1; % No more outputs needed
            end
        end

        % --- Prepare for Next Step ---
        %% Prepare for next step
        w = w_new;
        t = t_new;

    end % End while loop

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Final Output Formatting and Statistics                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Final output formatting and statistics
    % Trim unused preallocated space
    sol_out = sol_out(1:output_idx-1, :);
    t_out = t_out(1:output_idx-1);
    dt_history = dt_history(1:k);

    % --- Statistics ---
    %% Statistics
    stats.nsteps = k;
    stats.nfevals = nfevals;
    % stats.dt_history = dt_history; % Optionally return dt history

    fprintf('Integration finished at t = %.3f s after %d steps.\n', t_out(end), k);

end % Function end