%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +time/integrate_rk4_adaptive.m
%
% Purpose:
%   Solves a system of ODEs dw/dt = rhs_func(t, w, cfg) using the classic
%   explicit 4th-order Runge-Kutta (RK4) method with adaptive time stepping.
%   While RK4 itself is fixed-order, the time step `dt` is adapted at each step
%   based on the CFL condition to ensure numerical stability for hyperbolic problems.
%   Solution is stored at user-specified output times, and the step size is
%   adjusted to hit these times exactly.
%
% Syntax:
%   [sol_out, t_out, stats] = integrate_rk4_adaptive(rhs_func, tspan, w0, cfg)
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
%                stats.nsteps:   Total number of internal RK4 time steps taken.
%                stats.nfevals:  Total number of RHS evaluations (4 * nsteps for RK4).
%
% Dependencies:
%   - core.utils.calculate_dt_cfl.m (for adaptive time step)
%   - Progress bar utility (optional)
%
% References:
%   - Butcher, J. C. (2008). Numerical Methods for Ordinary Differential Equations (2nd ed.). Wiley.
%   - LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations.
%     SIAM.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol_out, t_out, stats] = integrate_rk4_adaptive(rhs_func, tspan, w0, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Validation and Setup                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 4
        error('integrate_rk4_adaptive:NotEnoughInputs', 'Not enough input arguments.');
    end
    if ~isfield(cfg, 'time') || ~isfield(cfg.time, 'cfl') || isempty(cfg.time.cfl) || cfg.time.cfl <= 0
        error('integrate_rk4_adaptive:MissingCFL', 'CFL number must be specified and positive in cfg.time.cfl');
    end
    if ~isvector(tspan) || ~issorted(tspan) || tspan(1) < 0 || length(tspan) < 2
        error('integrate_rk4_adaptive:InvalidTSPAN', 'TSPAN must be a monotonically increasing vector with at least two elements, starting from t0 >= 0.');
    end

    t0 = tspan(1);         % [s] Initial time
    tf = tspan(end);       % [s] Final time
    t_out_req = tspan(:)'; % Ensure requested output times is a row vector

    w = w0(:); % Ensure w0 is a column vector for internal calculations
    t = t0;    % [s] Current simulation time
    k = 0;     % Step counter
    nfevals = 0; % RHS evaluation counter

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preallocate Output Arrays                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_outputs = length(t_out_req);
    % Store solution as columns initially for easier concatenation if resize needed
    sol_out_internal = zeros(length(w0), num_outputs);
    t_out = zeros(1, num_outputs);
    % Estimate max steps for dt_history (can be resized if needed)
    max_diff_val = max(max(diff(t_out_req), 1e-6)); % Get the single maximum value
    estimated_steps = ceil(10 * (tf - t0) / max_diff_val) + 100; 
    dt_history = zeros(1, estimated_steps);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Store Initial Condition                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    output_idx = 1;
    if abs(t - t_out_req(output_idx)) < 1e-12 % Check if t0 is the first output time
        sol_out_internal(:, output_idx) = w;
        t_out(output_idx) = t;
        output_idx = output_idx + 1;
    end
    if num_outputs >= output_idx
         t_next_plot = t_out_req(output_idx);
    else
         t_next_plot = tf + 1; % No more plotting needed
    end

    fprintf('Starting adaptive RK4 integration from t=%.3f to t=%.3f\n', t0, tf);
    fprintf('Output requested at %d time points.\n', num_outputs);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Progress Reporting Setup                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    last_report_time = t0;
    num_reports = 10; % Default number of reports
    if isfield(cfg, 'time') && isfield(cfg.time, 'num_progress_reports') && cfg.time.num_progress_reports > 0
        num_reports = cfg.time.num_progress_reports;
    end
    report_interval = (tf - t0) / num_reports; % Report progress roughly num_reports times

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main Time Stepping Loop                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    max_internal_steps = 1e7; % Safety break
    while t < tf
        if k >= max_internal_steps
             warning('integrate_rk4_adaptive:MaxStepsExceeded', 'Maximum internal steps (%d) exceeded. Aborting.', max_internal_steps);
             break;
        end

        % --- Calculate Adaptive Time Step Based on CFL ---
        dt = core.utils.calculate_dt_cfl(w, cfg);

        % --- Adjust dt to Hit Output Times Exactly ---
        dt_to_tf = tf - t;
        dt_to_plot = t_next_plot - t;
        % Choose smallest of CFL dt, time to tf, time to next plot
        dt = min([dt, dt_to_tf, dt_to_plot]);

        % --- Safety Checks for dt ---
        if dt <= 1e-12 % Prevent excessively small steps
            if abs(t-tf) < 1e-9
                 fprintf('Reached final time tf=%.4f\n', tf);
                 break; % Exit loop if effectively at the end time
            else
                warning('integrate_rk4_adaptive:SmallDt', 'Time step dt=%.3e is too small at t=%.3f. Aborting integration.', dt, t);
                break;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Classic 4th-Order Runge-Kutta (RK4) Stages              %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Formula: w_{n+1} = w_n + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
        % where k_i are estimates of the slope at different points in the interval.
        k1 = rhs_func(t,             w,            cfg); % Slope at the beginning
        k2 = rhs_func(t + 0.5 * dt, w + 0.5 * dt*k1, cfg); % Slope at midpoint using k1
        k3 = rhs_func(t + 0.5 * dt, w + 0.5 * dt*k2, cfg); % Slope at midpoint using k2
        k4 = rhs_func(t + dt,       w + dt*k3,       cfg); % Slope at the end using k3
        nfevals = nfevals + 4; % Increment RHS evaluation count

        % --- Update Solution and Time ---
        w_new = w + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        t_new = t + dt;
        k = k + 1;

        % --- Store dt History (Resize if needed) ---
        if k > length(dt_history)
             dt_chunk_size = length(dt_history); % Double the current size
             warning('Time:Integrate:GrowStats', 'Growing dt_history size at step %d (t=%.3f)', k, t);
             dt_history(end+1 : end+dt_chunk_size) = 0; % Grow using direct indexing
        end
        dt_history(k) = dt;

        % --- Progress Reporting ---
        if report_interval > 0 && t_new >= last_report_time + report_interval
            fprintf('  t = %.3f s (%.1f%%), dt = %.3e s\n', t_new, 100*(t_new-t0)/(tf-t0), dt);
            last_report_time = t_new;
        end

        % --- Output Handling: Store Solution at Requested Times ---
        % Check if the new time step landed exactly on a requested output time.
        if abs(t_new - t_next_plot) < 1e-12
            sol_out_internal(:, output_idx) = w_new;
            t_out(output_idx) = t_new;
            output_idx = output_idx + 1;
            if output_idx <= num_outputs
                 t_next_plot = t_out_req(output_idx);
            else
                 t_next_plot = tf + 1; % No more outputs needed
            end
        end

        % --- Prepare for Next Step ---
        w = w_new;
        t = t_new;

    end % End while loop

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Final Output Formatting and Statistics                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Trim unused preallocated space
    sol_out_internal = sol_out_internal(:, 1:output_idx-1);
    t_out = t_out(1:output_idx-1);
    dt_history = dt_history(1:k);

    % Transpose solution to match expected output format [M x length(w0)]
    sol_out = sol_out_internal';

    % --- Statistics ---
    stats.nsteps = k;
    stats.nfevals = nfevals;
    % stats.dt_history = dt_history; % Optionally return dt history

    fprintf('Integration finished at t = %.3f s after %d steps.\n', t_out(end), k);

end % Function end