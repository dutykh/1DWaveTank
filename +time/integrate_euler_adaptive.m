%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +time/integrate_euler_adaptive.m
%
% Purpose:
%   Solves a system of ODEs dw/dt = rhs_func(t, w, cfg) using the explicit
%   Forward Euler method with adaptive time stepping. The time step dt is
%   determined dynamically at each step using a CFL condition (see
%   core.utils.calculate_dt_cfl). Solution is stored at user-specified output
%   times, and the step size is adjusted to hit these times exactly.
%
% Syntax:
%   [sol_out, t_out, stats] = integrate_euler_adaptive(rhs_func, t_span, w0, cfg)
%
% Inputs:
%   rhs_func - [function handle] RHS of the ODE system. Signature:
%                f = rhs_func(t, w, cfg)
%   t_span   - [1 x 2, double] Start and end times [t0, tf].
%   w0       - [vector, double] Initial state vector at t0.
%   cfg      - [struct] Configuration structure. Must contain:
%                cfg.time.dt_plot: Output interval [s]
%                cfg.time.cfl: CFL number for adaptive time stepping
%                cfg.mesh: Mesh configuration needed by calculate_dt_cfl
%                cfg.phys: Physical parameters needed by calculate_dt_cfl
%
% Outputs:
%   sol_out  - [M x length(w0), double] Solution matrix at output times.
%   t_out    - [M x 1, double] Column vector of output times.
%   stats    - [struct] Statistics:
%                stats.nsteps:   Total number of internal time steps taken
%                stats.nfevals:  Total number of RHS evaluations
%
% Dependencies:
%   - core.utils.calculate_dt_cfl.m (for adaptive time step)
%
% References:
%   - LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems.
%     Cambridge University Press. (Chapter 6)
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol_out, t_out, stats] = integrate_euler_adaptive(rhs_func, t_span, w0, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Preparation and Output Allocation                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(cfg.time, 'dt_plot') || isempty(cfg.time.dt_plot) || cfg.time.dt_plot <= 0
        warning('cfg.time.dt_plot not set or invalid, using default of 0.1s');
        cfg.time.dt_plot = 0.1;
    end

    t0 = t_span(1);
    tf = t_span(end);
    dt_plot = cfg.time.dt_plot; % Desired interval between output points

    % Estimate the maximum number of output steps needed
    max_output_steps = ceil((tf - t0) / dt_plot) + 2; % +1 for t0, +1 for potential tf rounding
    t_out = zeros(max_output_steps, 1);   % Preallocate time output vector
    num_vars = length(w0);
    w_out = zeros(max_output_steps, num_vars); % Preallocate solution output matrix

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Store Initial Condition                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = t0;
    w = w0(:)'; % Ensure w is a row vector for consistency
    output_count = 1;
    t_out(output_count) = t;
    w_out(output_count, :) = w;

    % Determine the next target time for output
    t_plot_next = t0 + dt_plot;
    final_time_target = tf;

    step = 0;
    k = 0; % Step counter
    dt_history = zeros(1, 10000); % Preallocate dt history
    TOL = 1e-9 * max(dt_plot, tf-t0); % Tolerance for floating point comparisons relative to timescale
    fprintf('Starting adaptive Euler integration from t=%.3f to t=%.3f, plotting every %.3f s\n', t0, tf, dt_plot);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Progress Reporting Setup                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_reports = 10; % Default number of reports
    if isfield(cfg, 'time') && isfield(cfg.time, 'num_progress_reports') && cfg.time.num_progress_reports > 0
        num_reports = cfg.time.num_progress_reports;
    end
    report_interval = (tf - t0) / num_reports; % Report progress roughly num_reports times
    next_report_time = t0 + report_interval;
    last_report_perc = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main Time Stepping Loop                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while t < final_time_target - TOL % Loop until very close to the final target time
        step = step + 1;

        % --- Calculate adaptive timestep based on current state w ---
        w = w(:); % Ensure column vector for CFL calculation
        dt_adaptive = core.utils.calculate_dt_cfl(w, cfg);

        % --- Determine timestep 'dt' for the current step ---
        dt_to_next_output = min(t_plot_next - t, final_time_target - t);
        dt_to_next_output = max(dt_to_next_output, 0); % Ensure non-negative
        hit_output_time = false;

        % Decide the step size
        if dt_adaptive >= dt_to_next_output - TOL
            dt = dt_to_next_output;
            hit_output_time = true;
        else
            dt = dt_adaptive;
            if abs((t + dt) - t_plot_next) < TOL || abs((t + dt) - final_time_target) < TOL
                hit_output_time = true;
            end
        end

        % --- Compute Forward Euler Step ---
        f_n = rhs_func(t, w, cfg);
        w_new = w + dt * f_n; % Note: f_n must be a column vector

        t_new = t + dt;
        k = k + 1;
        dt_history(k) = dt;

        % --- Progress Reporting ---
        if report_interval > 0 && t >= next_report_time
            perc_done = round(100 * (next_report_time - t0) / (tf - t0));
            % Avoid printing the same percentage twice if steps are small
            if perc_done > last_report_perc || perc_done == 0
               fprintf('\n  t = %.3f s (%.1f%%), dt = %.3e s', next_report_time, perc_done, dt); % Match RK4/BS format
               next_report_time = next_report_time + report_interval;
               % Ensure next report time doesn't slightly exceed t_end due to float arithmetic
               if next_report_time > tf
                   next_report_time = tf + 1; % Effectively disable further reports
               end
               last_report_perc = perc_done;
            end
         end
 
         % --- Output Handling: Store Solution at Requested Times ---
        if hit_output_time
            output_count = output_count + 1;
            t_out(output_count) = t_plot_next;
            w_out(output_count, :) = w_new';
            t_plot_next = t_plot_next + dt_plot;
        end

        % --- Prepare for Next Step ---
        w = w_new;
        t = t_new;
    end % End while loop

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Final State Storage Check                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ensure the final state at tf is stored if the loop terminated slightly before
    % or exactly at tf, but it wasn't recorded in the last iteration.
    if abs(t - tf) < TOL && abs(t_out(output_count) - tf) > TOL
        output_count = output_count + 1;
        t_out(output_count) = tf;
        w_out(output_count, :) = w';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Output Formatting and Statistics                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_out = t_out(1:output_count);
    sol_out = w_out(1:output_count, :);
    stats.nsteps = k;
    stats.nfevals = k; % Each step has one RHS evaluation

end