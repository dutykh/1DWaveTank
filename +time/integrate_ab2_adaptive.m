%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +time/integrate_ab2_adaptive.m
%
% Purpose:
%   Solves a system of ODEs dw/dt = rhs_func(t, w, cfg) using the explicit
%   2nd-order Adams-Bashforth (AB2) method with adaptive time stepping.
%   The time step dt is determined dynamically at each step using a CFL
%   condition (see core.utils.calculate_dt_cfl). The first step is performed
%   using Forward Euler to initialize the AB2 method.
%
% Syntax:
%   [sol_out, t_out, stats] = integrate_ab2_adaptive(rhs_func, t_span, w0, cfg)
%
% Inputs:
%   rhs_func - [function handle] RHS of the ODE system. Signature:
%                f = rhs_func(t, w, cfg)
%   t_span   - [1 x M, double] Output times [t0, t1, ..., tf].
%   w0       - [vector, double] Initial state vector at t0.
%   cfg      - [struct] Configuration structure. Must contain:
%                cfg.phys.g, cfg.time.cfl, cfg.mesh.N, cfg.mesh.dx,
%                cfg.time.num_progress_reports (for progress bar).
%
% Outputs:
%   sol_out  - [M x length(w0), double] Solution matrix at output times.
%   t_out    - [1 x M, double] Row vector of output times.
%   stats    - [struct] Statistics:
%                stats.nsteps:   Total number of internal time steps taken
%                stats.nfevals:  Total number of RHS evaluations
%
% Dependencies:
%   - core.utils.calculate_dt_cfl.m (for adaptive time step)
%   - Progress bar utility (optional)
%
% References:
%   - LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems.
%     Cambridge University Press. (Chapter 6)
%   - Standard ODE texts for Adams-Bashforth methods.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol_out, t_out, stats] = integrate_ab2_adaptive(rhs_func, t_span, w0, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Preparation and Output Allocation                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t0 = t_span(1);
    tf = t_span(end);
    t_out_req = t_span; % Requested output times
    num_out_points = length(t_out_req);

    % Preallocate output arrays conservatively (estimate steps)
    max_steps = 1e7; % Set a large number for safety break
    sol_out = zeros(num_out_points, length(w0));
    t_out = zeros(1, num_out_points);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialization                                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = t0;
    w = w0;
    k = 0; % Step counter
    output_count = 1; % Index for storing output
    nfevals = 0;
    f_prev = []; % Store previous RHS evaluation

    % Store initial condition
    sol_out(output_count,:) = w';
    t_out(output_count) = t;
    output_count = output_count + 1;
    if num_out_points > 1
        t_next_plot = t_out_req(output_count);
    else
        t_next_plot = tf + 1; % No intermediate plotting if only t0, tf requested
    end

    % Progress reporting setup
    report_interval = (tf - t0) / cfg.time.num_progress_reports;
    next_report_time = t0 + report_interval;
    fprintf('Starting adaptive AB2 integration from t=%.3f to t=%.3f\n', t0, tf);
    fprintf('Output requested at %d time points.\n', num_out_points);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main Time Stepping Loop                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while t < tf
        % --- Safety Break: Prevent Infinite Loops ---
        if k >= max_steps
            warning('Maximum number of steps (%d) exceeded. Aborting integration at t=%.3f s.', max_steps, t);
            break;
        end

        % --- Calculate Adaptive Time Step Based on CFL ---
        dt = core.utils.calculate_dt_cfl(w, cfg);

        % --- Prevent Overshooting Final Time ---
        if t + dt > tf
            dt = tf - t;
        end

        % --- Prevent Overshooting Next Output Time ---
        if t + dt > t_next_plot && t < t_next_plot
            dt = t_next_plot - t; % Step exactly to the plot time
        end

        % --- Check for Excessively Small Time Step ---
        if dt < 1e-12
            warning('Time step dt=%.3e is too small at t=%.3f. Aborting integration.', dt, t);
            break;
        end

        % --- Compute RHS (f_n) ---
        f_n = rhs_func(t, w, cfg); nfevals = nfevals + 1;

        % --- First Step: Use Forward Euler (AB2 needs two points) ---
        if k == 0
            w_new = w + dt * f_n;
        else
            % AB2 Formula: w_{n+1} = w_n + dt/2 * (3*f_n - f_{n-1})
            w_new = w + dt/2 * (3*f_n - f_prev);
        end

        % --- Advance Time ---
        t_new = t + dt;
        k = k + 1;

        % --- Progress Reporting ---
        if t_new >= next_report_time
            fprintf('  t = %.3f s (%.1f%%), dt = %.3e s\n', t_new, 100*(t_new-t0)/(tf-t0), dt);
            next_report_time = next_report_time + report_interval;
        end

        % --- Output Handling: Store Solution at Requested Times ---
        if t_new >= t_next_plot - 1e-12
            sol_out(output_count,:) = w_new';
            t_out(output_count) = t_next_plot;
            output_count = output_count + 1;
            if output_count <= num_out_points
                t_next_plot = t_out_req(output_count);
            else
                t_next_plot = tf + 1; % No more outputs required
            end
        end

        % --- Prepare for Next Step ---
        w = w_new;
        t = t_new;
        f_prev = f_n;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Output Statistics                                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stats.nsteps = k;
    stats.nfevals = nfevals;

end
