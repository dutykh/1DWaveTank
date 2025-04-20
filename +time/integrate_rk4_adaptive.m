function [sol_out, t_out, stats] = integrate_rk4_adaptive(rhs_func, tspan, w0, cfg)

    % INTEGRATE_RK4_ADAPTIVE Solves ODEs using the classic 4th-order Runge-Kutta method.
    %   Integrates the system using a fixed-order RK4 method, but adapts the
    %   time step based on the CFL condition for stability.
    %
    %   [T_OUT, SOL_OUT, K, DT_HISTORY] = INTEGRATE_RK4_ADAPTIVE(RHS_FUNC, TSPAN, W0, CFG)
    %   integrates the system defined by RHS_FUNC over TSPAN with initial
    %   condition W0.
    %
    %   Inputs:
    %   RHS_FUNC   - Function handle for the right-hand side, expected to have the
    %                signature rhs_func(t, w, cfg).
    %   TSPAN      - Vector specifying the time points for output [t0, t1, ..., tf].
    %                Must be monotonically increasing.
    %   W0         - Initial condition vector (column vector).
    %   CFG        - Configuration structure, must contain CFG.TIME.CFL.
    %
    %   Outputs:
    %   T_OUT      - Row vector of time points corresponding to the solution points.
    %   SOL_OUT    - Matrix of solution vectors at times in T_OUT. Each column
    %                corresponds to a time point.
    %   K          - Total number of time steps taken.
    %   DT_HISTORY - Vector containing the time step size used at each step.

    % --- Input Validation and Setup ---
    if nargin < 4
        error('Not enough input arguments.');
    end
    if ~isfield(cfg, 'time') || ~isfield(cfg.time, 'CFL')
        error('CFL number must be specified in cfg.time.CFL');
    end
    if ~isvector(tspan) || ~issorted(tspan) || tspan(1) < 0
        error('TSPAN must be a monotonically increasing vector starting from t0 >= 0.');
    end

    t0 = tspan(1);
    tf = tspan(end);
    t_out_req = tspan; % Requested output times

    w = w0(:); % Ensure w0 is a column vector
    t = t0;
    k = 0;     % Step counter

    % --- Preallocate Output Arrays --- 
    num_outputs = length(t_out_req);
    sol_out = zeros(length(w), num_outputs);
    t_out = zeros(1, num_outputs);
    dt_history = zeros(1, 10000); % Preallocate reasonable size, will grow if needed

    % --- Store Initial Condition --- 
    output_idx = 1;
    if t >= t_out_req(output_idx)
        sol_out(:, output_idx) = w;
        t_out(output_idx) = t;
        output_idx = output_idx + 1;
    end

    fprintf('Starting RK4 integration from t=%.3f to t=%.3f\n', t0, tf);
    fprintf('Output requested at %d time points.\n', num_outputs);

    last_report_time = t0;
    num_reports = 10; % Default number of reports
    if isfield(cfg, 'time') && isfield(cfg.time, 'num_progress_reports') && cfg.time.num_progress_reports > 0
        num_reports = cfg.time.num_progress_reports;
    end
    report_interval = (tf - t0) / num_reports; % Report progress roughly num_reports times

    % --- Time Stepping Loop --- 
    while t < tf
        % Calculate adaptive time step based on CFL
        dt = core.utils.calculate_dt_cfl(w, cfg);

        % Prevent dt from overshooting tf
        if t + dt > tf
            dt = tf - t; 
        end
        
        % Prevent dt from overshooting the next output time if very close
        if output_idx <= num_outputs && (t + dt > t_out_req(output_idx))
             dt_to_output = t_out_req(output_idx) - t;
             if dt_to_output < dt && dt_to_output > 1e-12 % Avoid tiny steps due to float precision
                 dt = dt_to_output; 
             end
        end

        % Check for invalid dt
        if dt <= 0 || isnan(dt)
            error('Invalid time step dt = %.4e at t = %.4f. Simulation unstable?', dt, t);
        end

        % --- RK4 Stages --- 
        k1 = dt * rhs_func(t,             w,            cfg);
        k2 = dt * rhs_func(t + 0.5 * dt, w + 0.5 * k1, cfg);
        k3 = dt * rhs_func(t + 0.5 * dt, w + 0.5 * k2, cfg);
        k4 = dt * rhs_func(t + dt,       w + k3,       cfg);

        % --- Update Solution and Time ---
        w = w + (k1 + 2*k2 + 2*k3 + k4) / 6;
        t = t + dt;
        k = k + 1;

        % Store dt history (expand array if needed)
        if k > length(dt_history)
            dt_history = [dt_history, zeros(1, length(dt_history))];
        end
        dt_history(k) = dt;

        % --- Store Output ---
        while output_idx <= num_outputs && t >= t_out_req(output_idx) - 1e-9 % Tolerate small floating point errors
            sol_out(:, output_idx) = w;
            t_out(output_idx) = t; % Store the actual time reached
            output_idx = output_idx + 1;
        end
        
        % --- Progress Reporting ---
        if report_interval > 0 && t - last_report_time >= report_interval
            fprintf('  t = %.3f s (%.1f%%), dt = %.3e s\n', t, (t/tf)*100, dt);
            last_report_time = t;
        end

        % Safety break for excessive steps
        if k > 1e7
            warning('Exceeded maximum number of steps (1e7). Stopping integration.');
            break;
        end
    end

    % Trim unused parts of history arrays
    dt_history = dt_history(1:k);
    % Ensure final point is included if loop terminated exactly at tf
    if output_idx <= num_outputs && abs(t - tf) < 1e-9
         sol_out(:, output_idx) = w;
         t_out(output_idx) = t;
    end
    % Trim unused output slots if loop finished early (e.g., due to instability)
    sol_out = sol_out(:, 1:find(t_out > 0, 1, 'last'));
    t_out = t_out(1:size(sol_out, 2));
    t_out = t_out(:)'; % Ensure row vector

    fprintf('RK4 integration finished at t = %.3f s after %d steps.\n', t, k);
    if abs(t - tf) > 1e-6
        warning('Integration did not reach the final time tf = %.3f. Stopped at t = %.3f.', tf, t);
    end

    % --- Transpose output to match expected format (Time x State) --- 
    sol_out = sol_out';
    stats = struct('nsteps', k, 'nfevals', k);

end % Function end