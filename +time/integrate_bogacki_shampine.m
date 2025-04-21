function [sol_out, t_out, stats] = integrate_bogacki_shampine(rhs_handle, tspan, w0, cfg)
%INTEGRATE_BOGACKI_SHAMPINE Integrates an ODE system using Bogacki-Shampine.
%   Integrates the system defined by rhs_handle over tspan = [t_start, t_end]
%   starting from the initial condition w0, using the adaptive BS3(2) method.
%
%   Args:
%       rhs_handle: Function handle for the RHS, expected signature f(t, w, cfg).
%       tspan:      Time interval [t_start, t_end].
%       w0:         Initial state vector at t_start.
%       cfg:        Configuration struct (passed to rhs_handle and for settings).
%
%   Returns:
%       sol_out:    Solution matrix (each row is state at time t_out(i)).
%       t_out:      Column vector of time points.
%       stats:      Struct with integration statistics (e.g., nsteps).

    % --- Configuration --- % 
    t_start = tspan(1);
    t_end   = tspan(end);
    dt_init = cfg.time.dt_init; % Initial time step guess
    tol     = cfg.time.adaptive_tol; % Desired tolerance
    safety  = 0.9;  % Safety factor for step size adjustment
    min_factor = 0.2; % Minimum step size reduction factor
    max_factor = 5.0; % Maximum step size increase factor
    order   = 3; % Order of the method for step size control (p+1 for p=2)
    max_steps = 1e6; % Maximum allowed steps to prevent infinite loops
    dt_min = 1e-12 * (t_end - t_start); % Minimum allowed dt
    dt_plot = cfg.time.dt_plot; % How often to store results

    % --- Bogacki-Shampine Coefficients --- %
    a = [0; 1/2; 3/4; 1];
    b = [ ...
        0,   0,   0,  0; ...
        1/2, 0,   0,  0; ...
        0,   3/4, 0,  0; ...
        2/9, 1/3, 4/9, 0 ...
    ];
    c3 = [2/9, 1/3, 4/9, 0];   % 3rd order coefficients
    c2 = [7/24, 1/4, 1/3, 1/8]; % 2nd order coefficients (embedding)

    % --- Initialization --- %
    t = t_start;
    w = w0;
    dt = dt_init;
    step = 0;
    
    % Preallocate output arrays (estimate size, will grow if needed)
    n_vars = length(w0);
    n_out_est = ceil((t_end - t_start) / dt_plot) + 1;
    t_out = zeros(n_out_est, 1);
    sol_out = zeros(n_out_est, n_vars);
    n_out_points = length(0:dt_plot:t_end) + (mod(t_end, dt_plot)~=0); % Calculate number of output points
    
    % Store initial condition
    out_idx = 1;
    t_out(out_idx) = t;
    sol_out(out_idx, :) = w'; % Store as row vector
    next_plot_time = t + dt_plot;

    % Progress reporting setup
    num_reports = cfg.time.num_progress_reports; % Get number of reports from config
    report_interval = (t_end - t_start) / num_reports;
    next_report_time = t_start + report_interval;
    last_report_perc = 0;

    fprintf('Starting Bogacki-Shampine integration from t=%.3f to t=%.3f (tol=%.1e)', t_start, t_end, tol);
    fprintf('\nOutput requested at %d time points.', n_out_points);
 
    % --- Main Time Loop --- %
    while t < t_end && step < max_steps
        step = step + 1;
        
        % Ensure dt doesn't overshoot t_end
        if t + dt > t_end
            dt = t_end - t;
        end

        % --- Adaptive Step Calculation (Inner Loop) --- %
        accepted = false;
        while ~accepted
            if dt < dt_min
                 warning('Time step (%g) below minimum (%g) at t=%g. Proceeding, but may lose accuracy or fail.', dt, dt_min, t);
                 dt = dt_min; % Prevent infinite loop in extreme cases
                 % Optionally: error('Time step too small...');
            end

            % Calculate stages (k_i)
            k = cell(4, 1);
            try
                k{1} = rhs_handle(t + a(1)*dt, w, cfg); % Call user's RHS function
                k{2} = rhs_handle(t + a(2)*dt, w + dt * (b(2,1)*k{1}), cfg);
                k{3} = rhs_handle(t + a(3)*dt, w + dt * (b(3,1)*k{1} + b(3,2)*k{2}), cfg);

                % Calculate 3rd order solution (higher order)
                w_next_3 = w + dt * (c3(1)*k{1} + c3(2)*k{2} + c3(3)*k{3});

                % Calculate the 4th stage using the 3rd order result
                k{4} = rhs_handle(t + a(4)*dt, w_next_3, cfg); % Use 3rd order result for k4

                % Calculate 2nd order solution (lower order embedding)
                w_next_2 = w + dt * (c2(1)*k{1} + c2(2)*k{2} + c2(3)*k{3} + c2(4)*k{4});
            catch ME
                fprintf('Error evaluating RHS function at t=%g with dt=%g', t, dt);
                rethrow(ME);
            end

            % Estimate error
            error_est = w_next_3 - w_next_2;
            
            % Calculate error norm (relative to solution magnitude)
            % Add small epsilon to avoid division by zero if w or w_next_3 is zero
            epsilon = 1e-10; 
            scale = tol * max(norm(w, Inf), norm(w_next_3, Inf)) + epsilon;
            error_norm = norm(error_est, Inf) / scale; 
            % Alternative Norms:
            % scale = tol * (abs(w) + abs(w_next_3))/2 + epsilon; % Weighted component-wise
            % error_norm = norm(error_est ./ scale, Inf);
            % error_norm = norm(error_est) / (tol * (1 + norm(w))); % RMS norm used before

            % Calculate optimal step size factor
            if error_norm == 0
                factor = max_factor;
            else
                factor = safety * (error_norm ^ (-1 / order));
            end
            % Limit the factor
            factor = min(max_factor, max(min_factor, factor));

            % Check if step is accepted
            if error_norm <= 1.0 || dt <= dt_min % Accept if error is okay OR dt is already minimum
                % Step accepted
                t_prev = t;
                w_prev = w;
                t = t + dt;
                w = w_next_3; % Use the higher-order result
                dt_new = dt * factor; % Suggested dt for the *next* step
                dt = min(dt_new, t_end - t); % Use suggested dt, but don't overshoot
                accepted = true;
                
                % --- Progress Reporting --- %
                if t >= next_report_time
                    perc_done = round(100 * (t - t_start) / (t_end - t_start));
                    % Avoid printing the same percentage twice if steps are small
                    if perc_done > last_report_perc || perc_done == 0 
                       fprintf('\n  t = %.3f s (%.1f%%), dt = %.3e s', t, perc_done, dt_new); % Match RK4 format
                       next_report_time = next_report_time + report_interval;
                       % Ensure next report time doesn't slightly exceed t_end due to float arithmetic
                       if next_report_time > t_end
                           next_report_time = t_end + 1; % Effectively disable further reports
                       end
                       last_report_perc = perc_done;
                    end
                end

                % Store results if past the next plot time
                while t >= next_plot_time && next_plot_time <= t_end
                    out_idx = out_idx + 1;
                    % Grow arrays if needed
                    if out_idx > size(t_out, 1)
                        t_out = [t_out; zeros(n_out_est, 1)]; %#ok<AGROW>
                        sol_out = [sol_out; zeros(n_out_est, n_vars)]; %#ok<AGROW>
                    end
                    % Interpolate if necessary (simple linear interp)
                    theta = (next_plot_time - t_prev) / (t - t_prev);
                    w_interp = w_prev + theta * (w - w_prev);
                    t_out(out_idx) = next_plot_time;
                    sol_out(out_idx, :) = w_interp';
                    next_plot_time = next_plot_time + dt_plot;
                end

            else
                % Step rejected, reduce step size for retry
                dt = dt * factor;
                % Check for stagnation or extremely small steps
                if dt < dt_min 
                    warning('Step size reduced below minimum (%g) after rejection at t=%g.', dt_min, t);
                    % We might accept the current step with the minimal dt next time
                    % or force dt = dt_min here if needed
                end
            end % end accept/reject check
        end % end while ~accepted (inner loop)
    end % end while t < t_end (main loop)

    % --- Finalization --- %
    
    % Store final state if not already stored exactly at t_end
    if abs(t_out(out_idx) - t_end) > 1e-10
       out_idx = out_idx + 1; 
       if out_idx > size(t_out, 1)
            t_out = [t_out; zeros(1, 1)]; %#ok<AGROW>
            sol_out = [sol_out; zeros(1, n_vars)]; %#ok<AGROW>
        end
       t_out(out_idx) = t; % t should be t_end here
       sol_out(out_idx, :) = w';
    end

    % Trim unused preallocated space
    t_out = t_out(1:out_idx);
    sol_out = sol_out(1:out_idx, :);

    % Set statistics
    stats.nsteps = step;
    stats.nfailed = step - out_idx; % Rough estimate, needs refinement if rejects are counted differently
    stats.integration_method = 'Bogacki-Shampine (BS3(2)) Adaptive';

    % Final progress report at 100%
    fprintf('\n'); % Newline before final message
 
    if t < t_end
        warning('Integration terminated early at t = %g (max_steps = %d reached).', t, max_steps);
    else
        fprintf('Integration finished at t = %.3f s after %d steps.', t, step); % Match RK4 format
    end
end