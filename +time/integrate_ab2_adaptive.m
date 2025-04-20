function [sol_out, t_out, stats] = integrate_ab2_adaptive(rhs_func, t_span, w0, cfg)

    % INTEGRATE_AB2_ADAPTIVE Solves ODE using 2nd-order Adams-Bashforth (AB2)
    % with adaptive time step.
    %
    %   Solves the system of ODEs dw/dt = rhs_func(t, w, cfg) over the time
    %   interval t_span = [t0, tf] with initial condition w0, using the
    %   explicit 2nd-order Adams-Bashforth method with an adaptive time step
    %   dt determined by the CFL condition.
    %
    %   The first step is performed using Forward Euler.
    %
    %   AB2 Formula: w_{n+1} = w_n + dt/2 * (3*f(t_n, w_n) - f(t_{n-1}, w_{n-1}))
    %
    %   Inputs:
    %       rhs_func : Function handle for the RHS of the ODE system.
    %                  Expected signature: f = rhs_func(t, w, cfg)
    %       t_span   : 1xM array defining output times [t0, t1, ..., tf].
    %                  The integration runs from t0 to tf.
    %       w0       : Initial condition vector (at t0).
    %       cfg      : Configuration structure containing parameters needed by
    %                  rhs_func and the time stepper (e.g., cfg.time.CFL).
    %
    %   Outputs:
    %       sol_out  : Solution matrix where each row corresponds to a time point
    %                  in t_out. Size is M x length(w0).
    %       t_out    : Row vector of time points corresponding to the rows of sol_out.
    %                  These are the requested output times from t_span.
    %       stats    : Struct containing statistics:
    %                  - nsteps: Total number of internal time steps taken.
    %                  - nfevals: Total number of RHS function evaluations.
    %
    %   The adaptive time step dt is calculated at each step using:
    %       dt = core.utils.calculate_dt_cfl(w, cfg);
    %
    %   Outputs are stored only at the time points specified in t_span.

    t0 = t_span(1);
    tf = t_span(end);
    t_out_req = t_span; % Requested output times
    num_out_points = length(t_out_req);

    % Preallocate output arrays conservatively (estimate steps)
    dt_est = core.utils.calculate_dt_cfl(w0, cfg); % Initial estimate
    est_steps = ceil((tf - t0) / dt_est) * 2; % Estimate * safety factor
    max_steps = max(est_steps, num_out_points * 10); % Ensure enough space

    sol_out = zeros(num_out_points, length(w0));
    t_out = zeros(1, num_out_points);

    % --- Initialization ---
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

    % --- Time Stepping Loop ---
    while t < tf
        % Check if max steps exceeded (safety break)
        if k >= max_steps
            warning('Maximum number of steps (%d) exceeded. Aborting integration at t=%.3f s.', max_steps, t);
            break;
        end

        % Calculate adaptive time step based on current state w
        dt = core.utils.calculate_dt_cfl(w, cfg);

        % Ensure dt does not step over the final time tf
        if t + dt > tf
            dt = tf - t;
        end

        % Prevent dt from stepping over the next required output time exactly
        % (small adjustment to ensure t_next_plot is reached *after* a step)
        if abs((t + dt) - t_next_plot) < 1e-10 % If step lands exactly on plot time
             % This logic might be unnecessary if comparison is t >= t_next_plot
        elseif t + dt > t_next_plot && t < t_next_plot % If step would cross plot time
            dt = t_next_plot - t; % Step exactly to the plot time
        end

        % --- Adams-Bashforth 2nd Order Step ---
        f_curr = rhs_func(t, w, cfg); % Calculate RHS at current step
        nfevals = nfevals + 1;

        if k == 0
            % First step: Use Forward Euler
            w_next = w + dt * f_curr;
            f_prev = f_curr; % Store f0 for the next (AB2) step
        else
            % Subsequent steps: Use AB2
            w_next = w + (dt / 2) * (3 * f_curr - f_prev);
            f_prev = f_curr; % Store current f for the next step
        end
        % -------------------------------------

        % Update state
        t = t + dt;
        w = w_next;
        k = k + 1;

        % --- Store Output --- 
        % Check if we have reached or passed the next required output time
        while output_count <= num_out_points && t >= t_out_req(output_count) - 1e-9 % Use tolerance for floating point comparison
            if output_count <= size(sol_out, 1)
                sol_out(output_count,:) = w'; % Store current state
                t_out(output_count) = t;      % Store current time
            else
                % Should not happen with preallocation, but as safety:
                warning('Output array size exceeded estimate.');
                sol_out = [sol_out; w'];
                t_out = [t_out, t];
            end
            
            % Update t_next_plot for the *next* required output time
            if output_count < num_out_points
                 t_next_plot = t_out_req(output_count + 1);
            end
            output_count = output_count + 1; % Move to next output slot
        end
        
        % --- Progress Report --- 
        if t >= next_report_time && tf > t0 % Avoid division by zero if t0=tf
            fprintf('  t = %.3f s (%.1f%%), dt = %.3e s\n', t, (t-t0)/(tf-t0)*100, dt);
            next_report_time = next_report_time + report_interval;
            if next_report_time > tf % Don't schedule reports past tf
                next_report_time = tf + 1;
            end
        end

    end % End while loop

    % Trim unused preallocated space (if loop finished before filling arrays)
    output_idx = find(t_out > 0 | abs(t_out - t0) < 1e-12, 1, 'last'); % Find last non-zero time (handle t0 case)
    if isempty(output_idx)
        output_idx = 1; % Ensure at least the initial condition is kept
    end
    t_out = t_out(1:output_idx);
    sol_out = sol_out(1:output_idx, :);
    t_out = t_out(:)'; % Ensure row vector

    % Debug prints (optional)
    % disp(['integrate_ab2_adaptive: t_out size: ', mat2str(size(t_out))]);
    % disp(['integrate_ab2_adaptive: sol_out size: ', mat2str(size(sol_out))]);

    fprintf('Integration finished at t = %.3f s after %d steps.\n', t_out(end), k);
    stats = struct('nsteps', k, 'nfevals', nfevals);

end
