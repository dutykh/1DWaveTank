function [sol_out, t_out, stats] = integrate_ssp2_adaptive(rhs_func, tspan, w0, cfg)

    % INTEGRATE_SSP2_ADAPTIVE Solves ODEs using adaptive SSP(2,2) method.
    %   Integrates the system of differential equations dw/dt = rhs_func(t, w, cfg)
    %   from time tspan(1) to tspan(end) with initial condition w0, using the
    %   second-order Strong Stability Preserving Runge-Kutta (SSP2) method with
    %   an adaptive time step determined by the CFL condition.
    %
    %   SSP(2,2) Scheme:
    %   w^(1) = w^n + dt * F(t^n, w^n, cfg)
    %   w^(n+1) = 0.5 * w^n + 0.5 * (w^(1) + dt * F(t^n + dt, w^(1), cfg))
    %
    %   Reference:
    %       Gottlieb, S., Shu, C.-W., & Tadmor, E. (2001). Strong Stability-Preserving
    %       High-Order Time Discretization Methods. SIAM Review, 43(1), 89-112.
    %
    %   [T_OUT, SOL_OUT, K, DT_HISTORY] = INTEGRATE_SSP2_ADAPTIVE(RHS_FUNC, TSPAN, W0, CFG)
    %   integrates the system. RHS_FUNC is a function handle. TSPAN is a vector
    %   specifying the time points where output is desired [t0, t1, ..., tf].
    %   W0 is the initial condition vector.
    %   CFG is a configuration structure required by RHS_FUNC and for CFL calculation.
    %
    %   Outputs:
    %   T_OUT      - Vector of time points corresponding to the solution points.
    %   SOL_OUT    - Matrix of solution vectors at the time points in T_OUT.
    %                Each row corresponds to a time point.
    %   K          - Total number of time steps taken.
    %   DT_HISTORY - Vector containing the dt value used at each time step.

    % Configuration
    t0 = tspan(1);
    tf = tspan(end);
    dt_plot = cfg.vis.dt_plot; % Time interval for storing/plotting results
    cfl_target = cfg.time.CFL;
    dx = (cfg.domain.xmax - cfg.domain.xmin) / cfg.mesh.N;

    % Initialization
    t = t0;
    w = w0;
    k = 0; % Step counter

    % Preallocate output arrays
    num_plots = length(tspan);
    t_out = zeros(1, num_plots);
    sol_out = zeros(num_plots, length(w0));
    t_out(1) = t;
    sol_out(1,:) = w';
    plot_idx = 2; % Index for next storage/plot time
    t_next_plot = tspan(plot_idx);

    % Preallocate dt history (adjust size if needed for very long simulations)
    max_steps_guess = round((tf - t0) / (cfl_target * dx / sqrt(cfg.phys.g * cfg.param.H0)) * 2); % Rough guess
    max_steps_guess = max(max_steps_guess, 10000); % Ensure a minimum size
    dt_history = zeros(1, max_steps_guess);

    fprintf('Starting adaptive SSP(2,2) integration from t=%.3f to t=%.3f, plotting every %.3f s\n', t0, tf, dt_plot);
    fprintf('Output requested at %d time points.\n', num_plots);

    last_report_time = t0;
    num_reports = 10; % Default number of reports
    if isfield(cfg, 'time') && isfield(cfg.time, 'num_progress_reports') && cfg.time.num_progress_reports > 0
        num_reports = cfg.time.num_progress_reports;
    end
    report_interval = (tf - t0) / num_reports; % Report progress roughly num_reports times

    % --- Time Stepping Loop ---
    while t < tf
        % Calculate adaptive time step using CFL condition (based on current state w)
        % The CFL number and dx are expected inside cfg
        dt = core.utils.calculate_dt_cfl(w, cfg);

        % Ensure dt does not overshoot the next plot time or final time
        if t + dt > t_next_plot
            dt = t_next_plot - t;
        end
        if t + dt > tf
            dt = tf - t;
        end

        % Prevent excessively small dt close to output times
        if dt < 1e-9
           dt = 1e-9; % Minimum step size to avoid stalling
           if t + dt > tf % If even min step overshoots, just finish
               dt = tf - t;
           elseif t + dt > t_next_plot
                dt = t_next_plot - t;
           end
        end

        % --- SSP(2,2) Step ---
        if dt > 0 % Proceed only if dt is positive
            % Stage 1
            F_n = feval(rhs_func, t, w, cfg);
            w1 = w + dt * F_n;

            % Stage 2
            F_1 = feval(rhs_func, t + dt, w1, cfg); % Evaluate RHS at t+dt, w1
            w_new = 0.5 * w + 0.5 * (w1 + dt * F_1);

            % Update solution and time
            w = w_new;
            t = t + dt;
            k = k + 1; % Increment step count

            % Record dt
            if k <= length(dt_history)
                dt_history(k) = dt;
            else
                % Resize dt_history if preallocation was insufficient
                dt_history = [dt_history, zeros(1, max_steps_guess)]; %#ok<AGROW>
                dt_history(k) = dt;
                max_steps_guess = max_steps_guess * 2; % Double the guess for next time
                warning('Resized dt_history, consider increasing initial guess.');
            end
        else
             % If dt became zero or negative (shouldn't happen with checks above)
             warning('Zero or negative dt encountered (%.2e). Breaking loop.', dt);
             break;
        end

        % Store solution at plot intervals
        if abs(t - t_next_plot) < 1e-9 || t >= t_next_plot
            if plot_idx <= num_plots
                 t_out(plot_idx) = t;
                 sol_out(plot_idx,:) = w';
                 plot_idx = plot_idx + 1;
                 if plot_idx <= num_plots
                     t_next_plot = tspan(plot_idx);
                 else
                     t_next_plot = tf + 1; % Ensure we don't try to plot beyond tf
                 end
            end
        end

        % --- Progress Reporting ---
        if report_interval > 0 && t - last_report_time >= report_interval
            fprintf('  t = %.3f s (%.1f%%), dt = %.3e s\n', t, (t/tf)*100, dt);
            last_report_time = t;
        end

        % Safety break for excessive steps
        if k > 1e7
            warning('Excessive steps (%d). Breaking loop.', k);
            break;
        end
    end % End while loop

    % Trim unused parts of output arrays and dt_history
    t_out = t_out(1:plot_idx-1);
    sol_out = sol_out(1:plot_idx-1,:);
    t_out = t_out(:)'; % Ensure row vector
    % dt_history = dt_history(1:k); % No longer returned

    fprintf('Integration finished at t = %.3f s after %d steps.\n', t, k);
    stats = struct('nsteps', k, 'nfevals', k);

end % Function end