function [sol_out, t_out, stats] = integrate_ssp3_adaptive(rhs_func, tspan, w0, cfg)

    % INTEGRATE_SSP3_ADAPTIVE Solves ODEs using adaptive SSP(3,3) method.
    %   Integrates the system of differential equations dw/dt = rhs_func(t, w, cfg)
    %   from time tspan(1) to tspan(end) with initial condition w0, using the
    %   third-order Strong Stability Preserving Runge-Kutta (SSP3) method with
    %   an adaptive time step determined by the CFL condition.
    %
    %   SSP(3,3) Scheme (Heun's 3rd Order):
    %   w^(1)   = w^n + dt * F(t^n, w^n, cfg)
    %   w^(2)   = (3/4)*w^n + (1/4)*( w^(1) + dt*F(t^n + dt, w^(1), cfg) )
    %   w^(n+1) = (1/3)*w^n + (2/3)*( w^(2) + dt*F(t^n + 0.5*dt, w^(2), cfg) )
    %
    %   Reference:
    %       Gottlieb, S., Shu, C.-W., & Tadmor, E. (2001). Strong Stability-Preserving
    %       High-Order Time Discretization Methods. SIAM Review, 43(1), 89-112.
    %
    %   [SOL_OUT, T_OUT, STATS] = INTEGRATE_SSP3_ADAPTIVE(RHS_FUNC, TSPAN, W0, CFG)
    %   integrates the system.
    %   Inputs:
    %       RHS_FUNC - Function handle for the right-hand side (e.g., @(t, w, cfg) core.rhs_...).
    %       TSPAN    - Vector specifying the time points where output is desired [t0, t1, ..., tf].
    %       W0       - Initial condition vector [h0; q0] (column vector).
    %       CFG      - Configuration structure required by RHS_FUNC and for time stepping.
    %                  Must contain: cfg.param.g, cfg.time.CFL, cfg.domain.xmin,
    %                              cfg.domain.xmax, cfg.mesh.N, cfg.vis.dt_plot.
    %                  Optional: cfg.time.num_progress_reports.
    %
    %   Outputs:
    %       SOL_OUT  - Matrix of solution vectors [h; q] at the time points in T_OUT.
    %                  Each row corresponds to a time point.
    %       T_OUT    - Row vector of time points corresponding to the solution points.
    %       STATS    - Structure containing statistics: STATS.nsteps (total steps),
    %                  STATS.nfevals (total RHS evaluations, 3*nsteps for SSP3).
    %
    %   Author: Denys Dutykh
    %   Date:   20 April 2025

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

    fprintf('Starting adaptive SSP(3,3) integration from t=%.3f to t=%.3f, plotting every %.3f s\n', t0, tf, dt_plot);
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

        % --- SSP(3,3) Step ---
        if dt > 0 % Proceed only if dt is positive
            % Stage 1
            F_n = feval(rhs_func, t, w, cfg);
            w1 = w + dt * F_n;

            % Stage 2
            F_1 = feval(rhs_func, t + dt, w1, cfg);
            w2 = (3/4)*w + (1/4)*(w1 + dt * F_1);

            % Stage 3
            F_2 = feval(rhs_func, t + 0.5*dt, w2, cfg);
            w_new = (1/3)*w + (2/3)*(w2 + dt * F_2);

            % Update solution and time
            w = w_new;
            t = t + dt;
            k = k + 1; % Increment step count

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

    % Trim unused parts of output arrays
    t_out = t_out(1:plot_idx-1);
    sol_out = sol_out(1:plot_idx-1,:);
    t_out = t_out(:)'; % Ensure row vector

    fprintf('Integration finished at t = %.3f s after %d steps.\n', t, k);
    stats = struct('nsteps', k, 'nfevals', 3*k); % SSP3 takes 3 RHS evals per step

end % Function end