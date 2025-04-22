%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +time/integrate_bogacki_shampine.m
%
% Purpose:
%   Integrates a system of ordinary differential equations (ODEs) using the
%   adaptive Bogacki-Shampine 3(2) embedded Runge-Kutta method. This method
%   provides a 3rd-order accurate solution and uses an embedded 2nd-order
%   solution to estimate the local error and adapt the time step size.
%
% Syntax:
%   [sol_out, t_out, stats] = integrate_bogacki_shampine(rhs_handle, tspan, w0, cfg)
%
% Inputs:
%   rhs_handle - [function handle] Function handle for the RHS of the ODE system.
%                Expected signature: dwdt = rhs_handle(t, w, cfg)
%   tspan      - [1x2 vector] Time interval for integration: [t_start, t_end].
%   w0         - [Nx1 vector] Initial state vector at time t_start.
%   cfg        - [struct] Configuration structure containing parameters passed
%                to the rhs_handle and controlling the integrator.
%                Required fields in cfg.time:
%                  dt_init:        Initial time step guess [s].
%                  adaptive_tol:   Desired relative tolerance for error control.
%                  dt_plot:        Time interval for storing output results [s].
%                Optional fields in cfg.time:
%                  safety_factor:  Safety factor for step size adjustment (default: 0.9).
%                  min_factor:     Minimum step size reduction factor (default: 0.2).
%                  max_factor:     Maximum step size increase factor (default: 5.0).
%                  max_steps:      Maximum allowed integration steps (default: 1e6).
%                  dt_min_factor:  Factor times tspan duration for min dt (default: 1e-12).
%
% Outputs:
%   sol_out    - [M x N matrix] Solution matrix, where M is the number of output
%                time points and N is the size of the state vector. Each row
%                sol_out(i,:) is the state vector at time t_out(i).
%   t_out      - [M x 1 vector] Column vector of time points corresponding to
%                the rows in sol_out.
%   stats      - [struct] Structure containing integration statistics:
%                stats.nsteps: Total number of successful time steps taken.
%                stats.nfailed: Total number of rejected steps (failed attempts).
%
% Dependencies:
%   - The RHS function specified by rhs_handle.
%   - Configuration settings within the cfg structure.
%
% References:
%   - Bogacki, P. and Shampine, L.F. (1989). "A 3(2) Pair of Runge-Kutta Formulas".
%     Applied Mathematics Letters, 2(4), pp. 321-325.
%   - Hairer, E., Nørsett, S. P., & Wanner, G. (1993). Solving Ordinary
%     Differential Equations I: Nonstiff Problems. Springer.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   22 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol_out, t_out, stats] = integrate_bogacki_shampine(rhs_handle, tspan, w0, cfg)

    % --- Configuration Parameters ---
    % Extract parameters from the configuration struct, providing defaults if missing.
    t_start = tspan(1);
    t_end   = tspan(end);
    dt_init = cfg.time.dt_init;             % Initial time step guess [s]
    tol     = cfg.time.adaptive_tol;        % Desired relative tolerance for error control
    dt_plot = cfg.time.dt_plot;             % How often to store results for output [s]

    % Safety factor ensures dt doesn't grow too aggressively after a good step.
    safety  = 0.9;
    if isfield(cfg, 'time') && isfield(cfg.time, 'safety_factor'), safety = cfg.time.safety_factor; end
    % min/max factors prevent dt from changing too drastically in one step.
    min_factor = 0.2; % Minimum step size reduction factor
    if isfield(cfg, 'time') && isfield(cfg.time, 'min_factor'), min_factor = cfg.time.min_factor; end
    max_factor = 5.0; % Maximum step size increase factor
    if isfield(cfg, 'time') && isfield(cfg.time, 'max_factor'), max_factor = cfg.time.max_factor; end
    % Order used in step size adaptation formula (p+1 where p is the lower order, here 2+1=3).
    order   = 3;
    % Maximum steps prevent infinite loops if dt becomes excessively small.
    max_steps = 1e6;
    if isfield(cfg, 'time') && isfield(cfg.time, 'max_steps'), max_steps = cfg.time.max_steps; end
    % Minimum allowed dt prevents issues near singularities or with excessive stiffness.
    dt_min_factor = 1e-12;
    if isfield(cfg, 'time') && isfield(cfg.time, 'dt_min_factor'), dt_min_factor = cfg.time.dt_min_factor; end
    dt_min = dt_min_factor * abs(t_end - t_start);


    % --- Bogacki-Shampine FSAL (First Same As Last) 3(2) Coefficients ---
    % These define the Butcher Tableau for the method.
    % a: Node points (fractions of dt)
    a = [0; 1/2; 3/4; 1];
    % b: Coefficients for calculating intermediate stages (k_i)
    b = [ ...
        0,   0,   0,  0; ... % k1 uses w at t
        1/2, 0,   0,  0; ... % k2 uses k1
        0,   3/4, 0,  0; ... % k3 uses k2
        2/9, 1/3, 4/9, 0 ...  % k4 uses k1, k2, k3
    ];
    % c3: Weights for the 3rd order accurate solution (y_{n+1})
    c3 = [2/9, 1/3, 4/9, 0];
    % c2: Weights for the embedded 2nd order solution (yhat_{n+1}) used for error estimation
    c2 = [7/24, 1/4, 1/3, 1/8]; % Note: k4 is used here, not just k1,k2,k3 as in c3

    % --- Initialization ---
    t = t_start;            % Current simulation time
    w = w0(:);              % Ensure w is a column vector
    dt = dt_init;           % Current time step size
    N = length(w0);         % Size of the state vector
    nsteps = 0;             % Counter for successful steps
    nfailed = 0;            % Counter for rejected steps (failed attempts)

    % Pre-allocate output arrays for efficiency. Estimate size based on dt_plot.
    % Add extra buffer space.
    num_out_est = ceil(abs(t_end - t_start) / dt_plot) + 100;
    sol_out = zeros(num_out_est, N);
    t_out   = zeros(num_out_est, 1);

    % Store initial condition
    store_idx = 1;
    sol_out(store_idx, :) = w';
    t_out(store_idx) = t;
    t_next_plot = t_start + dt_plot; % Time for the next output storage

    % Allocate space for the intermediate stage derivatives (k values)
    k_stages = zeros(N, 4); % Stores k1, k2, k3, k4

    % --- Adaptive Time Stepping Loop ---
    fprintf('Starting adaptive Bogacki-Shampine integration from t=%.3f to t=%.3f\n', t_start, t_end);
    progress_reports = linspace(t_start, t_end, cfg.time.num_progress_reports + 1);
    report_idx = 2; % Start checking from the first interval end

    while t < t_end && nsteps < max_steps
        % Ensure dt doesn't overshoot t_end and is not smaller than dt_min
        dt = max(min(dt, t_end - t), dt_min);

        step_accepted = false; % Flag to track if the current step is accepted

        % --- Inner loop for step rejection/retrial ---
        % This loop repeats a step with smaller dt if the error estimate is too large.
        while ~step_accepted
            if dt < dt_min
                warning('BS3(2): dt reached minimum allowed value dt_min=%.2e at t=%.3f. Stopping.', dt_min, t);
                break; % Exit inner and outer loops
            end

            % --- Calculate Stages (k_i) ---
            % k1 = f(t, w)
            k_stages(:, 1) = rhs_handle(t, w, cfg);
            % k2 = f(t + a(2)*dt, w + dt*b(2,1)*k1)
            k_stages(:, 2) = rhs_handle(t + a(2)*dt, w + dt*b(2,1)*k_stages(:,1), cfg);
            % k3 = f(t + a(3)*dt, w + dt*(b(3,1)*k1 + b(3,2)*k2))
            k_stages(:, 3) = rhs_handle(t + a(3)*dt, w + dt*(b(3,1)*k_stages(:,1) + b(3,2)*k_stages(:,2)), cfg);
            % k4 = f(t + a(4)*dt, w + dt*(b(4,1)*k1 + b(4,2)*k2 + b(4,3)*k3))
            % Note: FSAL property means k4 = f(t+dt, w_next_3), so if step is accepted,
            % k1 for the *next* step will be equal to k4 of the *current* step.
            % This implementation doesn't explicitly use FSAL to simplify logic,
            % but calculates k4 directly based on the formula.
            k_stages(:, 4) = rhs_handle(t + a(4)*dt, w + dt*(b(4,1)*k_stages(:,1) + b(4,2)*k_stages(:,2) + b(4,3)*k_stages(:,3)), cfg);

            % --- Calculate 3rd and 2nd Order Solutions ---
            % w_next_3: Higher order (3rd) estimate - used as the main solution if step accepted
            w_next_3 = w + dt * (k_stages * c3'); % Matrix multiplication k_stages * c3'
            % w_next_2: Lower order (2nd) embedded estimate - used for error calculation
            w_next_2 = w + dt * (k_stages * c2'); % Matrix multiplication k_stages * c2'

            % --- Error Estimation ---
            % Estimate the local truncation error by comparing the 3rd and 2nd order solutions.
            % Scale the error estimate relative to the solution magnitude for tolerance comparison.
            error_est = abs(w_next_3 - w_next_2);
            % Scale factor: Avoid division by zero, use absolute tolerance for small w.
            % Using a common approach: scale = atol + rtol * max(|w|, |w_next_3|)
            % Here, simplified using only rtol (passed as 'tol') and current w.
            scale = abs(w) + 1e-6; % Add small value to avoid zero scale
            error_norm = norm(error_est ./ scale) / sqrt(N); % RMS norm of relative error

            % --- Step Acceptance / Rejection ---
            if error_norm <= tol || dt <= dt_min % Accept step if error is within tolerance OR dt is already minimal
                step_accepted = true;
                t = t + dt; % Advance time
                w = w_next_3; % Update solution
                nsteps = nsteps + 1;

                % Store solution if dt_plot interval is reached
                if t >= t_next_plot || t >= t_end
                    store_idx = store_idx + 1;
                    % Resize output arrays if needed
                    if store_idx > size(sol_out, 1)
                        sol_out = [sol_out; zeros(100, N)]; %#ok<AGROW>
                        t_out   = [t_out; zeros(100, 1)]; %#ok<AGROW>
                    end
                    sol_out(store_idx, :) = w';
                    t_out(store_idx) = t;
                    % Determine next plot time, ensuring it doesn't get stuck if dt is large
                    while t_next_plot <= t
                        t_next_plot = t_next_plot + dt_plot;
                    end
                end

                % Calculate optimal dt for the *next* step based on current error
                if error_norm == 0 % Handle zero error case to avoid division by zero
                    optimal_factor = max_factor;
                else
                    % Standard formula: dt_new = dt * (tol / error_norm)^(1/(p+1))
                    % where p is the lower order (2), so p+1 = 3.
                    optimal_factor = safety * (tol / error_norm)^(1/order);
                end
                % Limit the change in dt
                dt = dt * min(max_factor, max(min_factor, optimal_factor));

            else % Reject step
                nfailed = nfailed + 1;
                step_accepted = false; % Ensure flag is false
                % Reduce dt based on error for the *retry*
                if error_norm == 0 % Avoid issues if error estimate is exactly zero
                     optimal_factor = min_factor;
                else
                    optimal_factor = safety * (tol / error_norm)^(1/order);
                end
                % Use a more conservative reduction for rejected steps, ensure it's less than 1
                dt = dt * min(0.9, max(min_factor, optimal_factor)); % Use min(0.9, ...) to guarantee reduction

                if nsteps > max_steps % Check max steps inside inner loop too
                    warning('BS3(2): Exceeded maximum allowed steps (%d). Stopping.', max_steps);
                    break; % Exit inner loop
                end
            end % End step acceptance check

        end % End inner loop (step retry)

        if dt < dt_min || nsteps >= max_steps % Check exit conditions from inner loop again
             break; % Exit outer loop
        end

        % --- Progress Reporting ---
        if t >= progress_reports(report_idx)
            fprintf('  t = %.3f s (%.1f%%), dt = %.3e s\n', t, (t-t_start)/(t_end-t_start)*100, dt);
            report_idx = report_idx + 1;
        end

    end % End main integration loop (while t < t_end)

    % --- Finalization ---
    fprintf('Integration finished at t = %.3f s after %d steps (%d failed attempts).\n', t, nsteps, nfailed);

    % Trim unused pre-allocated space
    sol_out = sol_out(1:store_idx, :);
    t_out   = t_out(1:store_idx);

    % Store statistics
    stats = struct();
    stats.nsteps = nsteps;
    stats.nfailed = nfailed;

    if t < t_end && nsteps >= max_steps
        warning('BS3(2): Maximum steps reached before reaching t_end.');
    elseif t < t_end && dt <= dt_min
        warning('BS3(2): Time step size reached minimum limit before reaching t_end.');
    end

end % End main function integrate_bogacki_shampine