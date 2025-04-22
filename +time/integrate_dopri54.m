%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +time/integrate_dopri54.m
%
% Purpose:
%   Solves a system of ordinary differential equations (ODEs) using the
%   Dormand-Prince 5(4) method (DOPRI54), an adaptive embedded Runge-Kutta
%   method. This method provides 5th-order accuracy with an embedded 4th-order
%   solution for error estimation. The Dormand-Prince method is optimized for
%   efficiency and is the basis for MATLAB's ode45 solver.
%
% Syntax:
%   [sol_out, t_out, stats] = integrate_dopri54(rhs_handle, tspan, w0, cfg)
%
% Inputs:
%   rhs_handle - [function handle] Function handle for the RHS of the ODE system.
%                Expected signature: dwdt = rhs_handle(t, w, cfg)
%   tspan      - [vector] Time points at which the solution should be output.
%                The integration proceeds from tspan(1) to tspan(end).
%   w0         - [column vector] Initial state vector at time tspan(1).
%   cfg        - [struct] Configuration structure containing parameters passed
%                to the rhs_handle and controlling the integrator.
%                Required fields in cfg.time:
%                  adaptive_tol:   Desired tolerance for error control.
%                Optional fields in cfg.time:
%                  dt_init:        Initial time step guess [s].
%                  safety_factor:  Safety factor for step size adjustment (default: 0.9).
%                  min_factor:     Minimum step size reduction factor (default: 0.2).
%                  max_factor:     Maximum step size increase factor (default: 5.0).
%                  max_steps:      Maximum allowed integration steps (default: 1e6).
%                  num_progress_reports: Number of progress updates to print (default: 10).
%
% Outputs:
%   sol_out    - [M x N matrix] Solution matrix, where M is the number of output
%                time points and N is the size of the state vector. Each row
%                sol_out(i,:) is the state vector at time t_out(i).
%   t_out      - [M x 1 vector] Column vector of time points corresponding to
%                the rows in sol_out.
%   stats      - [struct] Structure containing integration statistics:
%                stats.nsteps:   Total number of successful time steps taken.
%                stats.nfailed:  Total number of rejected steps (failed attempts).
%                stats.nfevals:  Total number of function evaluations.
%
% Dependencies:
%   - The RHS function specified by rhs_handle.
%   - Configuration settings within the cfg structure.
%
% References:
%   - Dormand, J. R., & Prince, P. J. (1980). "A family of embedded Runge-Kutta
%     formulae". Journal of Computational and Applied Mathematics, 6(1), 19-26.
%   - Hairer, E., Nørsett, S. P., & Wanner, G. (1993). Solving Ordinary
%     Differential Equations I: Nonstiff Problems. Springer.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   22 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol_out, t_out, stats] = integrate_dopri54(rhs_handle, tspan, w0, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Configuration Parameters and Initialization                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract parameters from the configuration struct and set defaults if missing
    t_start = tspan(1);
    t_end   = tspan(end);
    
    % Tolerance for error control
    if isfield(cfg.time, 'adaptive_tol')
        tol = cfg.time.adaptive_tol;
    else
        tol = 1e-4;
        warning('DOPRI54:DefaultTolerance', 'Using default adaptive_tol = 1e-4');
    end
    
    % Initial time step guess
    if isfield(cfg.time, 'dt_init')
        dt = cfg.time.dt_init;
    else
        dt = (t_end - t_start) / 100; % Default initial step is 1% of total span
        warning('DOPRI54:DefaultInitialStep', 'Using default dt_init = %.2e', dt);
    end
    
    % Safety factor for step size control (prevents stepping too close to stability limit)
    if isfield(cfg.time, 'safety_factor')
        safety = cfg.time.safety_factor;
    else
        safety = 0.9;
    end
    
    % Step size adjustment limits
    if isfield(cfg.time, 'min_factor')
        min_factor = cfg.time.min_factor;
    else
        min_factor = 0.2; % Minimum step size reduction factor
    end
    
    if isfield(cfg.time, 'max_factor')
        max_factor = cfg.time.max_factor;
    else
        max_factor = 5.0; % Maximum step size increase factor
    end
    
    % Maximum number of steps (to prevent infinite loops)
    if isfield(cfg.time, 'max_steps')
        max_steps = cfg.time.max_steps;
    else
        max_steps = 1e6;
    end
    
    % Progress reporting
    if isfield(cfg.time, 'num_progress_reports')
        num_reports = cfg.time.num_progress_reports;
    else
        num_reports = 10;
    end
    report_interval = (t_end - t_start) / num_reports;
    next_report_time = t_start + report_interval;
    
    % Problem dimensionality
    w = w0(:); % Ensure w0 is a column vector
    N = length(w);
    
    % Statistics initialization
    nsteps = 0;
    nfailed = 0;
    nfevals = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dormand-Prince 5(4) Butcher Tableau                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The Butcher tableau defines the coefficients for the method
    % Stage time fractions
    c = [0; 1/5; 3/10; 4/5; 8/9; 1; 1];
    
    % Runge-Kutta matrix (a_ij coefficients)
    a = zeros(7, 6);
    a(2,1) = 1/5;
    a(3,1:2) = [3/40, 9/40];
    a(4,1:3) = [44/45, -56/15, 32/9];
    a(5,1:4) = [19372/6561, -25360/2187, 64448/6561, -212/729];
    a(6,1:5) = [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656];
    a(7,1:6) = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84];
    
    % 5th order solution weights (b_i coefficients)
    b5 = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
    
    % 4th order solution weights for error estimation
    b4 = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    
    % Error estimator coefficients (difference between b5 and b4)
    e = b5 - b4;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preallocate Output Arrays                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine number of output points
    num_out_times = length(tspan);
    
    % Preallocate output arrays
    sol_out = zeros(num_out_times, N);
    t_out = zeros(num_out_times, 1);
    
    % Store initial condition at first output time
    output_idx = 1;
    sol_out(output_idx, :) = w';
    t_out(output_idx) = t_start;
    
    % Initialize next output time
    if num_out_times > 1
        output_idx = output_idx + 1;
        t_next_out = tspan(output_idx);
    else
        t_next_out = t_end + 1; % Ensure no more outputs are recorded
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize Stage Variables                                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preallocate arrays for RK stages
    k = zeros(N, 7);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main Time Stepping Loop                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Starting DOPRI54 integration from t=%.3f to t=%.3f\n', t_start, t_end);
    
    t = t_start;
    while t < t_end
        % Safety check for maximum steps
        if nsteps >= max_steps
            warning('DOPRI54:MaxStepsExceeded', 'Maximum number of steps (%d) exceeded.', max_steps);
            break;
        end
        
        % Adjust dt to hit output times exactly
        if t + dt > t_next_out && t < t_next_out
            dt = t_next_out - t;
        end
        
        % Ensure we don't go beyond t_end
        if t + dt > t_end
            dt = t_end - t;
        end
        
        % If the time step is too small, issue warning and proceed with minimum step
        if dt < 1e-10 * (t_end - t_start)
            warning('DOPRI54:TimeStepTooSmall', 'Time step is too small (dt=%.2e) at t=%.3f.', dt, t);
            dt = 1e-10 * (t_end - t_start);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute RK Stages                                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % First stage (RHS evaluation at current state)
        % For FSAL: if we had a previous successful step, we can reuse the last stage
        % But we don't implement FSAL here for simplicity and to ensure correctness
        k(:, 1) = rhs_handle(t, w, cfg);
        nfevals = nfevals + 1;
        
        % Subsequent stages
        for i = 2:7
            % Compute stage time
            stage_t = t + c(i) * dt;
            
            % Compute stage state using previous stages
            stage_w = w;
            for j = 1:(i-1)
                stage_w = stage_w + dt * a(i, j) * k(:, j);
            end
            
            % Evaluate RHS at this stage
            k(:, i) = rhs_handle(stage_t, stage_w, cfg);
            nfevals = nfevals + 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute Solutions and Error Estimate                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 5th order solution (main solution)
        w5 = w;
        for i = 1:7
            w5 = w5 + dt * b5(i) * k(:, i);
        end
        
        % 4th order solution (for error estimation)
        w4 = w;
        for i = 1:7
            w4 = w4 + dt * b4(i) * k(:, i);
        end
        
        % Compute error estimate
        error_est = w5 - w4;
        
        % Scale error relative to solution magnitude and tolerance
        scale = max(abs(w), abs(w5)) * tol + 1e-15;
        error_ratio = max(abs(error_est) ./ scale);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step Acceptance/Rejection and Step Size Adjustment          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine if step should be accepted based on error estimate
        if error_ratio <= 1.0
            % --- Step accepted ---
            nsteps = nsteps + 1;
            t = t + dt;
            w = w5; % Use the 5th order solution
            
            % Store solution at requested output times
            if abs(t - t_next_out) < 1e-10
                sol_out(output_idx, :) = w';
                t_out(output_idx) = t;
                
                % Determine next output time
                if output_idx < num_out_times
                    output_idx = output_idx + 1;
                    t_next_out = tspan(output_idx);
                else
                    t_next_out = t_end + 1; % Ensure no more outputs are recorded
                end
            end
            
            % Progress reporting
            if t >= next_report_time
                fprintf('  t = %.3f s (%.1f%%), dt = %.3e s\n', t, 100*(t-t_start)/(t_end-t_start), dt);
                next_report_time = min(next_report_time + report_interval, t_end);
            end
        else
            % --- Step rejected ---
            nfailed = nfailed + 1;
        end
        
        % Compute new step size based on error estimate
        % Formula: dt_new = dt * min(max_factor, max(min_factor, safety * (1/error_ratio)^(1/5)))
        % Where 1/5 comes from the 5th order method (p+1 = 5, where p = 4 for the embedded method)
        if error_ratio == 0
            dt_factor = max_factor;
        else
            dt_factor = safety * (1.0 / error_ratio)^(1/5);
            dt_factor = max(min_factor, min(max_factor, dt_factor));
        end
        
        dt = dt * dt_factor;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finalize Output                                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Trim output arrays to actual size if not all requested points were reached
    if output_idx < num_out_times
        sol_out = sol_out(1:output_idx, :);
        t_out = t_out(1:output_idx);
        warning('DOPRI54:IncompleteSolution', 'Integration ended before reaching all requested output times.');
    end
    
    % Populate statistics structure
    stats = struct('nsteps', nsteps, 'nfailed', nfailed, 'nfevals', nfevals);
    
    fprintf('Integration finished at t = %.3f s after %d steps (%d rejected).\n', t, nsteps, nfailed);
    fprintf('Total function evaluations: %d\n', nfevals);

end