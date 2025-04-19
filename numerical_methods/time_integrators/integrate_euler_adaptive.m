function [t_out, sol_out] = integrate_euler_adaptive(rhs_handle, tspan, w0, cfg)

    %INTEGRATE_EULER_ADAPTIVE Solves ODE system using adaptive Explicit Euler.
    %   [t_out, sol_out] = INTEGRATE_EULER_ADAPTIVE(rhs_handle, tspan, w0, cfg)
    %   integrates the system dw/dt = rhs_handle(t, w, cfg) from tspan(1) to
    %   tspan(end) using the initial condition w0. The time step dt is adapted
    %   at each step based on the CFL condition.
    %
    %   Inputs:
    %       rhs_handle - Function handle to the RHS function.
    %       tspan      - Vector of desired output times [t0, t1, ..., tEnd].
    %       w0         - Initial state vector (flattened).
    %       cfg        - Configuration structure (must contain mesh.dx, param.g,
    %                    and optionally time.CFL).
    %
    %   Outputs:
    %       t_out      - Vector of time points where solution is computed (matches tspan).
    %       sol_out    - Solution array where sol_out(k,:) is the solution at t_out(k).
    
    % --- Input Checks and Parameter Extraction ---
    if ~isfield(cfg.time, 'CFL')
        cfg.time.CFL = 0.9; % Default CFL number if not provided
        fprintf('Using default CFL number: %.2f\n', cfg.time.CFL);
    end
    
    N = cfg.mesh.N;
    dx = cfg.mesh.dx;
    g = cfg.param.g;
    CFL = cfg.time.CFL;
        
    t0 = tspan(1);
    tEnd = tspan(end);
    num_output_steps = length(tspan);
        
    % Initialize solution array and output time vector
    sol_out = zeros(num_output_steps, length(w0));
    t_out = tspan(:); % Ensure column vector
        
    % Set initial condition
    w_current = w0(:); % Ensure column vector
    t_current = t0;
        
    % Store initial condition if t0 is in tspan
    if abs(t_out(1) - t0) < 1e-10
        sol_out(1, :) = w_current';
        output_idx = 2; % Index for next output storage
    else
        output_idx = 1; % Start storing from the first index
    end
        
    fprintf('Starting adaptive Euler integration...\n');
    integration_step = 0;
        
    % --- Time Stepping Loop ---
    while t_current < tEnd
        
        integration_step = integration_step + 1;
            
        % --- Calculate Adaptive Time Step dt ---
            
        % 1. Estimate max wave speed S = |u| + sqrt(gH)
        w_state = [w_current(1:N), w_current(N+1:2*N)]; % Reshape to N x 2
        H = w_state(:, 1);
        HU = w_state(:, 2);
        U = zeros(N, 1);
        wet_indices = H > 1e-10;
        U(wet_indices) = HU(wet_indices) ./ H(wet_indices);
            
        max_wave_speed = max(abs(U) + sqrt(g * max(H, 0))); % Ensure H >= 0 for sqrt
        if max_wave_speed < 1e-10 % Handle completely dry or static case
            dt_cfl = inf; % Allow large step if no wave speed
        else
            dt_cfl = CFL * dx / max_wave_speed;
        end
            
        % 2. Determine time remaining until next output and end time
        t_remaining_total = tEnd - t_current;
        if output_idx <= num_output_steps
            t_to_next_output = t_out(output_idx) - t_current;
        else
            t_to_next_output = inf; % No more outputs needed
        end
            
        % 3. Choose dt: minimum of CFL limit, time to next output, time to end
        dt = min([dt_cfl, t_remaining_total, t_to_next_output]);
            
        % --- Perform Euler Step ---
        if dt < 1e-14 % Avoid tiny steps if already at tEnd or output time
            if t_to_next_output < 1e-14 % If exactly at output time
                t_current = t_out(output_idx); % Snap to exact time
            else % Likely means dt_cfl is extremely small or t_current is already tEnd
                t_current = tEnd; % Force end
                dt = 0; % Prevent further steps
                 end
            end
    
            if dt > 0
                % Evaluate RHS at current time and state
                dwdt = rhs_handle(t_current, w_current, cfg);
                
                % Update state using Explicit Euler
                w_next = w_current + dt * dwdt;
                
                % Update current state and time
                w_current = w_next;
                t_current = t_current + dt;
            end
    
            % --- Store Output ---
            % Check if we have reached or passed the next output time
            if output_idx <= num_output_steps && abs(t_current - t_out(output_idx)) < 1e-10
                sol_out(output_idx, :) = w_current';
                fprintf('  Stored output at t = %.4f s (Step %d, dt = %.4e s)\n', t_current, integration_step, dt);
                output_idx = output_idx + 1;
            elseif output_idx <= num_output_steps && t_current > t_out(output_idx)
                 warning('Time step overshot output time %.4f. Consider reducing CFL.', t_out(output_idx));
                 % Optional: Interpolate back to the exact output time if needed
                 sol_out(output_idx, :) = w_current'; % Store current state as approximation
                 output_idx = output_idx + 1;
            end
    
            % Safety break for potential infinite loops (e.g., if dt becomes NaN or Inf)
            if ~isfinite(dt) || dt <= 0 && t_current < tEnd
                 warning('Invalid time step dt=%.4e encountered at t=%.4f. Stopping integration.', dt, t_current);
                 break;
            end
            
        end % End while loop
        
        % Ensure the last state at tEnd is stored if tEnd was in tspan
        if abs(t_current - tEnd) < 1e-10 && abs(t_out(end) - tEnd) < 1e-10 && output_idx == num_output_steps + 1
            % Already stored by the loop check
        elseif abs(t_out(end) - tEnd) < 1e-10 && output_idx <= num_output_steps
            % If loop finished slightly before tEnd or didn't store last point
            sol_out(end, :) = w_current';
            fprintf('  Stored final output at t = %.4f s\n', tEnd);
        end
        
        fprintf('Integration finished at t = %.4f s after %d steps.\n', t_current, integration_step);
    end

end