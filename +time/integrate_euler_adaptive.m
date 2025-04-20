function [t_out, w_out, k, dt_history] = integrate_euler_adaptive(rhs_func, t_span, w0, cfg)

    % INTEGRATE_EULER_ADAPTIVE Solves ODE using Forward Euler with adaptive time step.
    %
    %   Integrates the system dw/dt = rhs_func(t, w) from t_span(1) to t_span(end)
    %   using the explicit Forward Euler method. The time step 'dt' is adapted
    %   at each step based on the CFL condition provided by core.utils.calculate_dt_cfl
    %   and ensures that output points specified by cfg.time.dt_plot are met exactly.
    %
    %   Inputs:
    %     rhs_func   - Function handle for the RHS of the ODE: dw/dt = f(t, w).
    %                  It must accept (t, w_flat) and return dwdt_flat.
    %     t_span     - 2-element vector [t_start, t_end] specifying the integration interval.
    %     w0         - Initial condition vector (flat array) at t_start.
    %     cfg        - Configuration structure containing:
    %                  cfg.time.dt_plot: Time interval for storing output.
    %                  cfg.time.cfl: CFL number for adaptive time stepping.
    %                  cfg.mesh: Mesh configuration needed by calculate_dt_cfl.
    %                  cfg.phys: Physical parameters needed by calculate_dt_cfl.
    %
    %   Outputs:
    %     t_out      - Column vector of time points where the solution is stored.
    %                  Includes t_start, t_start+dt_plot, ..., t_end.
    %     w_out      - Solution matrix where each row corresponds to the solution
    %                  vector (flat) at the time points in t_out (M_out x num_vars).
    %     k          - Total number of time steps taken.
    %     dt_history - Vector containing the dt value used at each time step.

    % --- Input Validation and Parameter Extraction ---
    if ~isfield(cfg.time, 'dt_plot') || isempty(cfg.time.dt_plot) || cfg.time.dt_plot <= 0
        warning('cfg.time.dt_plot not set or invalid, using default of 0.1s');
        cfg.time.dt_plot = 0.1;
    end

    t0 = t_span(1);
    tf = t_span(end);
    dt_plot = cfg.time.dt_plot; % Desired interval between output points

    % --- Setup Output Storage ---
    % Estimate the maximum number of output steps needed
    max_output_steps = ceil((tf - t0) / dt_plot) + 2; % +1 for t0, +1 for potential tf rounding
    t_out = zeros(max_output_steps, 1);   % Preallocate time output vector
    num_vars = length(w0);
    w_out = zeros(max_output_steps, num_vars); % Preallocate solution output matrix

    % Store initial condition
    t = t0;
    w = w0(:)'; % Ensure w is a row vector for consistency
    output_count = 1;
    t_out(output_count) = t;
    w_out(output_count, :) = w;

    % Determine the next target time for output
    t_plot_next = t0 + dt_plot;

    % Add final time to the list of target times implicitly
    final_time_target = tf;

    step = 0;
    k = 0; % Step counter
    dt_history = zeros(1, 10000); % Preallocate dt history
    TOL = 1e-9 * max(dt_plot, tf-t0); % Tolerance for floating point comparisons relative to timescale
    fprintf('Starting adaptive Euler integration from t=%.3f to t=%.3f, plotting every %.3f s\n', t0, tf, dt_plot);

    % --- Time Stepping Loop ---
    while t < final_time_target - TOL % Loop until very close to the final target time
        step = step + 1;

        % --- Calculate adaptive timestep based on current state w ---
        % This requires the state vector 'w' in its physical dimensions (e.g., [H; HU])
        % Ensure w is a column vector before passing to CFL calculation
        w = w(:);
        dt_adaptive = core.utils.calculate_dt_cfl(w, cfg); % Use fully qualified name

        % --- Determine timestep 'dt' for the current step ---
        % We need to hit t_plot_next and final_time_target exactly.
        dt_to_next_output = min(t_plot_next - t, final_time_target - t);
        dt_to_next_output = max(dt_to_next_output, 0); % Ensure non-negative

        hit_output_time = false;

        % Decide the step size
        if dt_adaptive >= dt_to_next_output - TOL
            % Adaptive step would overshoot (or land very close to) the next output time.
            % Force the step to land exactly on the output time.
            dt = dt_to_next_output;
            hit_output_time = true;
        else
            % Adaptive step is smaller than time remaining to next output time.
            % Take the adaptive step.
            dt = dt_adaptive;
            % Check if this adaptive step *incidentally* lands on an output time
            if abs((t + dt) - t_plot_next) < TOL || abs((t + dt) - final_time_target) < TOL
                hit_output_time = true;
            end
        end
        
        % Ensure dt is positive and we don't step past tf definitively
        dt = max(dt, 0);
        if t + dt > final_time_target + TOL % Allow tiny overshoot for tolerance
            dt = final_time_target - t;
            hit_output_time = true; % Mark final step to be stored
            dt = max(dt, 0); % Ensure dt isn't negative if t was already tf
        end
        
        % Prevent infinitely small steps if stuck
        if dt < TOL * 1e-3 
             % If dt is tiny because we are *very* close to the next output time,
             % just force dt to hit that time exactly.
             if t >= final_time_target - TOL
                  fprintf('Reached final time tf=%.4f\n', final_time_target);
                  break; % Exit loop naturally if already at the end
             else
                  warning('Timestep became extremely small (%.2e) at t=%.3f before reaching tf. Stepping directly to next output time.', dt, t);
                  dt = dt_to_next_output;
                  hit_output_time = true;
                  if dt < TOL * 1e-3 % If even that is too small, break
                       fprintf('Next output time is also too close. Ending simulation early at t=%.3f\n', t);
                       break;
                  end
             end
        end

        % --- Perform Euler Step ---
        dwdt = feval(rhs_func, t, w, cfg); % Calculate derivative
        w_new = w + dt * dwdt;             % Euler update
        w_new = w_new(:);                  % Ensure column vector
        t_new = t + dt;                    % Update time

        % --- Update State ---
        w = w_new;
        t = t_new;

        % --- Store Output if a mandatory output time was reached ---        
        % Check if current time 't' is very close to the next plot time or final time
        is_plot_time = abs(t - t_plot_next) < TOL;
        is_final_time = abs(t - final_time_target) < TOL;
        
        if is_plot_time || is_final_time 
            output_count = output_count + 1;
            if output_count > max_output_steps
                % Dynamically resize output arrays if preallocation was insufficient
                 increase_size = max(10, ceil(0.2 * max_output_steps)); % Increase by 20% or 10, whichever is larger
                 t_out = [t_out, zeros(1, increase_size)]; 
                 w_out = [w_out, zeros(num_vars, increase_size)];
                 max_output_steps = max_output_steps + increase_size;
                 fprintf('Warning: Resized output arrays at step %d (t=%.3f)\n', step, t);
            end
            t_out(output_count) = t; % Store the actual time reached
            w_out(output_count, :) = w(:)'; % Ensure w is column, then transpose for row assignment
            fprintf('  Stored output at t = %.3f s (Step %d, dt = %.3e s)\n', t, step, dt);

            % --- Update the next plot time target --- 
            if is_plot_time
                % We hit a plot time, find the *next* one strictly greater than current time
                current_plot_multiple = round(t / dt_plot);
                t_plot_next = (current_plot_multiple + 1) * dt_plot;
                 % Handle cases where t might be slightly over due to tolerance
                 while t_plot_next <= t + TOL 
                      t_plot_next = t_plot_next + dt_plot;
                 end
                t_plot_next = min(t_plot_next, final_time_target); % Don't overshoot tf
            end
        end

        % --- Store dt history ---
        k = k + 1;
        dt_history(k) = dt;

        % --- Check for Simulation End (redundant with while condition but safe) ---
        if t >= final_time_target - TOL
            break; % Exit loop
        end
    end % End while loop

    % --- Final State Storage Check ---
    % Ensure the final state at tf is stored if the loop terminated slightly before
    % or exactly at tf, but it wasn't recorded in the last iteration.
    if abs(t_out(output_count) - final_time_target) > TOL && t >= final_time_target - TOL
         % Only store if the *last* stored time wasn't already tf
         output_count = output_count + 1;
         if output_count > max_output_steps
             t_out = [t_out, zeros(1,1)]; % Add just one slot
             w_out = [w_out, zeros(num_vars, 1)];
         end
         t_out(output_count) = final_time_target; % Store exactly tf
         w_out(output_count, :) = w(:)'; % Ensure w is column, then transpose for row assignment
         fprintf('  Stored final output state explicitly at t = %.3f s\n', final_time_target);
    elseif abs(t_out(output_count) - final_time_target) <= TOL
         % If the last stored time IS tf, ensure the stored state is the *current* state
         w_out(output_count, :) = w(:)'; % Ensure w is column, then transpose for row assignment
    end

    % Trim unused preallocated space
    t_out = t_out(1:output_count);
    w_out = w_out(1:output_count, :);
    dt_history = dt_history(1:k);

    fprintf('Integration finished at t = %.3f s after %d steps.\n', t_out(end), k);

end