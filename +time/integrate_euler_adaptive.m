function [t_out, w_out] = integrate_euler_adaptive(rhs_func, t_span, w0, cfg)
    % Adaptive time-step forward Euler integrator for 1D NSW.
    % Ensures output is stored exactly at times t0, t0+dt_plot, t0+2*dt_plot, ..., tf

    % Imports and Setup
    import bc.* % Ensure BCs are visible if needed by RHS or CFL calc

    t0 = t_span(1);
    tf = t_span(end);
    fprintf('[DEBUG] Initial tf assigned from t_span(end) = %.6f\n', tf); % DEBUG
    dt_plot = cfg.time.dt_plot; % Get the plot interval
    if ~isfield(cfg.time, 'dt_plot') || isempty(dt_plot) || dt_plot <= 0
        warning('cfg.time.dt_plot not set or invalid, using default of 0.1s');
        dt_plot = 0.1;
    end
    if t0 > tf
        error('Start time t0 must be less than or equal to end time tf.');
    end

    fprintf('Starting adaptive Euler integration...\n  Plotting interval dt_plot = %.4f s\n', dt_plot);

    t = t0;
    w = w0;
    num_vars = length(w0);

    % --- Preallocate output arrays ---
    % Estimate number of outputs: initial + intermediate plot steps + final
    % Add a small buffer to avoid re-allocation in edge cases
    num_plot_steps = floor((tf - t0) / dt_plot);
    max_outputs = num_plot_steps + 5; % Initial + plot steps + final + buffer
    t_out = zeros(1, max_outputs);
    w_out = zeros(num_vars, max_outputs);

    % --- Store initial condition ---
    output_count = 1;
    t_out(output_count) = t;
    w_out(:, output_count) = w;
    fprintf('  Stored initial output at t = %.4f s\n', t);

    % --- Calculate the *next* target plot time ---
    % Find the first multiple of dt_plot strictly greater than t0, or tf if closer.
    if abs(t0) < 1e-12 % Handle t0=0 case
        t_plot_next = dt_plot;
    else
        % Find the next multiple of dt_plot >= t0
        current_plot_multiple = floor(t0 / dt_plot);
        t_plot_target = (current_plot_multiple + 1) * dt_plot;
        % If t0 is *exactly* on a plot time, target the *next* one
        if abs(t0 - current_plot_multiple * dt_plot) < 1e-9 * dt_plot
             t_plot_next = t_plot_target;
        else 
            % If t0 is between plot times, target the next one
            t_plot_next = ceil(t0 / dt_plot) * dt_plot; 
        end
    end
    t_plot_next = min(t_plot_next, tf); % Don't go past tf
    
    % Add final time to the list of target times implicitly
    final_time_target = tf;
    fprintf('[DEBUG] Initial final_time_target assigned from tf = %.6f\n', final_time_target); % DEBUG

    step = 0;
    TOL = 1e-9 * max(dt_plot, tf-t0); % Tolerance for floating point comparisons relative to timescale
    if tf == t0 % Handle zero duration case
        fprintf('Integration finished at t = %.4f s after %d steps (t0=tf).\n', t_out(end), step);
        t_out = t_out(1:output_count);
        w_out = w_out(:, 1:output_count);
        return;
    end

    % --- Time Stepping Loop ---
    while t < final_time_target - TOL % Loop until very close to tf
        % Restore original fprintf argument order
        fprintf('[DEBUG] Loop Start: t = %.6f | Target (final_time_target) = %.6f | Next Plot (t_plot_next) = %.6f\n', t, final_time_target, t_plot_next); % DEBUG - Original order
        step = step + 1;

        % --- Calculate adaptive timestep based on current state w ---
        dt_adaptive = core.utils.calculate_dt_cfl(w, cfg); % Use fully qualified name
        if dt_adaptive <= 0
            error('Non-positive adaptive timestep calculated at t=%.4f. Check CFL calculation or simulation state.', t);
        end

        % --- Determine the actual timestep dt to take ---
        % Time remaining until the *next mandatory output* time (either plot or final)
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
             fprintf('[DEBUG] Small dt detected! t=%.6f, dt_adaptive=%.4e, dt_to_next_output=%.4e, dt_chosen=%.4e\n', t, dt_adaptive, dt_to_next_output, dt); % DEBUG
             if t >= final_time_target - TOL
                  fprintf('Reached final time tf=%.4f\n', final_time_target);
                  break; % Exit loop naturally
             else
                  warning('Timestep became extremely small (%.2e) at t=%.4f before reaching tf. Stepping directly to next output time.', dt, t);
                  dt = dt_to_next_output;
                  hit_output_time = true;
                  if dt < TOL * 1e-3 % If even that is too small, break
                       fprintf('Next output time is also too close. Ending simulation early at t=%.4f\n', t);
                       break;
                  end
             end
        end

        % --- Perform Euler Step ---
        dwdt = feval(rhs_func, t, w, cfg); % Calculate derivative
        w_new = w + dt * dwdt;             % Euler update
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
            if output_count > size(t_out, 2)
                % Dynamically resize output arrays if preallocation was insufficient
                 increase_size = max(10, ceil(0.2 * size(t_out, 2))); % Increase by 20% or 10, whichever is larger
                 t_out = [t_out, zeros(1, increase_size)]; 
                 w_out = [w_out, zeros(num_vars, increase_size)];
                 fprintf('Warning: Resized output arrays at step %d (t=%.4f)\n', step, t);
            end
            t_out(output_count) = t; % Store the actual time reached
            w_out(:, output_count) = w;
            fprintf('  Stored output at t = %.4f s (Step %d, dt = %.4e s)\n', t, step, dt);

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
         if output_count > size(t_out, 2)
             t_out = [t_out, zeros(1,1)]; % Add just one slot
             w_out = [w_out, zeros(num_vars, 1)];
         end
         t_out(output_count) = final_time_target; % Store exactly tf
         w_out(:, output_count) = w; % Store the final computed state
         fprintf('  Stored final output state explicitly at t = %.4f s\n', final_time_target);
    elseif abs(t_out(output_count) - final_time_target) <= TOL
         % If the last stored time IS tf, ensure the stored state is the *current* state
         w_out(:, output_count) = w;
    end

    % Trim unused preallocated space
    t_out = t_out(1:output_count);
    w_out = w_out(:, 1:output_count)'; % Transpose to M_out x (2*N)

    fprintf('[DEBUG] integrate_euler_adaptive: Returning w_out size = [%d, %d]\n', size(w_out, 1), size(w_out, 2)); % DEBUG
    fprintf('Integration finished at t = %.4f s after %d steps.\n', t_out(end), step);
end