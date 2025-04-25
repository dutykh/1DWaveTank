%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +core/+utils/odetpbar.m
%
% Purpose:
%   An Output Function for MATLAB's ODE solvers (e.g., ode45, ode113) to
%   display a text-based progress bar in the command window during integration.
%   This function uses the persistent variables `timer_start`, `t_start_sim`, 
%   `t_end_sim`, and `is_initialized` to store the start time, start of simulation, 
%   end of simulation, and initialization status across calls.
%
% Syntax (used within odeset):
%   options = odeset('OutputFcn', @utils.odetpbar);
%   [t, y] = ode45(odefun, tspan, y0, options);
%
% Inputs (as called by ODE solver):
%   t      - [scalar or vector] Current time point(s) in the integration.
%   y      - [vector or matrix] Solution vector(s) at time(s) t. (Unused)
%   flag   - [char] String indicating the status of the integration:
%            'init': Called before the first integration step.
%            'done': Called after the last integration step.
%            [] (empty): Called after each successful integration step.
%   varargin- Optional arguments passed via odeset.
%
% Outputs:
%   status - [scalar] Return value required by ODE OutputFcn.
%            0: Continue integration.
%            1: Stop integration.
%            (This function always returns 0).
%
% Dependencies:
%   Requires the +core/+utils/textprogressbar.m function.
%
% References:
%   - MATLAB documentation for `odeset` and 'OutputFcn'.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function status = odetpbar(t, y, flag, varargin)
    %ODETPBAR Output function for ODE solvers to display a text progress bar.
    %
    %   Usage:
    %       options = odeset('OutputFcn', @utils.odetpbar);
    %       [T, Y] = ode45(odefun, tspan, y0, options);
    %
    %   Inputs:
    %       t       - Current time (scalar) or time interval [t0, tf] (for 'init').
    %       y       - Solution at time t (column vector) or empty (for 'init'/'done').
    %       flag    - String: 'init', 'done', or empty ''.
    %       varargin- Optional arguments passed via odeset.
    %
    %   Outputs:
    %       status  - Integer: 0 to continue, 1 to stop.
    %
    %   Relies on the helper function: utils.textprogressbar

    % Persistent variables to store timing, simulation span, and initialization status
    persistent timer_start t_start_sim t_end_sim is_initialized;

    if nargin < 3 || isempty(flag)
        % --- Regular Step Update --- 
        
        % Only proceed if the progress bar system has been initialized.
        if isempty(is_initialized) || ~is_initialized
            status = 0; 
            return; % Not initialized yet, do nothing and continue solver.
        end
        
        % Calculate progress percentage based on *simulation time*.
        if ~isempty(t) && ~isempty(t_end_sim) && t_end_sim > t_start_sim
            current_t = t(1); % Use the first element of t for progress calculation
            progress = max(0, min(100, (current_t - t_start_sim) / (t_end_sim - t_start_sim) * 100));
        else
            progress = 0; % Default progress if time values are invalid
        end
        
        % --- Update Progress Bar --- 
        % Only update if initialized and flag is empty (standard step call)
        
        if isempty(is_initialized) || ~is_initialized
            
            status = 0;
            return; % Exit if not initialized 
        end
        
        % Check if it's a standard progress update call (flag is empty)
        if isempty(flag)
            
            utils.textprogressbar(progress);
            % Update elapsed time
            elapsed_time = toc(timer_start);
        else
            % Don't update progress bar for 'init' or 'done' flags here
            
        end
        
        status = 0; % Return 0 to signal the solver to continue.

    else
        % --- Handle 'init' and 'done' flags --- 
        
        switch flag
            case 'init'
                % --- Initialize --- 
                
                sim_span = t; % In 'init', t contains the simulation time span [t0, tf]
                t_start_sim = sim_span(1);
                t_end_sim = sim_span(end);
                timer_start = tic; % Start the wall-clock timer
                % --- Now initialize with the title ---
                utils.textprogressbar('ODE integration: '); % Initialize the text bar display FIRST
                is_initialized = true; % Mark as initialized

            case 'done'
                % --- Finalize --- 
                % Only update/print if the system was actually initialized.
                if ~isempty(is_initialized) && is_initialized
                    utils.textprogressbar(100); % Ensure the bar shows 100%
                    % Print the total elapsed *wall-clock* time.
                    if ~isempty(timer_start)
                        fprintf('\n   Integration time: %.3f s', toc(timer_start)); % Remove trailing \n
                    end
                    % Call textprogressbar again with an empty string to reset its internal state.
                    utils.textprogressbar(''); 
                    fprintf('\n'); % Print the final newline after resetting the bar.
                end
                % Clear persistent variables to reset for the next potential run.
                clear timer_start t_start_sim t_end_sim is_initialized;

            otherwise
                % --- Handle unknown flags --- 
                warning('odetpbar:UnknownFlag', 'Unknown flag ''%s'' passed to odetpbar.', flag);
        end
        status = 0; % Required return value for OutputFcn (0 to continue/finish normally)
    end

end