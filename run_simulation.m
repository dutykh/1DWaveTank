%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run_simulation.m
%
% Purpose:
%   Main script to execute the 1D Wave Tank simulation. This script serves as
%   the primary entry point for configuring, running, and visualizing a
%   simulation defined by the settings in `+cfg/simulation_config.m`.
%
% Workflow:
%   1. Sets up the MATLAB environment (clears variables, closes figures, adds paths).
%   2. Loads the simulation configuration structure (`config`) by calling
%      `cfg.simulation_config()`.
%   3. Prints key configuration details to the command window.
%   4. Optionally creates an output directory if results saving is enabled.
%   5. Calls the core solver `core.solver(config)` to run the simulation.
%   6. Measures and prints the CPU time taken for the simulation.
%   7. If results are available, calculates global axis limits for consistent
%      plotting across all time frames.
%   8. Iterates through the simulation results, calling `vis.plot_state` at
%      each output time step to create an animation.
%   9. Prints final execution statistics.
%
% Usage:
%   Simply run this script from the MATLAB command window or editor:
%
%   >> run_simulation
%
%   To change the simulation setup, edit `+cfg/simulation_config.m`.
%
% Inputs:
%   (none) - Reads configuration from `+cfg/simulation_config.m`.
%
% Outputs:
%   (none) - Displays simulation progress, animation (if configured), and statistics
%            to the command window and figure plots.
%          - Creates a `results` directory and saves `.mat` file if configured.
%
% Dependencies:
%   - Requires all package directories (`+cfg`, `+core`, `+vis`, etc.) to be
%     on the MATLAB path (handled by `addpath(genpath(pwd))`).
%   - Requires the functions called within the script (e.g., `cfg.simulation_config`,
%     `core.solver`, `vis.plot_state`, `utils.calculate_dt_cfl`, etc.)
%     to exist and function correctly.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Environment Setup                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optional: Suppress OpenGL warnings if they occur on specific systems.
% warning('off', 'MATLAB:opengl:SwitchToSoftwareOpenGL');

% Start with a clean environment to avoid conflicts from previous runs.
clear; close all; format longE;
rehash toolboxcache; % Force MATLAB to update its cache for new functions/packages

% Add the project root directory and all its subdirectories to the MATLAB path.
% This ensures all package functions (e.g., `+core`, `+cfg`) are accessible.
addpath(genpath(pwd));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Simulation Configuration                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The `simulation_config` function defines different experimental setups
% and returns a complete configuration structure 'config' containing all
% necessary parameters and function handles for the chosen simulation.
fprintf('--- Loading Simulation Configuration ---\n');
config = cfg.simulation_config();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print Key Configuration Details                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display selected parameters to the command window for user verification.

% --- Boundary Conditions --- 
% Print left BC handle and parameters (if any)
left_bc_handle_str = func2str(config.bc.left.handle);
if isfield(config.bc.left, 'param') && ~isempty(fieldnames(config.bc.left.param))
    params = config.bc.left.param;
    param_names = fieldnames(params);
    % Use cellfun to format each parameter as 'name=value'
    param_strs = cellfun(@(name) sprintf('%s=%.3g', name, params.(name)), param_names, 'UniformOutput', false);
    param_str = strjoin(param_strs, ', '); % Join into a single string
    fprintf('  BC Left: %s (Params: %s)\n', left_bc_handle_str, param_str);
else
    fprintf('  BC Left: %s\n', left_bc_handle_str);
end

% Print right BC handle (assuming simple BCs like wall/open often don't need params printed)
right_bc_handle_str = func2str(config.bc.right.handle);
fprintf('  BC Right: %s\n', right_bc_handle_str);

% --- Numerics and Time --- 
fprintf('  Numerical Flux: %s\n', func2str(config.numFlux));
fprintf('  Time Stepper: %s\n', func2str(config.timeStepper));
if isfield(config, 'time') && isfield(config.time, 'cfl')
    fprintf('  Time Span: [%.2f, %.2f] s, CFL: %.2f\n', config.t0, config.tEnd, config.time.cfl);
else
    fprintf('  Time Span: [%.2f, %.2f] s (CFL not applicable)\n', config.t0, config.tEnd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Directory and File Setup (Optional)                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create output directory and define save path if saving results is enabled.
if isfield(config, 'save_results') && config.save_results
    if ~isfield(config, 'outputPath') || isempty(config.outputPath)
        config.outputPath = './results'; % Default output directory if not specified
        warning('Output path not specified in config, using default: %s', config.outputPath);
    end
    if ~isfolder(config.outputPath)
        mkdir(config.outputPath);
        fprintf('Created results directory: %s\n', config.outputPath);
    end
    % Generate a filename incorporating a timestamp for uniqueness.
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('results_%s.mat', timestamp);
    savePath = fullfile(config.outputPath, filename);
    fprintf('Results will be saved to: %s\n', savePath);
    % Note: Actual saving happens after the simulation completes.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the Core Solver                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The `core.solver` function encapsulates the main simulation logic.
% It takes the configuration structure `config` and executes the time stepping
% loop according to the specified numerical methods and parameters.
% It returns a `results` structure containing the simulation output (time vector,
% state variables H, HU, U, etc.) and statistics.

fprintf('--- Running Core Solver ---\n');
tic; % Start timer to measure solver execution time.

results = core.solver(config); % Execute the main simulation function.

cpu_time = toc; % Stop timer and get elapsed time.
fprintf('--- Core Solver Finished (CPU Time: %.3f s) ---\n', cpu_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization / Animation                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if the solver returned valid results before attempting to plot.
if isfield(results, 't') && ~isempty(results.t) && isfield(results, 'H') && ~isempty(results.H)
    
    fprintf('--- Starting Visualization ---\n');
    
    %% Get Bathymetry
    % Calculate bathymetry `h(x)` at cell centers `xc` using the function handle from config.
    % Ensure bathymetry is a row vector for consistent subtraction later.
    h_bathy = config.bathyHandle(config, config.mesh.xc); 
    if iscolumn(h_bathy); h_bathy = h_bathy'; end % Ensure row vector
    
    num_time_steps = length(results.t); % Number of output frames
    fig = []; % Initialize figure handle (plot_state will create if needed)

    %% Compute Global Axis Limits
    % Calculate axis limits *before* the loop to ensure consistent axes across all animation frames.
    % This prevents the axes from rescaling dynamically, which can be distracting.
    
    % --- Y-Limits for Surface/Bathy Plot --- 
    % For the sloping beach case, we want to explicitly set the y-limits
    % to ensure the free surface is at y=0 and the bottom is at y=-1
    if isfield(config, 'experiment_setup') && strcmp(config.experiment_setup, 'sloping_beach')
        % Fixed y-limits for sloping beach case
        y_limits = [-1.5, 0.5]; % Show from y=-1.5 to y=0.5
    else
        % Default calculation for other cases
        eta_all = results.H + h_bathy; % Free surface = Total Water Depth H + Bottom Elevation z_b(x)
        
        y_min_data = min(min(h_bathy(:)), min(eta_all(:))); % Min of bottom and surface
        y_max_data = max(max(h_bathy(:)), max(eta_all(:))); % Max of bottom and surface
        
        range = y_max_data - y_min_data;
        if range < 1e-6; range = 1; end % Avoid zero range for flat cases
        padding = 0.1 * range; % 10% padding
        
        y_limits = [y_min_data - padding, y_max_data + padding];
    end
    
    % --- X-Limits --- 
    x_limits = [min(config.mesh.xc), max(config.mesh.xc)];

    % --- Y-Limits for Velocity Plot --- 
    if isfield(results, 'U') && ~isempty(results.U) % Check if velocity was calculated
        u_min = min(results.U(:));
        u_max = max(results.U(:));
        if u_min == u_max % Handle case of zero or constant velocity
            delta = max(abs(u_min), 1e-2) * 0.1; % 10% margin or at least 0.001
            u_limits = [u_min - delta, u_max + delta];
        else
            u_margin = 0.1 * (u_max - u_min); % 10% margin
            u_limits = [u_min - u_margin, u_max + u_margin];
        end
    else
        u_limits = [-1, 1]; % Default limits if U is not available
    end

    %% Animation Loop
    % Iterate through each output time step stored in the results.
    for idx = 1:num_time_steps
        current_t = results.t(idx);         % [s] Time for the current frame
        current_H = results.H(idx, :)';     % [N x 1, m] Water depth (needs column vector for plot_state)
        if isfield(results, 'U') && ~isempty(results.U)
             current_U = results.U(idx, :)'; % [N x 1, m/s] Velocity (needs column vector)
        else
             current_U = nan(size(current_H)); % Use NaN if velocity is not plotted/available
        end
        current_h_bathy = h_bathy';         % [N x 1, m] Bathymetry (needs column vector)
        
        % Call the plotting function to draw/update the figure.
        % Pass the existing figure handle `fig` to update the same window.
        fig = vis.plot_state(config.mesh.xc, current_H, current_h_bathy, current_U, current_t, config, fig, x_limits, y_limits, u_limits);
        
        drawnow; % Force MATLAB to render the plot immediately.
        pause(0.05); % Pause briefly to control animation speed.
        
        % --- Optional: Save Frames for Movie --- 
        % Uncomment and configure this section to save each frame as an image file.
        % These frames can later be compiled into a movie using tools like ffmpeg.
        % movie_dir = fullfile(config.outputPath, 'frames');
        % if ~isfolder(movie_dir), mkdir(movie_dir); end
        % frame_filename = fullfile(movie_dir, sprintf('frame_%04d.png', idx));
        % saveas(fig, frame_filename); % Or use print() for more control
    end
    fprintf('--- Visualization Finished ---\n');
else
    fprintf('--- Skipping Visualization (No valid results data) ---\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Results (Optional)                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the `results` and `config` structures to a .mat file if requested.
if isfield(config, 'save_results') && config.save_results
    try
        save(savePath, 'results', 'config', '-v7.3'); % Use v7.3 for potentially large files
        fprintf('Results and config saved successfully to: %s\n', savePath);
    catch ME
        warning('run_simulation:SaveError', 'Could not save results to %s. Error: %s', savePath, ME.message);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print Execution Statistics                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--- Simulation Statistics ---\n');
fprintf('  Mesh Cells (N) : %d\n', config.mesh.N);
if isfield(results, 'total_steps') && ~isnan(results.total_steps)
    fprintf('  Total Steps    : %d\n', results.total_steps);
else
    fprintf('  Total Steps    : N/A (Using MATLAB ODE Solver or info missing)\n');
end
fprintf('  CPU Time       : %.3f s\n', cpu_time);
fprintf('-----------------------------\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleanup (Optional)                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove added paths if desired (might be useful in some contexts)
% rmpath(genpath(pwd));

% End of script