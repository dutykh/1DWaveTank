% +run_simulation.m

% RUN_SIMULATION Main script to execute the 1D Wave Tank simulation.
%
%   This script performs the following steps:
%   1. Clears workspace and command window, closes figures.
%   2. Adds necessary paths (current directory and subdirectories).
%   3. Loads the simulation configuration using cfg.simulation_config.
%   4. Runs the core solver (core.solver) with the loaded configuration.
%   5. Visualizes the simulation results using vis.plot_state.

warning('off', 'MATLAB:opengl:SwitchToSoftwareOpenGL');

% Start with a clean environment
clear; close all; format longE;

% Add the project root directory and all subdirectories to the MATLAB path
addpath(genpath(pwd));

% --- Load Simulation Configuration ---
% The simulation_config function defines different experimental setups.
% It returns a configuration structure 'config' containing all parameters.
% --- Call the main configuration function ---
config = cfg.simulation_config();
% --- Configuration loaded ---

% --- Run the Solver ---
% The core.solver function takes the configuration 'config' and runs the
% numerical simulation according to the specified model, schemes, and parameters.
% It returns a 'results' structure containing the time vector, H, HU, etc.
results = core.solver(config); % Call the main solver function

% --- Save Results (if requested) ---
if isfield(config, 'save_results') && config.save_results
    % Generate a descriptive filename (timestamp only)
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('results_%s.mat', timestamp); % Simpler filename
    % Construct the full path using the DIRECTORY path from config and the generated filename
    savePath = fullfile(config.outputPath, filename); % config.outputPath is the SUBDIRECTORY path
    fprintf('Results saved to: %s\n', savePath); % Print the FULL path
    save(savePath, 'results', 'config');
end

% --- Visualization ---
if isfield(results, 't') && ~isempty(results.t)
    
    % Get bathymetry (still water depth) corresponding to cell centers
    h_bathy = config.bathyHandle(config.mesh.xc, config);

    num_time_steps = length(results.t);

    fig = []; % Initialize figure handle

        % Compute global axis limits for all frames before animation
    eta_all = results.H - h_bathy'; % Each row: eta at a frame
    eta_min = min(eta_all(:));
    eta_max = max(eta_all(:));
    bathy_min = min(-h_bathy(:));
    bathy_max = max(-h_bathy(:));
    y_min = min(eta_min, bathy_min);
    y_max = max(eta_max, bathy_max);
    margin = 0.1 * (y_max - y_min);
    y_limits = [y_min - margin, y_max + margin];
    x_limits = [min(config.mesh.xc), max(config.mesh.xc)];

    % --- Compute global velocity axis limits for all frames ---
    u_min = min(results.U(:));
    u_max = max(results.U(:));
    u_margin = 0.1 * max(abs([u_min, u_max]));
    u_limits = [u_min - u_margin, u_max + u_margin];

    for idx = 1:num_time_steps
        current_t = results.t(idx);
        current_H = results.H(idx, :)'; % Ensure column vector
        current_U = results.U(idx, :)'; % Ensure column vector (if needed)
        % fprintf('Stored output at t = %.3f s (Step %d, dt = %.3f s)\n', current_t, idx, current_t - results.t(max(1, idx-1))); % Removed to avoid duplicate output
        fig = vis.plot_state(config.mesh.xc', current_H, h_bathy, current_U, current_t, config, fig, x_limits, y_limits, u_limits);
        pause(0.05); % Pause briefly between plots for animation effect
        % --- Code to save frames for a movie (requires uncommenting and setup) ---
        % movie_dir = fullfile(cfg.outputPath, 'frames');
        % if ~isfolder(movie_dir), mkdir(movie_dir); end
        % frame_filename = fullfile(movie_dir, sprintf('frame_%04d.png', idx));
        % saveas(fig, frame_filename);
    end

end