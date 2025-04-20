% +run_simulation.m

warning('off', 'MATLAB:opengl:SwitchToSoftwareOpenGL');

% Run a 1D wave tank simulation and visualize results.

clear; close all; format longE;
addpath(genpath(pwd)); % Add all subfolders to path


fprintf('Loading simulation configuration...\n');
% --- Call the main configuration function ---
config = cfg.simulation_config();
% --- Configuration loaded ---

fprintf('Starting simulation: %s\n', config.caseName);
results = core.solver(config); % Call the main solver function
fprintf('Simulation finished.\n');

% --- Visualization ---
if isfield(results, 't') && ~isempty(results.t)
    fprintf('Visualizing results...\n');

    % Get bathymetry (still water depth) corresponding to cell centers
    h_bathy = config.bathyHandle(config.mesh.xc, config);

    num_time_steps = length(results.t);
    vis_indices = unique([1, round(num_time_steps/2), num_time_steps]); % Indices for initial, middle, final plots

    fig = []; % Initialize figure handle

    for i = 1:length(vis_indices)
        idx = vis_indices(i);
        current_t = results.t(idx);
        current_H = results.H(idx, :)'; % Ensure column vector
        current_U = results.U(idx, :)'; % Ensure column vector (if needed)

        % Call the visualization function
        fig = vis.plot_state(config.mesh.xc, current_H, h_bathy, current_t, config, fig);

        % Optional: Pause for viewing or saving frames for a movie
        if i < length(vis_indices)
            pause(0.5); % Pause briefly between plots
        end

        % --- Code to save frames for a movie (requires uncommenting and setup) ---
        % movie_dir = fullfile(cfg.outputPath, 'frames');
        % if ~isfolder(movie_dir), mkdir(movie_dir); end
        % frame_filename = fullfile(movie_dir, sprintf('frame_%04d.png', idx));
        % saveas(fig, frame_filename);
        % fprintf('Saved frame: %s\n', frame_filename);
        % --- End Movie Frame Saving ---

    end

    % --- Code to loop through ALL frames for an animation ---
    % disp('Playing animation of all frames...');
    % fig_anim = [];
    % for idx = 1:num_time_steps
    %     current_t = results.t(idx);
    %     current_H = results.H(idx, :)';
    %     fig_anim = vis.plot_state(cfg.mesh.xc, current_H, h_bathy, current_t, cfg, fig_anim);
    %     pause(0.01); % Adjust pause for animation speed
    % end
    % --- End Animation Loop ---

    fprintf('Visualization finished.\n');
else
    fprintf('No results found to visualize.\n');
end

fprintf('Run script finished.\n');