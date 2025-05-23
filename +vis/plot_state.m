%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +vis/plot_state.m
%
% Purpose:
%   Visualizes the state (water surface elevation, bathymetry, and optionally
%   velocity) of the 1D wave tank simulation at a specific time instant.
%   Designed to be called repeatedly within an animation loop (e.g., by
%   `run_simulation.m`) or for plotting single frames of the simulation results.
%   It handles figure creation/updating, subplot management for optional
%   velocity plotting, and consistent styling.
%
% Syntax:
%   fig_handle = plot_state(xc, H, h, U, t, cfg, fig_handle, x_limits, y_limits, u_limits)
%
% Inputs:
%   xc         - [1 x N, double] Row vector of cell center coordinates [m]. These are
%                the x-locations where the simulation state is defined.
%   H          - [N x 1, double] Column vector of total water depth (H) at each
%                cell center `xc` [m].
%   h          - [N x 1, double] Column vector of bottom elevation $z_b(x)$ (relative to the $z=0$ datum) at each cell center `xc` [m].
%   U          - [N x 1, double] Column vector of depth-averaged velocity (U) at each
%                cell center `xc` [m/s].
%   t          - [scalar, double] Current simulation time corresponding to the
%                provided state vectors H, h, U [s].
%   cfg        - [struct] Configuration structure containing visualization and
%                boundary condition settings. Relevant fields:
%                cfg.vis.plot_velocity: [logical] If true, plot velocity in a lower subplot.
%                cfg.vis.show_legend: [logical] If true, display legend(s).
%                cfg.bc.left_handle: [function handle] Used to determine marker for left BC.
%                cfg.bc.right_handle: [function handle] Used to determine marker for right BC.
%   fig_handle - [graphics handle, optional] Handle to an existing figure window.
%                If provided and valid, the function will update this figure.
%                If empty or invalid, a new figure window is created. Default: [].
%   x_limits   - [1 x 2, double] Pre-calculated X-axis limits [xmin, xmax] to be used
%                for all plots, ensuring consistency across animation frames [m].
%   y_limits   - [1 x 2, double] Pre-calculated Y-axis limits [ymin, ymax] for the
%                main plot (eta/bathy), ensuring consistency [m].
%   u_limits   - [1 x 2, double] Pre-calculated Y-axis limits [umin, umax] for the
%                velocity subplot, ensuring consistency [m/s].
%
% Outputs:
%   fig_handle - [graphics handle] Handle to the figure used for plotting. This can
%                be passed back into the function on the next call to update the same figure.
%
% Dependencies:
%   Relies on MATLAB's standard plotting functions (figure, subplot, plot, fill, etc.).
%   Assumes standard structure for the `cfg` input.
%
% References:
%   - Standard MATLAB plotting practices.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fig_handle = plot_state(xc, H, h, U, t, cfg, fig_handle, x_limits, y_limits, u_limits)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Validation                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Ensure inputs have the expected type, size, and properties for robust plotting.
    % Check variable types, dimensions (row/column), and numeric properties.
    validateattributes(xc, {'numeric'}, {'row', 'finite', 'nonempty'}, mfilename, 'xc', 1);
    validateattributes(H, {'numeric'}, {'column', 'finite', 'nonempty'}, mfilename, 'H', 2);
    validateattributes(h, {'numeric'}, {'column', 'finite', 'numel', size(H,1)}, mfilename, 'h', 3); % Must match H size
    validateattributes(U, {'numeric'}, {'column', 'finite', 'numel', size(H,1)}, mfilename, 'U', 4); % Must match H size
    validateattributes(t, {'numeric'}, {'scalar', 'finite', 'nonnegative'}, mfilename, 't', 5);
    validateattributes(cfg, {'struct'}, {'scalar'}, mfilename, 'cfg', 6);
    validateattributes(x_limits, {'numeric'}, {'row', 'numel', 2, 'finite'}, mfilename, 'x_limits', 8);
    validateattributes(y_limits, {'numeric'}, {'row', 'numel', 2, 'finite'}, mfilename, 'y_limits', 9);
    validateattributes(u_limits, {'numeric'}, {'row', 'numel', 2, 'finite'}, mfilename, 'u_limits', 10);

    %% Determine plotting options from configuration structure `cfg`.
    % Use logical flags for clarity in subsequent conditional plotting.
    plot_velocity = isfield(cfg, 'vis') && isfield(cfg.vis, 'plot_velocity') && cfg.vis.plot_velocity;
    legend_shown = isfield(cfg, 'vis') && isfield(cfg.vis, 'show_legend') && cfg.vis.show_legend;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figure and Axes Handling                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Reuse existing figure handle if provided and valid, otherwise create a new one.
    % This enables updating the same figure window during an animation loop.
    if nargin < 7 || isempty(fig_handle) || ~isgraphics(fig_handle, 'figure')
        fig_handle = figure('Name', '1DWaveTank Simulation', 'NumberTitle', 'off', 'Color', 'w');
        % Adjust default figure position for better layout depending on plot configuration.
        if plot_velocity
             set(fig_handle, 'Position', [461 444 1100 650]); % Taller figure needed for two subplots
        else
             set(fig_handle, 'Position', [651 587 1100 450]); % Shorter figure sufficient for one plot
        end
    else
        figure(fig_handle); % Bring the existing figure window to the front if it was obscured.
    end
    clf(fig_handle); % Clear the figure content to draw the new animation frame.
    set(fig_handle, 'color', 'w'); % Ensure background remains white (in case user defaults changed).

    %% Setup Axes based on whether velocity plot is requested.
    % This determines whether one or two subplots are created.
    if plot_velocity
        % Create two subplots stacked vertically.
        ax1 = subplot(2, 1, 1); % Handle for the top plot (water surface/bathymetry)
        ax2 = subplot(2, 1, 2); % Handle for the bottom plot (velocity)
    else
        % Use a single axes object occupying the entire figure area.
        ax1 = gca; % Get handle to the current axes (which is the only one)
        ax2 = [];  % Assign empty handle for consistency in later checks.
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot 1: Water Surface and Bathymetry (Axes ax1)             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold(ax1, 'on'); % Allow multiple plot commands on the first axes object.

    %% Calculate Derived Variables for Plotting
    % Calculate free surface elevation eta(x) = H(x) - h(x).

    %% Calculate Absolute Free Surface Elevation
    % For sloping beach experiment, we need to ensure the free surface is at y=0
    % and the bottom is at y=-1 for the flat region
    if isfield(cfg, 'experiment_setup') && strcmp(cfg.experiment_setup, 'sloping_beach')
        % For sloping beach, we want to ensure:
        % 1. Free surface at y=0 for lake at rest
        % 2. Bottom at y=-1 for the flat region
        
        % Calculate free surface elevation
        % Since h is already the bottom elevation (negative for underwater),
        % and H is the total water depth, the free surface is H + h
        eta_surf = H + h;
        
        % No vertical shift applied - using the actual values
        % to ensure the bottom doesn't move during visualization
    else
        % Standard calculation for other cases
        eta_surf = H + h; % [m] Absolute free surface elevation
    end

    %% Plot Water Volume using `fill`
    % Create vertices for a closed polygon representing the water area.
    % The x-coordinates go forward along the surface (xc) and back along the bottom (flipud(xc)).
    % The y-coordinates go forward along the surface (eta_surf) and back along the bottom (h).
    fill_x = [xc(:); flipud(xc(:))]; % [2N x 1] X-coordinates for the fill polygon
    fill_y = [eta_surf(:); flipud(h(:))]; % [2N x 1] Y-coordinates for the fill polygon (top: eta_surf, bottom: h)
    % Use `fill` for a visually clear and appealing representation of the water body.
    % Set 'DisplayName' for potential inclusion in the legend.
    fill(ax1, fill_x, fill_y, [0.3 0.6 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Water');
    hold(ax1, 'on');
    % Plot bathymetry (bottom)
    plot(ax1, xc, h, 'Color', [0.4 0.2 0], 'LineWidth', 2.5, 'DisplayName', 'Bathymetry'); % Plot bottom elevation h = z_b(x)
    % Plot free surface
    plot(ax1, xc, eta_surf, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Free surface'); % Plot free surface elevation H+h

    %% Add Boundary Condition Indicators
    % Place markers near the bottom of the plot to visually indicate the type of BC applied.
    % This helps quickly identify the simulation setup from the plot.
    % Calculate y-position for boundary markers
    % For sloping beach, place markers at the bottom of the tank
    if isfield(cfg, 'experiment_setup') && strcmp(cfg.experiment_setup, 'sloping_beach')
        % For left boundary, use the bottom elevation at the left edge
        left_bottom = h(1);  % Bottom elevation at left boundary
        right_bottom = h(end);  % Bottom elevation at right boundary
        
        % Use the bottom elevation for marker positions
        marker_y_left = left_bottom;
        marker_y_right = right_bottom;
    else
        % Default behavior for other cases
        y_range = y_limits(2) - y_limits(1);
        marker_y_left = y_limits(1) + 0.05 * y_range; % Position markers 5% up from the bottom axis limit
        marker_y_right = marker_y_left;
    end

    % --- Left Boundary Marker --- 
    if isfield(cfg, 'bc') && isfield(cfg.bc, 'left') && isfield(cfg.bc.left, 'handle') % Check field existence first
        bc_left_name = func2str(cfg.bc.left_handle); % Get the function name as a string
        if contains(bc_left_name, 'wall', 'IgnoreCase', true)
            plot(ax1, x_limits(1), marker_y_left, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'DisplayName', 'Left: Wall');
        elseif contains(bc_left_name, 'generating', 'IgnoreCase', true)
            plot(ax1, x_limits(1), marker_y_left, 'b>', 'MarkerFaceColor', 'b', 'MarkerSize', 12, 'DisplayName', 'Left: Generating');
        elseif contains(bc_left_name, 'open', 'IgnoreCase', true)
             plot(ax1, x_limits(1), marker_y_left, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Left: Open');
        % Add other BC types here using 'elseif contains(...)' if needed
        end
    end
    
    % --- Right Boundary Marker --- 
     if isfield(cfg, 'bc') && isfield(cfg.bc, 'right') && isfield(cfg.bc.right, 'handle') % Check field existence first
        bc_right_name = func2str(cfg.bc.right_handle);
        if contains(bc_right_name, 'wall', 'IgnoreCase', true)
            plot(ax1, x_limits(2), marker_y_right, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'DisplayName', 'Right: Wall');
        elseif contains(bc_right_name, 'generating', 'IgnoreCase', true)
            plot(ax1, x_limits(2), marker_y_right, 'b<', 'MarkerFaceColor', 'b', 'MarkerSize', 12, 'DisplayName', 'Right: Generating');
         elseif contains(bc_right_name, 'open', 'IgnoreCase', true)
             plot(ax1, x_limits(2), marker_y_right, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Right: Open');
        % Add other BC types here using 'elseif contains(...)' if needed
        end
    end

    %% Formatting for Axes ax1 (Surface/Bathy Plot)
    ylabel(ax1, '$z$ (m)', 'Interpreter', 'latex', 'FontSize', 16); % Y-label with units
    title(ax1, sprintf('Wave Tank State at $t = %.2f$ s', t), 'Interpreter', 'latex', 'FontSize', 18); % Title with time
    xlim(ax1, x_limits); % Apply pre-calculated global x-limits for consistency
    ylim(ax1, y_limits); % Apply pre-calculated global y-limits for eta/bathy
    grid(ax1, 'off');     % Add grid lines
    box(ax1, 'on');      % Draw box around the plot
    % Set font size, line width for axes, grid transparency, and use LaTeX for tick labels
    set(ax1, 'FontSize', 14, 'LineWidth', 1.5, 'GridAlpha', 0.3, 'TickLabelInterpreter', 'latex');
    
    % Use at most four integer ticks including endpoints
    x_ticks = unique([x_limits(1), round(linspace(x_limits(1), x_limits(2), 4)), x_limits(2)]);
    set(ax1, 'XTick', x_ticks);
    if ~plot_velocity
        % Only add the x-axis label to this plot if it's the *only* plot.
        xlabel(ax1, '$x$ (m)', 'Interpreter', 'latex', 'FontSize', 16);
    else
        % If velocity subplot exists below, remove x-tick labels from this top plot
        % to avoid redundancy and make the combined plot cleaner.
        set(ax1, 'XTickLabel', '');
    end
    hold(ax1, 'off'); % Release the hold on the axes

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot 2: Velocity (Axes ax2) - Optional                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot velocity profile if requested in the configuration
    if plot_velocity
        plot(ax2, xc, U, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Velocity'); % Plot U vs xc
        hold(ax2, 'on'); % Hold axes for potential future additions (e.g., zero line)
        %% Formatting for Axes ax2 (Velocity Plot)
        ylabel(ax2, '$U$ (m/s)', 'Interpreter', 'latex', 'FontSize', 16); % Y-label with units
        xlabel(ax2, '$x$ (m)', 'Interpreter', 'latex', 'FontSize', 16); % X-label always needed here
        xlim(ax2, x_limits); % Use same x-limits as plot 1 for spatial alignment
        ylim(ax2, u_limits); % Apply pre-calculated global velocity limits for consistency
        grid(ax2, 'off');     % Add grid lines
        box(ax2, 'on');      % Draw box around the plot
        % Set font size, line width, grid transparency, LaTeX ticks
        set(ax2, 'FontSize', 14, 'LineWidth', 1.5, 'GridAlpha', 0.3, 'TickLabelInterpreter', 'latex');
        
        % Use at most four integer ticks including endpoints
        x_ticks = unique([x_limits(1), round(linspace(x_limits(1), x_limits(2), 4)), x_limits(2)]);
        set(ax2, 'XTick', x_ticks);
        hold(ax2, 'off'); % Release the hold on the axes
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Legend Handling                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Display legend only if requested and configured.
    % Placing the legend outside avoids obscuring plot data, especially in animations.
    if legend_shown
        % Create legend for the first axes (eta/bathy) based on DisplayName properties set during plotting.
        lgd1 = legend(ax1, 'show', 'Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 13, 'Box','off');
        
        if plot_velocity
             % Create a separate legend for the second axes (velocity).
             lgd2 = legend(ax2, 'show', 'Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 13, 'Box','off');
             % Adjust subplot positions slightly to prevent legends overlapping the plot area.
             % These normalized units ([left bottom width height]) might need tuning based
             % on figure size and specific legend content length.
             set(ax1, 'Position', [0.10 0.55 0.75 0.38]); % Example position
             set(ax2, 'Position', [0.10 0.10 0.75 0.38]); % Example position
        else
             % If only one plot exists, adjust its position to accommodate its legend.
             set(ax1, 'Position', [0.10 0.15 0.75 0.75]); % Example position
        end
    end

end