function fig_handle = plot_state(xc, H, h, t, cfg, fig_handle, x_limits, y_limits)
% PLOT_STATE Visualizes the state of the 1D wave tank at a given time.
%
%   This function plots the bathymetry, free surface elevation, and optionally
%   the velocity profile for the 1D shallow water simulation. It handles figure
%   creation and updates for animation.
%
%   Inputs:
%     xc         - Vector (1xN) of cell center coordinates.
%     H          - Vector (1xN) of water depth H at cell centers at time t.
%     h          - Vector (1xN) of bathymetry depth h(x) at cell centers (positive down).
%     t          - Current simulation time.
%     cfg        - Configuration structure containing visualization and mesh info.
%                  cfg.vis.plot_vars: Cell array of variables to plot ('H', 'U', 'HU').
%                  cfg.mesh.dx: Cell width (for calculating velocity if needed).
%                  cfg.bc: Boundary condition handles (for display).
%     fig_handle - Handle to the figure window. If empty or invalid, creates a new figure.
%     x_limits   - 2-element vector [x_min, x_max] for the x-axis.
%     y_limits   - 2-element vector [y_min, y_max] for the y-axis.
%
%   Outputs:
%     fig_handle - Handle to the figure used for plotting.

    % --- Figure Handling ---
    if isempty(fig_handle) || ~isvalid(fig_handle)
        fig_handle = figure('Color', 'w', 'Position', [651 587 1173 476]);
    else
        figure(fig_handle); % Bring to front
    end
    set(fig_handle, 'color', 'w'); % Set background to white
    if ~isempty(fig_handle) && isvalid(fig_handle)
        clf(fig_handle); % Clear the existing figure for the new frame
        set(fig_handle, 'color', 'w'); % Ensure background remains white
    end
    hold on; % Hold on to plot multiple lines on the same axes

    % --- Calculate Derived Variables ---
    eta = H - h; % Free surface elevation (relative to z=0 datum)
    % Note: Bathymetry h(x) is defined as depth below z=0, so h is typically >= 0.
    % Water depth H is always positive in wet cells.
    % Free surface eta = H - h represents the water level relative to z=0.

    % --- Fancy Plotting ---
    % Fill water region (between eta and -h)
    fill_x = [xc(:); flipud(xc(:))];
    fill_y = [eta(:); flipud(-h(:))];
    fill(fill_x, fill_y, [0.3 0.6 1.0], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); % Light blue water

    % Plot Bathymetry
    plot(xc, -h, 'Color', [0.4 0.2 0], 'LineWidth', 2.5); % Brownish for bottom

    % Plot Free Surface
    plot(xc, eta, 'b-', 'LineWidth', 1.5); % Thinner blue for surface

    % Plot Velocity (Optional, if requested in cfg.vis.plot_vars)
    if any(strcmpi(cfg.vis.plot_vars, 'U')) || any(strcmpi(cfg.vis.plot_vars, 'HU'))
        % Need HU to calculate U. Assume it's available via results if needed,
        % but this function signature only receives H.
        % For now, let's comment out velocity plotting as HU isn't passed directly.
        % If HU were available:
        % U = zeros(size(H));
        % wet_indices = H > 1e-6; % Define wet cells
        % U(wet_indices) = HU(wet_indices) ./ H(wet_indices);
        % yyaxis right; % Use right y-axis for velocity
        % plot(xc, U, 'r--', 'LineWidth', 1);
        % ylabel('Velocity U (m/s)');
        % ax = gca;
        % ax.YAxis(2).Color = 'r';
        % yyaxis left; % Switch back to left axis
    end

    % --- Boundary Condition Indicators (Visual Cue) ---
    % Add markers or lines at boundaries to indicate BC type (optional enhancement)
    y_range = y_limits(2) - y_limits(1);
    marker_y = y_limits(1) + 0.05 * y_range; % Position markers near the bottom
    if contains(func2str(cfg.bc.left_handle), 'wall', 'IgnoreCase', true)
        plot(x_limits(1), marker_y, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 10); % Black square for wall
    elseif contains(func2str(cfg.bc.left_handle), 'generating', 'IgnoreCase', true)
        plot(x_limits(1), marker_y, 'b>', 'MarkerFaceColor', 'b', 'MarkerSize', 12); % Blue triangle for generating
    end
    if contains(func2str(cfg.bc.right_handle), 'wall', 'IgnoreCase', true)
        plot(x_limits(2), marker_y, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 10); % Black square for wall
    elseif contains(func2str(cfg.bc.right_handle), 'generating', 'IgnoreCase', true)
        plot(x_limits(2), marker_y, 'b<', 'MarkerFaceColor', 'b', 'MarkerSize', 12); % Blue triangle for generating
    end

    % --- Plot Formatting ---
    xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$z$ (m)', 'Interpreter', 'latex', 'FontSize', 16);
    title(sprintf('Wave Tank State at $t = %.2f$ s', t), 'Interpreter', 'latex', 'FontSize', 18);
    xlim(x_limits);
    ylim(y_limits);
    grid on;
    box on;
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'GridAlpha', 0.3, 'TickLabelInterpreter', 'latex');

    % Add legend
    legend({'Water', 'Bathymetry', 'Free surface'}, 'Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 13, 'Box','off');

    % Set Aspect Ratio - make the tank look wider
    ax = gca; % Get current axes
    pbaspect(ax, [4 1 1]); % Set plot box aspect ratio (width:height:depth)

    hold off;
    drawnow; % Update the figure window immediately

end