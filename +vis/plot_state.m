function fig_handle = plot_state(xc, H, h, t, cfg, fig_handle)

    %PLOT_STATE Visualizes the wave tank state at a specific time.
    %   fig_handle = PLOT_STATE(xc, H, h, t, cfg, fig_handle) plots the
    %   bathymetry and free surface.
    %
    %   Inputs:
    %       xc          - Vector of cell centre coordinates.
    %       H           - Vector of water depth H at time t.
    %       h           - Vector of bathymetry (still water depth).
    %       t           - Current simulation time.
    %       cfg         - Configuration structure (used for BC info, limits).
    %       fig_handle  - Handle to the figure to plot on (optional, creates new if empty).
    %
    %   Outputs:
    %       fig_handle  - Handle to the figure used/created.

    if nargin < 6 || isempty(fig_handle) || ~isgraphics(fig_handle)
        fig_handle = figure; % Create a new figure if needed
    else
        figure(fig_handle); % Bring figure to front
        clf; % Clear the figure for the new frame
    end

    ax = axes('Parent', fig_handle);
    hold(ax, 'on');

    % Calculate bottom elevation zb and free surface elevation eta
    % Assuming h is still water depth, so zb = -h relative to z=0 (SWL)
    zb = -h;
    eta = H - h; % Free surface elevation relative to z=0

    % Plot the water volume using 'fill'
    water_color = [0.5, 0.7, 1.0]; % Light blue
    x_fill = [xc(:); flipud(xc(:))]; % Create closed polygon x-coords
    z_fill = [zb(:); flipud(eta(:))]; % Create closed polygon z-coords
    fill(ax, x_fill, z_fill, water_color, 'EdgeColor', 'none', 'FaceAlpha', 0.7);

    % Plot the free surface line
    plot(ax, xc, eta, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Free Surface');

    % Plot the bottom bathymetry line
    plot(ax, xc, zb, 'k-', 'LineWidth', 2, 'DisplayName', 'Bottom');

    % --- Visualize Boundaries ---
    ymin = min(zb) - 0.1*cfg.param.H0; % Add some margin below bottom
    ymax = max(eta) + 0.2*cfg.param.H0; % Add margin above highest expected surface
    % Adjust if waves are much larger than H0
    ymax = max(ymax, cfg.param.H0 * 0.5); % Ensure some space above SWL even if flat

    bc_line_options = {'LineWidth', 1.5};
    % Left Boundary
    bc_style_L = core.utils.get_bc_style(cfg.bc.left.handle);
    line(ax, [cfg.domain.xmin, cfg.domain.xmin], [ymin, ymax], 'Color', bc_style_L.color, 'LineStyle', bc_style_L.style, bc_line_options{:}, 'DisplayName', ['Left BC: ' func2str(cfg.bc.left.handle)]);

    % Right Boundary
    bc_style_R = core.utils.get_bc_style(cfg.bc.right.handle);
    line(ax, [cfg.domain.xmax, cfg.domain.xmax], [ymin, ymax], 'Color', bc_style_R.color, 'LineStyle', bc_style_R.style, bc_line_options{:}, 'DisplayName', ['Right BC: ' func2str(cfg.bc.right.handle)]);

    % --- Appearance ---
    hold(ax, 'off');
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, 'Position x (m)');
    ylabel(ax, 'Elevation z (m)');
    title(ax, sprintf('Wave Tank State at t = %.3f s', t));
    ylim(ax, [ymin, ymax]);
    xlim(ax, [cfg.domain.xmin, cfg.domain.xmax]);
    legend(ax, 'show', 'Location', 'southeast');

    drawnow; % Update the figure window

end