function fig_handle = plot_state(xc, H, h, U, t, cfg, fig_handle, x_limits, y_limits, u_limits)
% PLOT_STATE Visualizes the state of the 1D wave tank at a given time.
% Includes option to plot velocity in a subpanel.
%
% Inputs:
%   xc         - Cell center coordinates (1 x N)
%   H          - Water depth at cell centers (N x 1)
%   h          - Bathymetry depth at cell centers (N x 1)
%   U          - Velocity at cell centers (N x 1)
%   t          - Current time (scalar)
%   cfg        - Configuration structure
%   fig_handle - Handle to the existing figure (optional)
%   x_limits   - X-axis limits [xmin, xmax]
%   y_limits   - Y-axis limits [ymin, ymax] for the main plot
%   u_limits   - Y-axis limits [umin, umax] for velocity subplot (global for all frames)
%
% Outputs:
%   fig_handle - Handle to the figure used for plotting

    validateattributes(xc, {'numeric'}, {'row', 'finite'});
    validateattributes(H, {'numeric'}, {'column', 'finite'});
    validateattributes(h, {'numeric'}, {'column', 'finite'});
    validateattributes(U, {'numeric'}, {'column', 'finite'}); 
    validateattributes(t, {'numeric'}, {'scalar', 'finite'});
    validateattributes(cfg, {'struct'}, {'scalar'});
    validateattributes(x_limits, {'numeric'}, {'row', 'numel', 2, 'finite'});
    validateattributes(y_limits, {'numeric'}, {'row', 'numel', 2, 'finite'});

    plot_velocity = isfield(cfg, 'vis') && isfield(cfg.vis, 'plot_velocity') && cfg.vis.plot_velocity;

    % --- Figure Handling ---
    if isempty(fig_handle) || ~isvalid(fig_handle)
        fig_handle = figure('Color', 'w', 'Position', [651 587 1173 476]);
    else
        figure(fig_handle); % Bring to front
    end
    clf(fig_handle); % Clear figure for new frame
    set(fig_handle, 'color', 'w'); % Ensure background remains white

    % --- Main Plotting Logic ---
    if plot_velocity
        % --- Subplot 1: Water Surface and Bathymetry ---
        ax1 = subplot(2, 1, 1);
    else
        % --- Single Plot: Water Surface and Bathymetry ---
        ax1 = gca; % Use current axes directly
    end

    hold(ax1, 'on');
    % --- Calculate Derived Variables ---
    eta = H - h; % Free surface elevation

    % --- Fancy Plotting (Water Surface/Bathy) ---
    fill_x = [xc(:); flipud(xc(:))];
    fill_y = [eta(:); flipud(-h(:))];
    fill(ax1, fill_x, fill_y, [0.3 0.6 1.0], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); % Light blue water
    plot(ax1, xc, -h, 'Color', [0.4 0.2 0], 'LineWidth', 2.5); % Brownish for bottom
    plot(ax1, xc, eta, 'b-', 'LineWidth', 1.5); % Thinner blue for surface

    % --- Boundary Condition Indicators ---
    y_range = y_limits(2) - y_limits(1);
    marker_y = y_limits(1) + 0.05 * y_range;
    if contains(func2str(cfg.bc.left_handle), 'wall', 'IgnoreCase', true)
        plot(ax1, x_limits(1), marker_y, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 10);
    elseif contains(func2str(cfg.bc.left_handle), 'generating', 'IgnoreCase', true)
        plot(ax1, x_limits(1), marker_y, 'b>', 'MarkerFaceColor', 'b', 'MarkerSize', 12);
    end
    if contains(func2str(cfg.bc.right_handle), 'wall', 'IgnoreCase', true)
        plot(ax1, x_limits(2), marker_y, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 10);
    elseif contains(func2str(cfg.bc.right_handle), 'generating', 'IgnoreCase', true)
        plot(ax1, x_limits(2), marker_y, 'b<', 'MarkerFaceColor', 'b', 'MarkerSize', 12);
    end

    % --- Plot Formatting (Water Surface/Bathy) ---
    ylabel(ax1, '$z$ (m)', 'Interpreter', 'latex', 'FontSize', 16);
    title(ax1, sprintf('Wave Tank State at $t = %.2f$ s', t), 'Interpreter', 'latex', 'FontSize', 18);
    xlim(ax1, x_limits);
    ylim(ax1, y_limits);
    grid(ax1, 'on');
    box(ax1, 'on');
    set(ax1, 'FontSize', 14, 'LineWidth', 1.5, 'GridAlpha', 0.3, 'TickLabelInterpreter', 'latex');
    grid(ax1, 'off');

    % Always place legend outside for water tank plot
    xlabel(ax1, '$x$ (m)', 'Interpreter', 'latex', 'FontSize', 16);
    legend_shown = isfield(cfg, 'vis') && isfield(cfg.vis, 'show_legend') && cfg.vis.show_legend;
    % --- Subplot 2: Velocity (if requested) ---
    if plot_velocity
        ax2 = subplot(2, 1, 2);
    end
    if legend_shown
        lgd = legend(ax1, {'Water', 'Bathymetry', 'Free surface'}, 'Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 13, 'Box','off');
        % Set figure position for consistent output
        set(gcf, 'Position', [461 444 1331 714]);
        % Use normalized units for perfect alignment
        set(ax1, 'Units', 'normalized');
        if plot_velocity
            set(ax2, 'Units', 'normalized');
        end
        % Increase right margin for legend
        left = 0.08; right = 0.20; top = 0.06; bottom = 0.10; vgap = 0.04;
        width = 1 - left - right;
        if plot_velocity
            total_height = 1 - top - bottom - vgap;
            ax1_height = total_height * (2/3); % wave tank: 2/3
            ax2_height = total_height * (1/3); % velocity: 1/3
            % Top panel (wave tank)
            set(ax1, 'Position', [left, bottom + ax2_height + vgap, width, ax1_height]);
            % Bottom panel (velocity)
            set(ax2, 'Position', [left, bottom, width, ax2_height]);
        else
            height = 1 - top - bottom;
            set(ax1, 'Position', [left, bottom, width, height]);
        end
    else
        % Set figure position for consistent output
        set(gcf, 'Position', [461 444 1331 714]);
        % Use normalized units for perfect alignment
        set(ax1, 'Units', 'normalized');
        left = 0.08; right = 0.04; top = 0.06; bottom = 0.10;
        width = 1 - left - right;
        if plot_velocity
            vgap = 0.04;
            total_height = 1 - top - bottom - vgap;
            ax1_height = total_height * (2/3); % wave tank: 2/3
            ax2_height = total_height * (1/3); % velocity: 1/3
            % Top panel (wave tank)
            set(ax1, 'Position', [left, bottom + ax2_height + vgap, width, ax1_height]);
            set(ax2, 'Units', 'normalized');
            % Bottom panel (velocity)
            set(ax2, 'Position', [left, bottom, width, ax2_height]);
        else
            height = 1 - top - bottom;
            set(ax1, 'Position', [left, bottom, width, height]);
        end
    end

    % --- Set Aspect Ratio (Top or Single Plot) ---
    pbaspect(ax1, [4 1 1]);
    hold(ax1, 'off');

    % --- Subplot 2: Velocity (if requested) ---
    if plot_velocity
        % ax2 already created above

        hold(ax2, 'on');

        % Plot Velocity
        plot(ax2, xc, U, '-', 'Color', [0.1 0.4 0.9], 'LineWidth', 1.5); % Blue tone for velocity

        % Plot Formatting (Velocity)
        xlabel(ax2, '$x$ (m)', 'Interpreter', 'latex', 'FontSize', 16);
        ylabel(ax2, '$u$ (m/s)', 'Interpreter', 'latex', 'FontSize', 16);
        xlim(ax2, x_limits);
        ylim(ax2, u_limits);
        grid(ax2, 'off');
        box(ax2, 'on');
        set(ax2, 'FontSize', 14, 'LineWidth', 1.5, 'GridAlpha', 0.3, 'TickLabelInterpreter', 'latex');
        hold(ax2, 'off');

        % Link x-axes
        linkaxes([ax1, ax2], 'x');
    end

    drawnow; % Update the figure window immediately
end