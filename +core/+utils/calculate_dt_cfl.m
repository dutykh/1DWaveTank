function dt = calculate_dt_cfl(w, cfg)
% CALCULATE_DT_CFL Computes the maximum stable time step based on the CFL condition.
%
%   Calculates the time step 'dt' for the 1D Non-Linear Shallow Water (NSW)
%   equations using the Courant-Friedrichs-Lewy (CFL) condition:
%
%     dt = cfl * min(dx / (|U| + sqrt(g*H)))
%
%   where the minimum is taken over all cells in the domain.
%
%   Inputs:
%     w   - State vector (2N x 1) containing [H; HU] at the current time.
%     cfg - Configuration structure containing:
%           cfg.mesh.N:   Number of cells.
%           cfg.time.cfl: The desired CFL number (typically <= 1).
%           cfg.mesh.dx:  The spatial cell width.
%           cfg.phys.g:   Acceleration due to gravity.
%
%   Outputs:
%     dt  - The calculated maximum stable time step according to the CFL condition.

    % --- Input Checks and Parameter Extraction ---
    if ~isfield(cfg.time, 'cfl') || cfg.time.cfl <= 0
        error('cfg.time.cfl must be a positive value.');
    end
    if ~isfield(cfg.mesh, 'dx') || cfg.mesh.dx <= 0
        error('cfg.mesh.dx must be a positive value.');
    end
    if ~isfield(cfg.phys, 'g') || cfg.phys.g <= 0
        error('cfg.phys.g must be a positive value.');
    end
    if ~isfield(cfg.mesh, 'N') || cfg.mesh.N <= 0
        error('cfg.mesh.N must be a positive integer.');
    end

    N = cfg.mesh.N;
    if length(w) ~= 2*N
        error('Input state vector w must have length 2*N. Got length(w) = %d, N = %d, size(w) = [%d %d]', length(w), N, size(w,1), size(w,2));
    end

    g = cfg.phys.g;     % Acceleration due to gravity
    dx = cfg.mesh.dx;   % Spatial step size
    cfl = cfg.time.cfl; % CFL number

    % --- Extract H and HU from state vector w ---
    H = w(1:N);         % Water depth (first N elements)
    HU = w(N+1:2*N);    % Discharge (next N elements)

    % --- Calculate Time Step --- 
    % Avoid division by zero or issues with dry cells (H near zero)
    H_min = 1e-6; % Minimum water depth threshold for stable velocity calculation
    H_floor = max(H, H_min); % Ensure H is above the minimum threshold

    % Calculate velocity U = HU / H, handling potential division by zero
    U = zeros(size(H)); % Initialize velocity vector
    wet_indices = H > H_min; % Find indices of wet cells
    U(wet_indices) = HU(wet_indices) ./ H(wet_indices);

    % Calculate the wave speed sqrt(g*H)
    wave_speed = sqrt(g * H_floor);

    % Calculate the maximum characteristic speed (|U| + sqrt(g*H)) in each cell
    max_speed = abs(U) + wave_speed;

    % Find the maximum speed across all cells
    max_speed_global = max(max_speed);

    % Calculate the time step based on the CFL condition
    if max_speed_global > 1e-9 % Avoid division by zero if speeds are negligible
        dt = cfl * dx / max_speed_global;
    else
        % If speeds are essentially zero (e.g., lake at rest), return a large/default dt
        % This prevents dt from becoming infinite or NaN.
        % Use dt_plot as a reasonable upper bound, or a fixed large value.
        dt = cfg.time.dt_plot; % Or potentially a large fixed value like 1.0
    end

    % Ensure dt is positive
    dt = max(dt, 0);

end
