%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +core/+utils/calculate_dt_cfl.m
%
% Purpose:
%   Computes the maximum stable time step for the 1D Non-Linear Shallow Water
%   equations using the Courant-Friedrichs-Lewy (CFL) condition. Ensures
%   numerical stability of explicit time-stepping schemes by limiting the
%   time step according to the fastest wave speed in the domain.
%
% Syntax:
%   dt = calculate_dt_cfl(w, cfg)
%
% Inputs:
%   w   - [2N x 1, double] State vector at the current time, containing water depth (H)
%         and discharge (HU) for all N cells: w = [H; HU].
%   cfg - [struct] Configuration structure. Required fields:
%         cfg.mesh.N:   [integer] Number of cells.
%         cfg.time.cfl: [double] Desired CFL number (typically <= 1).
%         cfg.mesh.dx:  [double] Spatial cell width [m].
%         cfg.phys.g:   [double] Acceleration due to gravity [m/s^2].
%
% Outputs:
%   dt  - [double] Calculated maximum stable time step according to the CFL condition [s].
%
% Dependencies:
%   None (utility function, but expects correct cfg structure).
%
% References:
%   - LeVeque, R.J. (2002). Finite Volume Methods for Hyperbolic Problems.
%   - Courant, R., Friedrichs, K., & Lewy, H. (1928). "On the Partial Difference Equations of Mathematical Physics."
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dt = calculate_dt_cfl(w, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Checks and Parameter Extraction                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    g = cfg.phys.g;     % [m/s^2] Acceleration due to gravity
    dx = cfg.mesh.dx;   % [m] Spatial step size
    cfl = cfg.time.cfl; % [unitless] CFL number

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract H and HU from state vector w                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H = w(1:N);         % [m] Water depth (first N elements)
    HU = w(N+1:2*N);    % [m^2/s] Discharge (next N elements)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Time Step by CFL Condition                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The CFL condition for explicit FV schemes is:
    %   dt <= cfl * dx / max(|U| + c)
    % where U = HU/H (velocity), c = sqrt(g*H) (wave speed)
    %
    % To avoid instability in dry or nearly dry cells, a minimum threshold
    % is enforced for H. This prevents division by zero and ensures that
    % the wave speed is always well-defined.
    H_min = 1e-6; % [m] Minimum water depth threshold for stable velocity calculation
    H_floor = max(H, H_min); % Ensure H is above the minimum threshold

    % Calculate velocity U = HU / H, handling potential division by zero
    U = zeros(size(H)); % [m/s] Velocity vector
    wet_indices = H > H_min; % Find indices of wet cells
    U(wet_indices) = HU(wet_indices) ./ H(wet_indices);

    % Calculate the wave speed c = sqrt(g*H)
    wave_speed = sqrt(g * H_floor); % [m/s]

    % The maximum signal speed in each cell is |U| + c
    max_speed = abs(U) + wave_speed; % [m/s]
    max_signal_speed = max(max_speed); % [m/s] Maximum over all cells

    % If all cells are dry (H ~ 0), set dt to a large value (simulation will stop)
    if max_signal_speed < 1e-8
        dt = 1e6; % [s] Effectively disables time stepping (simulation should halt)
        return;
    end

    % CFL time step restriction
    dt = cfl * dx / max_signal_speed; % [s]

end