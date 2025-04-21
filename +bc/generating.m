%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +bc/generating.m
%
% Purpose:
%   Implements a wave-generating boundary condition for the 1DWaveTank code.
%   This boundary injects a prescribed time-dependent water surface elevation
%   (typically a sine wave) at the left or right boundary, and sets the
%   corresponding discharge using Riemann invariants (characteristic theory).
%
% Syntax:
%   w_padded = generating(w_padded, t, side, cfg, num_ghost_cells)
%
% Inputs:
%   w_padded        - [N+2*num_ghost_cells x 2] array. State vector including ghost cells.
%   t               - [scalar] Current simulation time (seconds).
%   side            - [char] 'left' or 'right'. Which boundary to apply.
%   cfg             - [struct] Configuration structure (simulation parameters).
%   num_ghost_cells - [integer] Number of ghost cells to fill (usually 1).
%
% Outputs:
%   w_padded        - [N+2*num_ghost_cells x 2] array. State vector with ghost cells filled.
%
% Required cfg fields:
%   cfg.bc.(side).param.a  - Amplitude of boundary wave (default: 0.1 m)
%   cfg.bc.(side).param.T  - Period of boundary wave (default: 2*pi s)
%   cfg.phys.g             - Gravitational acceleration (m/s^2)
%   cfg.param.H0           - Reference still water depth (m)
%
% Dependencies:
%   None (except for correct configuration setup).
%
% References:
%   - Stoker, J.J. (1957). Water Waves: The Mathematical Theory with Applications.
%   - LeVeque, R.J. (2002). Finite Volume Methods for Hyperbolic Problems.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w_padded = generating(w_padded, t, side, cfg, num_ghost_cells)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input validation & warnings %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if num_ghost_cells ~= 1
        warning('Generating BC typically implemented for 1 ghost cell. Using only the outermost.');
    end

    g = cfg.phys.g;           % [m/s^2] Gravitational acceleration
    N = cfg.mesh.N;           % Number of interior cells

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Select side and extract BC parameters/indices %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(side,'left')
        params = cfg.bc.left.param;
        interior_cell_idx = num_ghost_cells + 1;      % First interior cell
        ghost_cell_idx = num_ghost_cells;             % Outermost ghost cell (left)
    elseif strcmp(side,'right')
        params = cfg.bc.right.param;
        interior_cell_idx = N + num_ghost_cells;      % Last interior cell
        ghost_cell_idx = N + num_ghost_cells + 1;     % Outermost ghost cell (right)
    else
        error('Invalid side specified. Use ''left'' or ''right''.');
    end

    % Extract amplitude and period for boundary wave
    a = params.a;         % [m] Amplitude of boundary wave
    T = params.T;         % [s] Period of boundary wave
    H0 = cfg.param.H0;    % [m] Reference still water depth

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BoundaryValue: Prescribe water depth at the boundary    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function can be modified for other waveforms.
    function H_boundary = BoundaryValue(time)
        % Computes the target water depth at the boundary at a given time.
        omega = 2*pi/T;                       % [rad/s] Angular frequency
        H_boundary = H0 + a * sin(omega * time); % [m] Water depth
        H_boundary = max(H_boundary, 1e-6);   % Avoid negative depth (dry state)
    end
    % --- End of BoundaryValue ---

    % Set the prescribed water depth in the ghost cell
    H_ghost = BoundaryValue(t);
    w_padded(ghost_cell_idx, 1) = H_ghost;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute discharge (HU) in the ghost cell using Riemann invariants    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The goal is to inject a wave with prescribed H while maintaining
    % compatibility with outgoing characteristics (Riemann invariants).
    %
    % For the shallow water equations, the two Riemann invariants are:
    %   C- = u - 2*c   (left-going)
    %   C+ = u + 2*c   (right-going)
    %
    % For the left boundary (inflow):
    %   - C- is taken from the interior cell (information leaving domain)
    %   - C+ is set by the prescribed H_ghost (information entering domain)
    %
    % For the right boundary (inflow):
    %   - C+ is taken from the interior cell (information leaving domain)
    %   - C- is set by the prescribed H_ghost (information entering domain)
    %
    % This ensures correct wave injection and minimizes spurious reflections.

    H_interior = w_padded(interior_cell_idx, 1);      % [m] Water depth in interior cell
    HU_interior = w_padded(interior_cell_idx, 2);     % [m^2/s] Discharge in interior cell
    if H_interior > 1e-6
        U_interior = HU_interior / H_interior;         % [m/s] Velocity in interior cell
        c_interior = sqrt(g * H_interior);            % [m/s] Wave speed in interior cell
    else
        U_interior = 0;
        c_interior = 0;
    end

    c_ghost = sqrt(g * H_ghost);                      % [m/s] Wave speed in ghost cell

    if strcmp(side, 'left')
        % --- LEFT BOUNDARY (wave enters from left) ---
        % Use outgoing C- from interior, incoming C+ from boundary
        Riemann_minus = U_interior - 2 * c_interior;              % Outgoing charac.
        U_ghost = Riemann_minus + 2 * c_ghost;                    % Solve for U_ghost
        U_ghost = max(U_ghost, 0);                               % Enforce inflow (no outflow)
    elseif strcmp(side, 'right')
        % --- RIGHT BOUNDARY (wave enters from right) ---
        % Use outgoing C+ from interior, incoming C- from boundary
        Riemann_plus = U_interior + 2 * c_interior;               % Outgoing charac.
        U_ghost = Riemann_plus - 2 * c_ghost;                     % Solve for U_ghost
        U_ghost = min(U_ghost, 0);                               % Enforce inflow (no outflow)
    end

    % Set discharge in the ghost cell
    HU_ghost = H_ghost * U_ghost;                                 % [m^2/s] Discharge
    w_padded(ghost_cell_idx, 2) = HU_ghost;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extrapolate to other ghost cells if num_ghost_cells > 1 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if num_ghost_cells > 1
        if strcmp(side,'left')
             for i = 1:(num_ghost_cells-1)
                  w_padded(i,:) = w_padded(ghost_cell_idx,:);
             end
        elseif strcmp(side,'right')
             for i = (ghost_cell_idx+1):(N + 2*num_ghost_cells)
                  w_padded(i,:) = w_padded(ghost_cell_idx,:);
             end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End of generating.m       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end