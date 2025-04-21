%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +bc/wall.m
%
% Purpose:
%   Implements a solid wall (impermeable, perfectly reflecting) boundary
%   condition for the 1DWaveTank code. This BC enforces zero normal velocity
%   at the wall by mirroring the water height (H) and flipping the sign of
%   the discharge (HU) in the ghost cells. This ensures that the momentum
%   flux at the wall is reversed, as required by the physics of a solid wall.
%
% Syntax:
%   w_padded = wall(w_padded, t, side, cfg, num_ghost_cells)
%
% Inputs:
%   w_padded        - [N+2*num_ghost_cells x 2] array. State vector including ghost cells.
%   t               - [scalar] Current simulation time (seconds). (Unused)
%   side            - [char] 'left' or 'right'. Which boundary to apply.
%   cfg             - [struct] Configuration structure (only cfg.mesh.N used).
%   num_ghost_cells - [integer] Number of ghost cells to fill (>=1).
%
% Outputs:
%   w_padded        - [N+2*num_ghost_cells x 2] array. State vector with ghost cells filled.
%
% Dependencies:
%   Requires cfg.mesh.N to be set correctly.
%
% References:
%   - LeVeque, R.J. (2002). Finite Volume Methods for Hyperbolic Problems.
%   - Toro, E.F. (2001). Shock-Capturing Methods for Free-Surface Shallow Flows.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w_padded = wall(w_padded, t, side, cfg, num_ghost_cells)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solid wall BC: Mirror H, flip sign of HU in ghost cells %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This enforces u=0 at the wall (impermeable boundary),
    % reflecting all incident waves. This is the standard
    % approach for rigid, vertical walls in shallow water models.

    N = cfg.mesh.N; % [integer] Number of interior (physical) cells

    if strcmp(side, 'left')
        %--------------------------------------------------------------
        % LEFT boundary: fill ghost cells 1:num_ghost_cells
        % Each ghost cell i mirrors the state of the i-th interior cell
        interior_idx_start = num_ghost_cells + 1; % First interior cell
        for i = 1:num_ghost_cells
            ghost_idx = num_ghost_cells - i + 1;       % Outermost ghost cell first
            interior_idx = interior_idx_start + i - 1; % Corresponding interior cell
            % Mirror water height (H) and flip discharge (HU)
            w_padded(ghost_idx, 1) =  w_padded(interior_idx, 1);   % [m] H_ghost = H_interior
            w_padded(ghost_idx, 2) = -w_padded(interior_idx, 2);   % [m^2/s] HU_ghost = -HU_interior
        end
    elseif strcmp(side, 'right')
        %--------------------------------------------------------------
        % RIGHT boundary: fill ghost cells N+num_ghost_cells+1 : N+2*num_ghost_cells
        % Each ghost cell i mirrors the state of the i-th interior cell (from the end)
        interior_idx_end = N + num_ghost_cells; % Last interior cell
        for i = 1:num_ghost_cells
            ghost_idx = interior_idx_end + i;           % Outermost ghost cell
            interior_idx = interior_idx_end - i + 1;    % Corresponding interior cell
            % Mirror water height (H) and flip discharge (HU)
            w_padded(ghost_idx, 1) =  w_padded(interior_idx, 1);   % [m] H_ghost = H_interior
            w_padded(ghost_idx, 2) = -w_padded(interior_idx, 2);   % [m^2/s] HU_ghost = -HU_interior
        end
    else
        error('Invalid side specified. Use ''left'' or ''right''.');
    end

    % Note: This BC assumes a perfectly vertical, impermeable wall.
    % For sloping beaches or partially transmitting boundaries, a different
    % approach is required.

end