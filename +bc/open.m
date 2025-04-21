%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +bc/open.m
%
% Purpose:
%   Implements an "open" (zeroth-order extrapolation) boundary condition
%   for the 1DWaveTank code. This BC simply copies the state from the
%   nearest interior cell into the ghost cell(s), allowing waves to exit
%   the domain with minimal reflection. It is sometimes called a "non-reflecting"
%   or "outflow" boundary, but in practice, it is only approximately non-reflecting.
%
% Syntax:
%   w_padded = open(w_padded, t, side, cfg, num_ghost_cells)
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
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w_padded = open(w_padded, t, side, cfg, num_ghost_cells)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Zeroth-order extrapolation (open/outflow BC) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This BC is intended to let waves exit the domain freely by copying the
    % state from the nearest interior cell into the ghost cells. This is a
    % simple, robust choice for outflow boundaries, but may not be perfectly
    % non-reflecting for all wave types (especially for nonlinear or dispersive waves).

    N = cfg.mesh.N; % [integer] Number of interior (physical) cells

    if strcmp(side, 'left')
        %---------------------------------------------
        % LEFT boundary: fill ghost cells 1:num_ghost_cells
        % Each ghost cell copies the state of the first interior cell
        interior_cell_idx = num_ghost_cells + 1; % Index of first interior cell
        for i = 1:num_ghost_cells
            ghost_idx = i;
            w_padded(ghost_idx, :) = w_padded(interior_cell_idx, :);
        end
    elseif strcmp(side, 'right')
        %---------------------------------------------
        % RIGHT boundary: fill ghost cells N+num_ghost_cells+1 : N+2*num_ghost_cells
        % Each ghost cell copies the state of the last interior cell
        interior_cell_idx = N + num_ghost_cells; % Index of last interior cell
        for i = 1:num_ghost_cells
            ghost_idx = interior_cell_idx + i;
            w_padded(ghost_idx, :) = w_padded(interior_cell_idx, :);
        end
    else
        error('Invalid side specified. Use ''left'' or ''right''.');
    end

    % Note: This BC does not account for incoming waves or information from outside
    % the domain. For more accurate non-reflecting boundaries, more sophisticated
    % methods (e.g., characteristic-based, sponge layers) may be used.

end