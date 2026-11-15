%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +bc/open.m
%
% Purpose:
%   Implements an "open" (zeroth-order extrapolation) boundary condition
%   for the 1DWaveTank code. The ghost cells receive the free surface
%   eta = h + z and the discharge of the nearest interior cell, allowing waves
%   to exit the domain with minimal reflection. It is sometimes called a
%   "non-reflecting" or "outflow" boundary, but in practice, it is only
%   approximately non-reflecting. The ghost cells keep the bottom function
%   beyond the end of the tank, and their depth is the extrapolated free
%   surface minus that bottom,
%       h_g = max(0, h_e + z_e - z_g),   (hu)_g = (hu)_e,
%   e being the end cell and g a ghost cell (a dry end cell gives dry ghost
%   cells). Over a flat bottom this is the copy of the state of the end cell;
%   over a sloping one, a lake at rest (eta constant, u = 0) finds the same
%   free surface in the ghost cells as in the tank and stays at rest, and the
%   bottom goes on smoothly beyond the end.
%
% Syntax:
%   w_padded = open(w_padded, t, side, cfg, num_ghost_cells)
%
% Inputs:
%   w_padded        - [N+2*num_ghost_cells x 2] array. State vector including ghost cells.
%   t               - [scalar] Current simulation time (seconds). (Unused)
%   side            - [char] 'left' or 'right'. Which boundary to apply.
%   cfg             - [struct] Configuration structure: cfg.mesh.N, cfg.mesh.dx,
%                     cfg.mesh.xc, cfg.phys.dry_tolerance, cfg.bathyHandle.
%   num_ghost_cells - [integer] Number of ghost cells to fill (>=1).
%
% Outputs:
%   w_padded        - [N+2*num_ghost_cells x 2] array. State vector with ghost cells filled.
%
% Dependencies:
%   Requires cfg.mesh and the bottom function cfg.bathyHandle.
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
    % This BC is intended to let waves exit the domain freely by extending
    % the free surface and the discharge of the nearest interior cell into
    % the ghost cells. This is a simple, robust choice for outflow boundaries,
    % but may not be perfectly non-reflecting for all wave types (especially
    % for nonlinear or dispersive waves).

    N = cfg.mesh.N;                    % [integer] Number of interior (physical) cells
    dx = cfg.mesh.dx;                  % [m] Cell width
    dry_tol = cfg.phys.dry_tolerance;  % [m] Depth below which a cell is dry

    if strcmp(side, 'left')
        %---------------------------------------------
        % LEFT boundary: fill ghost cells 1:num_ghost_cells, outermost first
        interior_cell_idx = num_ghost_cells + 1;          % First interior cell
        ghost_idx = 1:num_ghost_cells;
        x_end = cfg.mesh.xc(1);                            % [m] Centre of the end cell
        x_ghost = x_end - (num_ghost_cells:-1:1) * dx;     % [m] Centres of the ghost cells
    elseif strcmp(side, 'right')
        %---------------------------------------------
        % RIGHT boundary: fill ghost cells N+num_ghost_cells+1 : N+2*num_ghost_cells
        interior_cell_idx = N + num_ghost_cells;           % Last interior cell
        ghost_idx = interior_cell_idx + (1:num_ghost_cells);
        x_end = cfg.mesh.xc(end);
        x_ghost = x_end + (1:num_ghost_cells) * dx;
    else
        error('Invalid side specified. Use ''left'' or ''right''.');
    end

    h_end = w_padded(interior_cell_idx, 1);    % [m] Depth of the end cell
    hu_end = w_padded(interior_cell_idx, 2);   % [m^2/s] Discharge of the end cell
    z_end = cfg.bathyHandle(cfg, x_end);       % [m] Bottom of the end cell
    z_ghost = cfg.bathyHandle(cfg, x_ghost);   % [m] Bottom of the ghost cells

    for k = 1:num_ghost_cells
        if h_end > dry_tol
            % Free surface of the end cell over the bottom of the ghost cell
            h_ghost = max(0, h_end + (z_end - z_ghost(k)));
        else
            h_ghost = h_end;   % a dry end: dry ghost cells
        end
        if h_ghost > dry_tol
            w_padded(ghost_idx(k), :) = [h_ghost, hu_end];
        else
            w_padded(ghost_idx(k), :) = [h_ghost, 0];
        end
    end

    % Note: This BC does not account for incoming waves or information from outside
    % the domain. For more accurate non-reflecting boundaries, more sophisticated
    % methods (e.g., characteristic-based, sponge layers) may be used.

end