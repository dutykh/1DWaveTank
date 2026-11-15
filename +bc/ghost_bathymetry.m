% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +bc/ghost_bathymetry.m
%
% Purpose:
%   Gives the ghost cells at each end of the tank the bottom elevation that
%   matches the state their boundary condition puts in them, so that the
%   free surface eta = h + z of a ghost cell is that of the water it stands
%   for. The right-hand sides of +core call it after evaluating the bottom at
%   the centres of all the cells, ghost cells included:
%     - solid wall (bc.wall): the mirror image of the bottom, as the depth is
%       mirrored (bc.wall_bathymetry);
%     - wave-making end (bc.generating): the bottom of the cell next to the
%       end, copied into every ghost cell, as the condition puts one state,
%       a depth H0 + a sin(2 pi t/T), into all its ghost cells: the tank is
%       continued by a flat bottom beyond its end, and the free surface of the
%       ghost cells oscillates about the still-water level 0 when H0 is the
%       depth -z_e of the end cell at rest;
%     - periodic tank (bc.periodic): the periodic image of the bottom, as of
%       the state;
%     - open end (bc.open): the bottom function beyond the end, which bc.open
%       reads itself to give its ghost cells the free surface of the end cell.
%   The bottom function evaluated beyond a wall or a wave-making end would
%   instead give a ghost cell holding the depth h_1 of the first cell the free
%   surface
%       eta_0 = h_1 + z(x_0) = eta_1 + (z(x_0) - z(x_1)),
%   off by about |dz/dx| dx over a slope: the hydrostatic reconstruction on
%   the end would then see two different free surfaces, and a lake at rest
%   would be set in motion from that end. With the bottom given here, a lake
%   at rest (u = 0, eta constant) gives a ghost free surface equal to that of
%   the tank at every end (at a wave-making end, when cfg.param.H0 = -z_e);
%   the hydrostatic reconstruction then finds equal states on both sides of
%   the end, which is a steady state of every scheme of +core.
%
% Syntax:
%   z_padded = ghost_bathymetry(z_padded, cfg, num_ghost_cells)
%
% Inputs:
%   z_padded        - [1 x (N+2*num_ghost_cells)] bottom elevation at the cell
%                     centres, ghost cells included [m].
%   cfg             - [struct] configuration: cfg.mesh.N, cfg.bc.left.handle,
%                     cfg.bc.right.handle.
%   num_ghost_cells - [integer] number of ghost cells on each side.
%
% Outputs:
%   z_padded        - the same values, those of the ghost cells replaced as
%                     described above. An open end, or an end with another
%                     condition, keeps the bottom function.
%
% References:
%   - Audusse, E., Bouchut, F., Bristeau, M.-O., Klein, R., Perthame, B. (2004).
%     A fast and stable well-balanced scheme with hydrostatic reconstruction
%     for shallow water flows. SIAM J. Sci. Comput. 25(6), 2050-2065.
%   - LeVeque, R.J. (2002). Finite Volume Methods for Hyperbolic Problems.
%     Cambridge University Press, Section 7.3 (boundary conditions).
%
% Date:   15 November 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z_padded = ghost_bathymetry(z_padded, cfg, num_ghost_cells)

    N = cfg.mesh.N;          % [integer] Number of interior (physical) cells
    ng = num_ghost_cells;    % [integer] Ghost cells on each side

    left = cfg.bc.left.handle;
    right = cfg.bc.right.handle;

    % Periodic tank: the ghost cells beyond each end are images of the cells
    % at the other end, for the bottom as for the state (bc.periodic)
    if isequal(left, @bc.periodic) && isequal(right, @bc.periodic)
        z_padded(1:ng) = z_padded(N+1:N+ng);
        z_padded(N+ng+1:N+2*ng) = z_padded(ng+1:2*ng);
        return;
    end

    % Solid walls: mirror image of the bottom
    z_padded = bc.wall_bathymetry(z_padded, cfg, ng);

    % Wave-making LEFT end: the bottom of the first cell
    if isequal(left, @bc.generating)
        z_padded(1:ng) = z_padded(ng+1);
    end

    % Wave-making RIGHT end: the bottom of the last cell
    if isequal(right, @bc.generating)
        z_padded(N+ng+1:N+2*ng) = z_padded(N+ng);
    end

end
