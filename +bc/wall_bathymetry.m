% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +bc/wall_bathymetry.m
%
% Purpose:
%   Gives the ghost cells of each end of the tank closed by a solid wall the
%   mirror image of the bottom elevation, the companion of bc.wall, which
%   mirrors the depth and reverses the velocity; the first-order and hybrid
%   right-hand sides of +core call it. A wall at x = x_w is a plane
%   of symmetry: data extended by
%       h(x_w + s) = h(x_w - s),  u(x_w + s) = -u(x_w - s),  z(x_w + s) = z(x_w - s)
%   stay symmetric under the shallow water equations, so that the velocity of
%   the extended solution vanishes at x_w (the argument of LeVeque 2002,
%   Section 7.3.3, with the bottom extended like the depth). With the ghost
%   cell 0 equal to the mirror image of the cell 1 next to a left wall, the
%   hydrostatic reconstruction of Audusse et al. (2004) at the interface on the
%   wall gives
%       z* = max(z_0, z_1) = z_1,   h*_- = h*_+ = h_1,   u_- = -u_1,  u_+ = u_1,
%   a Riemann problem symmetric about the wall, whose solution has zero
%   velocity on the wall. For a lake at rest (h + z constant, u = 0) both
%   states equal (h_1, 0), the flux is (0, g h_1^2/2) by consistency, and the
%   lake at rest is an exact steady state of the first-order scheme next to a
%   wall over any bottom. A ghost bottom taken from the bottom function beyond
%   the wall would leave depths differing by about |dz/dx| dx on the two sides
%   of the wall, hence a flow through the wall from rest.
%
% Syntax:
%   z_padded = wall_bathymetry(z_padded, cfg, num_ghost_cells)
%
% Inputs:
%   z_padded        - [1 x (N+2*num_ghost_cells)] bottom elevation at the cell
%                     centres, ghost cells included [m].
%   cfg             - [struct] configuration: cfg.mesh.N, cfg.bc.left.handle,
%                     cfg.bc.right.handle.
%   num_ghost_cells - [integer] number of ghost cells on each side.
%
% Outputs:
%   z_padded        - the same values, those of the ghost cells of each wall
%                     replaced by the mirror image of the cells next to it.
%
% References:
%   - Audusse, E., Bouchut, F., Bristeau, M.-O., Klein, R., Perthame, B. (2004).
%     A fast and stable well-balanced scheme with hydrostatic reconstruction
%     for shallow water flows. SIAM J. Sci. Comput. 25(6), 2050-2065.
%   - LeVeque, R.J. (2002). Finite Volume Methods for Hyperbolic Problems.
%     Cambridge University Press, Section 7.3.3 (solid walls).
%
% Date:   15 November 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z_padded = wall_bathymetry(z_padded, cfg, num_ghost_cells)

    N = cfg.mesh.N;          % [integer] Number of interior (physical) cells
    ng = num_ghost_cells;    % [integer] Ghost cells on each side

    % LEFT wall: the ghost cell ng+1-k mirrors the interior cell ng+k, k = 1..ng
    if isequal(cfg.bc.left.handle, @bc.wall)
        z_padded(1:ng) = z_padded(2*ng:-1:ng+1);
    end

    % RIGHT wall: the ghost cell N+ng+k mirrors the interior cell N+ng+1-k, k = 1..ng
    if isequal(cfg.bc.right.handle, @bc.wall)
        z_padded(N+ng+1:N+2*ng) = z_padded(N+ng:-1:N+1);
    end

end
