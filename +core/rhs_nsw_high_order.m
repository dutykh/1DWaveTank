%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +core/rhs_nsw_high_order.m
%
% Purpose:
%   Computes the right-hand side (RHS) for the 1D Nonlinear Shallow Water (NSW)
%   equations by the second-order well-balanced finite volume scheme of
%   Audusse, Bouchut, Bristeau, Klein and Perthame (2004, Section 3). The
%   depth h, the free surface h + z and the discharge hu are reconstructed at
%   the edges of the cells by the reconstruction of the configuration, the
%   bottom at each edge being the difference of the first two; the
%   hydrostatic reconstruction gives the two states passed to the numerical
%   flux at each interface, and the bed slope is the sum of the hydrostatic
%   terms at the interfaces and of the centred term g (h_l + h_r)/2 (z_l - z_r)
%   of each cell (core.hydrostatic_reconstruction_2nd_order). The scheme is
%   well balanced: a lake at rest (u = 0, h + z constant where there is water)
%   is an exact steady state over any bottom, at a shoreline, next to a wall
%   and at an open or wave-making end, the ghost cells holding the free
%   surface of the water they stand for (bc.ghost_bathymetry: the mirror
%   image of the bottom at a wall; bc.open: the free surface of the end cell
%   over the bottom function).
%   With a limited reconstruction (MUSCL) the depth stays nonnegative, and
%   the scheme is second-order accurate on smooth solutions.
%
% Syntax:
%   dwdt_flat = rhs_nsw_high_order(t, w_flat, cfg)
%
% Inputs:
%   t       - [scalar, double] Current simulation time [s].
%   w_flat  - [2N x 1, double] Flattened state vector [H1;...;HN; HU1;...;HUN].
%   cfg     - [struct] Configuration structure.
%
% Outputs:
%   dwdt_flat - [2N x 1, double] Flattened time derivative vector [dH/dt; dHU/dt].
%
% References:
%   - Audusse, E., Bouchut, F., Bristeau, M.-O., Klein, R., Perthame, B. (2004).
%     A fast and stable well-balanced scheme with hydrostatic reconstruction
%     for shallow water flows. SIAM J. Sci. Comput. 25(6), 2050-2065.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   April 24, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dwdt_flat = rhs_nsw_high_order(t, w_flat, cfg)
    % Extract parameters
    ng = cfg.bc.num_ghost_cells;
    N = cfg.mesh.N;
    g = cfg.phys.g;
    dry_tol = cfg.phys.dry_tolerance;
    dx = cfg.mesh.dx;

    % Reshape state vector from flat to structured
    w_flat = w_flat(:);
    w = [w_flat(1:N), w_flat(N+1:2*N)];

    % Apply boundary conditions with ghost cells
    w_padded = zeros(N + 2*ng, 2);
    w_padded(ng+1 : N+ng, :) = w;
    w_padded = apply_boundary_conditions(w_padded, t, cfg);

    % Get bathymetry at all cell centers
    x_cell_centers = get_all_cell_centers(cfg, ng);
    z_b = cfg.bathyHandle(cfg, x_cell_centers);
    z_b = z_b(:)';  % Ensure row vector
    % The ghost cells take the bottom that matches the state of their boundary
    % condition (mirror image at a wall, copy of the end cell at a wave-making
    % end, periodic image, bottom function at an open end, whose ghost depth
    % bc.open sets from the free surface), so that their free surface is that
    % of the water they stand for
    z_b = bc.ghost_bathymetry(z_b, cfg, ng);

    % Second-order hydrostatic reconstruction: interface states and bed slope
    [wL_hydro, wR_hydro, S_bed] = core.hydrostatic_reconstruction_2nd_order(w_padded, z_b, cfg);

    % Calculate numerical fluxes
    F_num = cfg.numFlux(wL_hydro, wR_hydro, cfg);

    % Compute flux divergence term: -(F_{i+1/2} - F_{i-1/2})/dx
    dwdt_flux = -(F_num(2:N+1,:) - F_num(1:N,:)) / dx;

    % Bed slope, in the momentum equation only
    dwdt_source = zeros(N, 2);
    dwdt_source(:, 2) = S_bed;

    % Add friction source term if specified
    if isfield(cfg.phys, 'friction_model') && ~isempty(cfg.phys.friction_model)
        H = w(:, 1);
        HU = w(:, 2);
        wet_indices = H > dry_tol;

        if any(wet_indices)
            friction_term = cfg.phys.friction_model(H(wet_indices), HU(wet_indices), g, cfg);
            dwdt_source(wet_indices, 2) = dwdt_source(wet_indices, 2) + friction_term(:);
        end
    end

    % Combine flux and source terms
    dwdt = dwdt_flux + dwdt_source;

    % Flatten output for ODE solver
    dwdt_flat = [dwdt(:,1); dwdt(:,2)];
end

function x_centers = get_all_cell_centers(cfg, ng)
    % Get coordinates of cell centers including ghost cells
    dx = cfg.mesh.dx;
    x_interior = cfg.mesh.xc;

    % Ghost cells to the left
    x_left_ghost = zeros(1, ng);
    for i = 1:ng
        x_left_ghost(i) = x_interior(1) - (ng-i+1) * dx;
    end

    % Ghost cells to the right
    x_right_ghost = zeros(1, ng);
    for i = 1:ng
        x_right_ghost(i) = x_interior(end) + i * dx;
    end

    % Combine all cell centers
    x_centers = [x_left_ghost, x_interior, x_right_ghost];
end

function w_padded = apply_boundary_conditions(w_padded, t, cfg)
    % Apply boundary conditions to fill ghost cells
    is_periodic = isequal(cfg.bc.left.handle, @bc.periodic) && isequal(cfg.bc.right.handle, @bc.periodic);
    ng = cfg.bc.num_ghost_cells;

    if is_periodic
        w_padded = cfg.bc.left.handle(w_padded, t, 'both', cfg, ng);
    else
        if isfield(cfg.bc.left, 'handle') && ~isempty(cfg.bc.left.handle)
            w_padded = cfg.bc.left.handle(w_padded, t, 'left', cfg, ng);
        end
        if isfield(cfg.bc.right, 'handle') && ~isempty(cfg.bc.right.handle)
            w_padded = cfg.bc.right.handle(w_padded, t, 'right', cfg, ng);
        end
    end
end
