%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +bc/open.m
%
% Purpose:
%   Implements an "open" boundary condition for the 1DWaveTank code, through
%   which waves leave the tank: the ghost cells receive the state that the
%   characteristics of the shallow water equations carry there. Over a flat
%   bottom the Riemann invariants
%       R+ = u + 2 c,   R- = u - 2 c,   c = sqrt(g h),
%   are constant along the characteristics dx/dt = u + c and u - c. At a
%   right end (subcritical flow), R+ leaves the tank and is taken from the
%   end cell, while R- enters it and is taken from the water at rest outside
%   the tank, whose free surface is the still-water level eta0:
%       R+ = u_e + 2 sqrt(g h_x),   R- = - 2 sqrt(g h_0),
%       u_g = (R+ + R-)/2,          c_g = (R+ - R-)/4,
%   with, over the bottom z_g of the ghost cell,
%       h_x = max(0, eta_e - z_g)   (the free surface of the end cell over z_g),
%       h_0 = max(0, eta0 - z_g)    (the water at rest outside),
%   so that the invariants are compared over the same bottom; at a left end
%   R- leaves and R+ enters, R+ = 2 sqrt(g h_0). The depth of the ghost cell
%   is h_g = h_0 (c_g/sqrt(g h_0))^2 = c_g^2/g, and its discharge h_g u_g. A
%   lake at rest (eta_e = eta0, u_e = 0) gives u_g = 0 and h_g = h_0
%   exactly: the free surface of the ghost cells is eta0 and the lake stays
%   at rest, over any bottom. A wave reaching the end leaves the tank with
%   its outgoing invariant, while the incoming one stays that of still
%   water, so that the end reflects much less than the copy of the end cell
%   (zeroth-order extrapolation). Where the flow leaves the tank faster than
%   the waves (u_e ⩾ c at a right end, u_e ⩽ -c at a left end), both
%   invariants come from the tank and the ghost cells receive the free
%   surface and the velocity of the end cell; a dry end cell gives dry ghost
%   cells, and so does c_g ⩽ 0. The ghost cells keep the bottom function
%   beyond the end of the tank (bc.ghost_bathymetry).
%
%   The still-water level eta0 is cfg.bc.left.param.eta0 or
%   cfg.bc.right.param.eta0 [m]; cfg.validate_config sets it, when it is not
%   given, to the free surface of the end cell in the initial state. Without
%   it, the ghost cells receive the free surface and the discharge of the
%   end cell, h_g = max(0, h_e + z_e - z_g), (hu)_g = (hu)_e.
%
% Syntax:
%   w_padded = open(w_padded, t, side, cfg, num_ghost_cells)
%
% Inputs:
%   w_padded        - [N+2*num_ghost_cells x 2] array. State vector including ghost cells.
%   t               - [scalar] Current simulation time (seconds). (Unused)
%   side            - [char] 'left' or 'right'. Which boundary to apply.
%   cfg             - [struct] Configuration structure: cfg.mesh.N, cfg.mesh.dx,
%                     cfg.mesh.xc, cfg.phys.g, cfg.phys.dry_tolerance,
%                     cfg.bathyHandle, and cfg.bc.(side).param.eta0.
%   num_ghost_cells - [integer] Number of ghost cells to fill (>=1).
%
% Outputs:
%   w_padded        - [N+2*num_ghost_cells x 2] array. State vector with ghost cells filled.
%
% Dependencies:
%   Requires cfg.mesh and the bottom function cfg.bathyHandle.
%
% References:
%   - Stoker, J.J. (1957). Water Waves: The Mathematical Theory with Applications.
%   - LeVeque, R.J. (2002). Finite Volume Methods for Hyperbolic Problems.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w_padded = open(w_padded, t, side, cfg, num_ghost_cells)

    N = cfg.mesh.N;                    % [integer] Number of interior (physical) cells
    dx = cfg.mesh.dx;                  % [m] Cell width
    g = cfg.phys.g;                    % [m/s^2] Gravitational acceleration
    dry_tol = cfg.phys.dry_tolerance;  % [m] Depth below which a cell is dry

    if strcmp(side, 'left')
        %---------------------------------------------
        % LEFT boundary: fill ghost cells 1:num_ghost_cells, outermost first
        interior_cell_idx = num_ghost_cells + 1;          % First interior cell
        ghost_idx = 1:num_ghost_cells;
        x_end = cfg.mesh.xc(1);                            % [m] Centre of the end cell
        x_ghost = x_end - (num_ghost_cells:-1:1) * dx;     % [m] Centres of the ghost cells
        outward = -1;                                      % the tank lies to the right
    elseif strcmp(side, 'right')
        %---------------------------------------------
        % RIGHT boundary: fill ghost cells N+num_ghost_cells+1 : N+2*num_ghost_cells
        interior_cell_idx = N + num_ghost_cells;           % Last interior cell
        ghost_idx = interior_cell_idx + (1:num_ghost_cells);
        x_end = cfg.mesh.xc(end);
        x_ghost = x_end + (1:num_ghost_cells) * dx;
        outward = 1;                                       % the tank lies to the left
    else
        error('Invalid side specified. Use ''left'' or ''right''.');
    end

    % Still-water level outside the tank, when the configuration gives one
    params = cfg.bc.(side).param;
    have_eta0 = isfield(params, 'eta0') && ~isempty(params.eta0);

    h_end = w_padded(interior_cell_idx, 1);    % [m] Depth of the end cell
    hu_end = w_padded(interior_cell_idx, 2);   % [m^2/s] Discharge of the end cell
    z_end = cfg.bathyHandle(cfg, x_end);       % [m] Bottom of the end cell
    z_ghost = cfg.bathyHandle(cfg, x_ghost);   % [m] Bottom of the ghost cells

    eta_end = h_end + z_end;                   % [m] Free surface of the end cell

    for k = 1:num_ghost_cells
        if h_end <= dry_tol
            % A dry end: dry ghost cells
            w_padded(ghost_idx(k), :) = [h_end, 0];
            continue;
        end

        if ~have_eta0
            % Zeroth-order extrapolation of the free surface and the discharge
            h_ext = max(0, h_end + (z_end - z_ghost(k)));
            if h_ext > dry_tol
                w_padded(ghost_idx(k), :) = [h_ext, hu_end];
            else
                w_padded(ghost_idx(k), :) = [h_ext, 0];
            end
            continue;
        end

        % Free surface of the end cell over the bottom of the ghost cell
        h_ext = max(0, eta_end - z_ghost(k));
        if h_ext <= dry_tol
            w_padded(ghost_idx(k), :) = [h_ext, 0];
            continue;
        end

        u_end = hu_end / h_end;                       % [m/s] Velocity of the end cell
        c_ext = sqrt(g * h_ext);                      % [m/s] Its celerity over z_g
        if outward * u_end >= c_ext
            % Supercritical outflow: both invariants leave the tank
            w_padded(ghost_idx(k), :) = [h_ext, h_ext * u_end];
            continue;
        end

        h_rest = max(0, params.eta0 - z_ghost(k));     % [m] Still water over z_g
        c_rest = sqrt(g * h_rest);
        % Outgoing invariant from the tank, incoming invariant of still water
        R_out = u_end + outward * 2 * c_ext;
        R_in = -outward * 2 * c_rest;
        u_ghost = (R_out + R_in) / 2;
        c_ghost = outward * (R_out - R_in) / 4;
        if c_ghost <= 0
            w_padded(ghost_idx(k), :) = [0, 0];
            continue;
        end
        if c_rest > 0
            h_ghost = h_rest * (c_ghost / c_rest)^2;  % h_rest exactly when c_ghost = c_rest
        else
            h_ghost = c_ghost^2 / g;
        end
        if h_ghost > dry_tol
            w_padded(ghost_idx(k), :) = [h_ghost, h_ghost * u_ghost];
        else
            w_padded(ghost_idx(k), :) = [h_ghost, 0];
        end
    end

end
