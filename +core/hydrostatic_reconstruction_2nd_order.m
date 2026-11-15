% Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +core/hydrostatic_reconstruction_2nd_order.m
%
% Purpose:
%   Second-order hydrostatic reconstruction of Audusse, Bouchut, Bristeau,
%   Klein and Perthame (2004, Section 3), used by the hybrid and high-order
%   right-hand sides: from the cell averages it computes the two states
%   passed to the numerical flux at every interface, and the bed slope of
%   every cell, so that the semi-discrete scheme
%       dx dU_i/dt + F_{i+1/2} - F_{i-1/2} = (0, S_i)
%   keeps a lake at rest exactly, keeps the depth nonnegative at a shoreline,
%   and is second-order accurate on smooth solutions.
%
%   Reconstruction. The depth h, the free surface eta = h + z and the
%   discharge hu are reconstructed at both edges of every cell, component by
%   component, by the reconstruction of the configuration; the bottom at
%   each edge is then
%       z_{i,l} = eta_{i,l} - h_{i,l},     z_{i,r} = eta_{i,r} - h_{i,r},
%   and the velocity u = hu/h (zero where h is below the dry tolerance).
%   Of the three quantities h, z and h + z, two are reconstructed and the
%   third follows; the pair (h, h + z) is the one choice that keeps both the
%   lake at rest and the depth nonnegative at an interface between a wet and
%   a dry cell (Audusse et al. 2004, Section 3 and Figure 1): a limited
%   reconstruction of h stays nonnegative, and that of eta is flat in a wet
%   cell at rest next to a dry one.
%
%   Where the configuration asks for a reconstruction in characteristic
%   variables (the Riemann invariants u -/+ 2 sqrt(g h) of MUSCL, or the
%   characteristic fields of WENO5, PPM, MP5 and THINC), the velocity u is
%   reconstructed instead of the discharge, component by component as well.
%   Characteristic variables cannot be used here: the invariants vary with
%   the depth, hence with the bottom, in a lake at rest, so that their
%   reconstruction is not flat at rest, and a depth taken from them would
%   not match the reconstructed free surface. A limited velocity stays
%   between the velocities of the neighbouring cells, which keeps a wetting
%   front from running ahead of the water, where the velocity hu/h of a
%   reconstructed discharge may grow without bound as the reconstructed
%   depth vanishes.
%
%   A cell whose depth is below the dry tolerance is dry: its depth at both
%   edges is zero, so that the thin layer of water it may hold stays in place
%   until water arrives, and its bed slope vanishes.
%
%   Hydrostatic reconstruction (equations (3.2) to (3.5)). At the interface
%   i+1/2, between the right edge of cell i and the left edge of cell i+1,
%       z_{i+1/2}  = max(z_{i,r}, z_{i+1,l}),
%       h_{i+1/2-} = max(0, eta_{i,r}   - z_{i+1/2}),
%       h_{i+1/2+} = max(0, eta_{i+1,l} - z_{i+1/2}),
%       U_{i+1/2-} = (h_{i+1/2-}, h_{i+1/2-} u_{i,r}),
%       U_{i+1/2+} = (h_{i+1/2+}, h_{i+1/2+} u_{i+1,l}),
%   and F_{i+1/2} is the numerical flux of the two states U_{i+1/2-} and
%   U_{i+1/2+}. Since 0 ⩽ h_{i+1/2-} ⩽ h_{i,r} and 0 ⩽ h_{i+1/2+} ⩽
%   h_{i+1,l}, a forward Euler step (and a strong-stability-preserving
%   Runge-Kutta step, a convex combination of such steps) keeps the depth
%   nonnegative under half the CFL condition of the first-order scheme,
%   whenever the reconstructed depths are nonnegative and have the mean h_i
%   in cell i, as those of MUSCL do (Theorem 3.1).
%
%   Bed slope (equations (3.6) to (3.8)). In the momentum equation of cell i,
%       S_i = g/2 (h_{i+1/2-}^2 - h_{i,r}^2) + g/2 (h_{i,l}^2 - h_{i-1/2+}^2)
%             + g (h_{i,l} + h_{i,r})/2 (z_{i,l} - z_{i,r}).
%   The first two terms vanish where the bottom is continuous across the
%   interfaces; the last one is the centred bed slope - g h z_x dx.
%
%   Lake at rest. If u = 0 and eta = eta_0 in every wet cell, the
%   reconstruction gives u = 0 at every edge and eta = eta_0 at every edge
%   of a wet cell (a limited reconstruction of a constant is that constant,
%   and in a wet cell next to a dry one the limiter takes the zero slope of
%   eta on the wet side). The two states at every interface are then equal,
%   h_{i+1/2-} = h_{i+1/2+} = max(0, eta_0 - z_{i+1/2}), and consistency of
%   the numerical flux gives F_{i+1/2} = (0, g h_{i+1/2-}^2/2). Since
%   z_{i,l} - z_{i,r} = h_{i,r} - h_{i,l},
%       S_i = g/2 (h_{i+1/2-}^2 - h_{i-1/2+}^2) = F^{hu}_{i+1/2} - F^{hu}_{i-1/2},
%   and every cell is at rest. In a dry cell next to the shoreline, the
%   depth reconstructed at the edges is zero, and the free surface
%   reconstructed at its wet edge lies at or above eta_0, the limited slope
%   being at most twice the smaller difference of eta (a bound that every
%   limiter of +reconstruct/+limiters satisfies): the bottom there is at or
%   above eta_0, the depths on both sides of the shoreline interface are
%   zero, and no water moves.
%
% Syntax:
%   [wL_hydro, wR_hydro, S_bed] = hydrostatic_reconstruction_2nd_order(w_padded, z_padded, cfg)
%
% Inputs:
%   w_padded - [(N+2*ng) x 2, double] depth and discharge [H, HU] of all the
%              cells, ghost cells filled by the boundary conditions.
%   z_padded - [1 x (N+2*ng)] or [(N+2*ng) x 1] bottom elevation at the cell
%              centres, ghost cells given by bc.ghost_bathymetry [m].
%   cfg      - [struct] configuration: cfg.mesh.N, cfg.mesh.dx, cfg.phys.g,
%              cfg.phys.dry_tolerance, cfg.bc.num_ghost_cells, and
%              cfg.reconstruct.handle (a function of +reconstruct; without
%              it the reconstruction is constant in each cell and the scheme
%              is the first-order hydrostatic reconstruction scheme).
%
% Outputs:
%   wL_hydro - [(N+1) x 2] the states U_{i+1/2-} at the N+1 interfaces,
%              from the left boundary of the tank to the right one.
%   wR_hydro - [(N+1) x 2] the states U_{i+1/2+} at the same interfaces.
%   S_bed    - [N x 1] the bed slope S_i/dx of each cell [m^2/s^2], to be
%              added to d(HU)/dt.
%
% References:
%   - Audusse, E., Bouchut, F., Bristeau, M.-O., Klein, R., Perthame, B. (2004).
%     A fast and stable well-balanced scheme with hydrostatic reconstruction
%     for shallow water flows. SIAM J. Sci. Comput. 25(6), 2050-2065,
%     Section 3.
%
% Date:   15 November 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [wL_hydro, wR_hydro, S_bed] = hydrostatic_reconstruction_2nd_order(w_padded, z_padded, cfg)

    N = cfg.mesh.N;                    % [integer] Number of cells
    ng = cfg.bc.num_ghost_cells;       % [integer] Ghost cells on each side
    dx = cfg.mesh.dx;                  % [m] Cell width
    g = cfg.phys.g;                    % [m/s^2] Gravity
    dry_tol = cfg.phys.dry_tolerance;  % [m] Depth below which a state is dry

    if isfield(cfg, 'reconstruct') && isfield(cfg.reconstruct, 'handle') && ~isempty(cfg.reconstruct.handle)
        reconstruct_handle = cfg.reconstruct.handle;
    else
        reconstruct_handle = @reconstruct.none;
    end
    [cfg_comp, componentwise] = componentwise_config(cfg, reconstruct_handle);

    z_padded = z_padded(:);
    eta_padded = w_padded(:, 1) + z_padded;   % free surface of every cell

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Depth, discharge and free surface at the edges of the cells %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Row k of the outputs of a reconstruction refers to the interface k,
    % between the cells k-1 and k: the first array holds the right edge of
    % cell k-1 (the side "-"), the second the left edge of cell k (side "+").
    [w_minus, w_plus] = reconstruct_handle(w_padded, cfg_comp);

    % The free surface is reconstructed by the same method, applied component
    % by component to the pair (1, eta): the first column, a depth of one
    % metre in every cell, is never dry, so that the treatment of dry states
    % inside the reconstruction leaves the second column, eta, untouched.
    [eta_minus, eta_plus] = reconstruct_handle([ones(size(eta_padded)), eta_padded], cfg_comp);
    eta_minus = eta_minus(:, 2);
    eta_plus = eta_plus(:, 2);

    h_minus = w_minus(:, 1);   % h_{i,r}
    h_plus = w_plus(:, 1);     % h_{i+1,l}

    % A cell whose depth is below the dry tolerance is dry: its depth at both
    % edges is zero, so that the thin layer of water left on dry ground stays
    % where it is until water arrives, and the bed slope of the cell vanishes
    dry = w_padded(:, 1) <= dry_tol;          % cells and ghost cells
    h_minus(dry(ng:ng+N)) = 0;               % cell to the left of each interface
    h_plus(dry(ng+1:ng+N+1)) = 0;            % cell to the right of each interface

    % Velocity at the edges: from the discharge reconstructed with the depth,
    % or, where the configuration asks for a reconstruction in characteristic
    % variables, the velocity itself, reconstructed by the same method
    u_minus = zeros(N+1, 1);
    u_plus = zeros(N+1, 1);
    if componentwise
        wet = w_minus(:, 1) > dry_tol;
        u_minus(wet) = w_minus(wet, 2) ./ w_minus(wet, 1);   % u_{i,r}
        wet = w_plus(:, 1) > dry_tol;
        u_plus(wet) = w_plus(wet, 2) ./ w_plus(wet, 1);      % u_{i+1,l}
    else
        u_padded = zeros(size(eta_padded));
        wet = w_padded(:, 1) > dry_tol;
        u_padded(wet) = w_padded(wet, 2) ./ w_padded(wet, 1);
        [u_rec_minus, u_rec_plus] = reconstruct_handle([ones(size(u_padded)), u_padded], cfg_comp);
        wet = h_minus > dry_tol;
        u_minus(wet) = u_rec_minus(wet, 2);                  % u_{i,r}
        wet = h_plus > dry_tol;
        u_plus(wet) = u_rec_plus(wet, 2);                    % u_{i+1,l}
    end

    % Bottom at the edges of the cells: z = eta - h
    z_minus = eta_minus - h_minus;   % z_{i,r}
    z_plus = eta_plus - h_plus;      % z_{i+1,l}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hydrostatic reconstruction at the interfaces, (3.2)-(3.5)    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z_interface = max(z_minus, z_plus);                % z_{i+1/2}
    hs_minus = max(0, eta_minus - z_interface);        % h_{i+1/2-}
    hs_plus = max(0, eta_plus - z_interface);          % h_{i+1/2+}

    wL_hydro = [hs_minus, hs_minus .* u_minus];        % U_{i+1/2-}
    wR_hydro = [hs_plus, hs_plus .* u_plus];           % U_{i+1/2+}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bed slope of each cell, (3.6)-(3.8)                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cell i has its left edge on the interface i (side "+") and its right
    % edge on the interface i+1 (side "-").
    h_l = h_plus(1:N);         % h_{i,l}
    h_r = h_minus(2:N+1);      % h_{i,r}
    z_l = z_plus(1:N);         % z_{i,l}
    z_r = z_minus(2:N+1);      % z_{i,r}

    S_right = 0.5 * g * (hs_minus(2:N+1).^2 - h_r.^2);   % S_{i+1/2-}
    S_left = 0.5 * g * (h_l.^2 - hs_plus(1:N).^2);       % S_{i-1/2+}
    S_centre = g * 0.5 * (h_l + h_r) .* (z_l - z_r);     % S_{c,i}

    S_bed = (S_right + S_left + S_centre) / dx;

end

function [cfg_comp, componentwise] = componentwise_config(cfg, reconstruct_handle)
    % The configuration with every reconstruction of +reconstruct set to work
    % component by component, and whether the given one does so already
    % (MUSCL, WENO5 and THINC read cfg.reconstruct.characteristic; PPM and MP5
    % read ppm_mode and mp5_mode, whose default is 'characteristic')
    if isfield(cfg, 'reconstruct') && isstruct(cfg.reconstruct)
        rc = cfg.reconstruct;
    else
        rc = struct();
    end
    componentwise = ~(isfield(rc, 'characteristic') && ~isempty(rc.characteristic) && rc.characteristic);
    if isequal(reconstruct_handle, @reconstruct.ppm)
        componentwise = isfield(rc, 'ppm_mode') && strcmpi(rc.ppm_mode, 'component');
    elseif isequal(reconstruct_handle, @reconstruct.mp5)
        componentwise = isfield(rc, 'mp5_mode') && strcmpi(rc.mp5_mode, 'component');
    elseif isequal(reconstruct_handle, @reconstruct.none) || isequal(reconstruct_handle, @reconstruct.eno2) ...
            || isequal(reconstruct_handle, @reconstruct.uno2)
        componentwise = true;
    end
    rc.characteristic = false;
    rc.ppm_mode = 'component';
    rc.mp5_mode = 'component';
    cfg_comp = cfg;
    cfg_comp.reconstruct = rc;
end
