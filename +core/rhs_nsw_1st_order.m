%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +core/rhs_nsw_1st_order.m
%
% Purpose:
%   Computes the right-hand side (RHS) for the 1D Nonlinear Shallow Water (NSW)
%   equations using a first-order finite volume (FV) scheme. This version
%   implements a well-balanced scheme using hydrostatic reconstruction as
%   described in Audusse et al. (based on the provided PDF document).
%   Handles boundary conditions, numerical fluxes, and source terms
%   (well-balanced bed slope, friction).
%
% Syntax:
%   dwdt_flat = rhs_nsw_1st_order(t, w_flat, cfg)
%
% Inputs:
%   t       - [scalar, double] Current simulation time [s].
%   w_flat  - [2N x 1, double] Flattened state vector [H1;...;HN; HU1;...;HUN].
%   cfg     - [struct] Configuration structure. Required fields:
%               cfg.mesh.N:   [integer] Number of spatial cells
%               cfg.mesh.dx:  [double] Cell width [m]
%               cfg.mesh.xc:  [1xN double] Cell center coordinates
%               cfg.phys.g:   [double] Acceleration due to gravity [m/s^2]
%               cfg.phys.dry_tolerance: [double] Tolerance for dry cells
%               cfg.numFlux:  [function handle] Numerical flux function
%               cfg.bc.left.handle:  [function handle] Left BC
%               cfg.bc.right.handle: [function handle] Right BC
%               cfg.bathyHandle: [function handle] Bathymetry
%             Optional:
%               cfg.phys.friction_model: [function handle] Friction model
%               cfg.phys.Cf: [double] Legacy friction coefficient
%
% Outputs:
%   dwdt_flat - [2N x 1, double] Flattened time derivative vector [dH/dt; dHU/dt].
%
% References:
%   - Based on well-balanced scheme concepts (e.g., Audusse et al.).
%   - LeVeque, R.J. (2002). Finite Volume Methods for Hyperbolic Problems.
%   - Toro, E.F. (2001). Shock-Capturing Methods for Free-Surface Shallow Flows.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   May 14, 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dwdt_flat = rhs_nsw_1st_order(t, w_flat, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Preparation and Parameter Extraction                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w_flat = w_flat(:); % Ensure w_flat is a column vector

    N = cfg.mesh.N;         % [integer] Number of spatial cells
    dx = cfg.mesh.dx;       % [m] Cell width
    g = cfg.phys.g;         % [m/s^2] Gravity
    dry_tol = cfg.phys.dry_tolerance; % Tolerance for dry cells

    % Reshape the flattened state vector w_flat into an N x 2 array [H, HU]
    w = [w_flat(1:N), w_flat(N+1:2*N)]; % H = w(:,1), HU = w(:,2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ghost Cells and Boundary Conditions                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_ghost_cells = 1; % Sufficient for 1st order fluxes and reconstruction
    Ng = num_ghost_cells;

    w_padded = zeros(N + 2*Ng, 2);
    w_padded(Ng+1 : N+Ng, :) = w;

    is_periodic = isequal(cfg.bc.left.handle, @bc.periodic) && isequal(cfg.bc.right.handle, @bc.periodic);

    if is_periodic
        w_padded = cfg.bc.left.handle(w_padded, t, 'both', cfg, Ng);
    else
        if isfield(cfg.bc.left, 'handle') && ~isempty(cfg.bc.left.handle)
            w_padded = cfg.bc.left.handle(w_padded, t, 'left', cfg, Ng);
        else
            warning('core:rhs:NoLeftBC', 'No left boundary condition handle specified.');
        end
        if isfield(cfg.bc.right, 'handle') && ~isempty(cfg.bc.right.handle)
            w_padded = cfg.bc.right.handle(w_padded, t, 'right', cfg, Ng);
        else
            warning('core:rhs:NoRightBC', 'No right boundary condition handle specified.');
        end
    end

    H_padded = w_padded(:,1);
    HU_padded = w_padded(:,2);
    U_padded = zeros(size(H_padded));
    wet_idx_padded = H_padded > dry_tol;
    U_padded(wet_idx_padded) = HU_padded(wet_idx_padded) ./ H_padded(wet_idx_padded);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bathymetry at Cell Centers (including ghost cells)         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_physical_centers = cfg.mesh.xc; % 1xN
    x_ghost_left_centers = zeros(1, Ng);
    for k_ghost = 1:Ng
        x_ghost_left_centers(k_ghost) = x_physical_centers(1) - (Ng - k_ghost + 1) * dx;
    end
    x_ghost_right_centers = zeros(1, Ng);
    for k_ghost = 1:Ng
        x_ghost_right_centers(k_ghost) = x_physical_centers(end) + k_ghost * dx;
    end
    x_all_cell_centers = [x_ghost_left_centers, x_physical_centers, x_ghost_right_centers]; % 1x(N+2*Ng)
    z_bc = cfg.bathyHandle(cfg, x_all_cell_centers); % Bathymetry at all cell centers
    z_bc = z_bc(:)'; % Ensure row vector 1x(N+2*Ng)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hydrostatic Reconstruction & Numerical Flux Calculation     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wL_recon = zeros(N+1, 2);
    wR_recon = zeros(N+1, 2);

    for k_int = 1:N+1 % Loop over N+1 interfaces
        % Cell i (left of interface k_int) and cell i+1 (right of interface k_int)
        idx_i = Ng + k_int -1; % Padded index for cell i (to the left of interface k_int)
        idx_ip1 = Ng + k_int;  % Padded index for cell i+1 (to the right of interface k_int)

        h_i = H_padded(idx_i);
        u_i = U_padded(idx_i);
        z_i = z_bc(idx_i);

        h_ip1 = H_padded(idx_ip1);
        u_ip1 = U_padded(idx_ip1);
        z_ip1 = z_bc(idx_ip1);
        
        max_z_interface = max(z_i, z_ip1);

        % Hydrostatic reconstruction for U_G (left state at interface) based on cell i
        h_G_k = max(0, h_i + z_i - max_z_interface);
        hu_G_k = h_G_k * u_i;
        wL_recon(k_int, :) = [h_G_k, hu_G_k];

        % Hydrostatic reconstruction for U_D (right state at interface) based on cell i+1
        h_D_k = max(0, h_ip1 + z_i - max_z_interface); % Error in PDF fixed: uses z_i instead of z_{i+1} from cell i for consistency with U_G
                                                      % Corrected based on standard hydrostatic reconstruction: h_D_k = max(0, h_ip1 + z_ip1 - max_z_interface);
        h_D_k_corrected = max(0, h_ip1 + z_ip1 - max_z_interface); % Corrected version
        hu_D_k = h_D_k_corrected * u_ip1;
        wR_recon(k_int, :) = [h_D_k_corrected, hu_D_k];
    end
    
    F_num = cfg.numFlux(wL_recon, wR_recon, cfg); % [(N+1) x 2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Flux Divergence (Spatial Derivative)                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dwdt_flux = -(F_num(2:N+1,:) - F_num(1:N,:)) / dx; % [N, 2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Source Terms (Well-Balanced Bed Slope, Friction)            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dwdt_source = zeros(N, 2);

    % --- Well-Balanced Bed Slope Source Term (Audusse et al. style) ---
    S_b_momentum = zeros(N,1);
    for i = 1:N % Loop over physical cells
        idx_cell_i = Ng + i; % Padded index for current cell i

        h_i_cell = H_padded(idx_cell_i);
        z_i_cell = z_bc(idx_cell_i);
        
        % Bathymetry of neighboring cells
        z_im1_cell = z_bc(idx_cell_i - 1); % z_{i-1}
        z_ip1_cell = z_bc(idx_cell_i + 1); % z_{i+1}

        % Reconstructed height at right interface of cell i (h_{i+1/2G})
        h_star_R_i = max(0, h_i_cell + z_i_cell - max(z_i_cell, z_ip1_cell));
        
        % Reconstructed height at left interface of cell i (h_{i-1/2D})
        % This is h at interface i-1/2, based on cell i's data looking towards cell i-1
        h_star_L_i = max(0, h_i_cell + z_i_cell - max(z_im1_cell, z_i_cell));
        
        S_b_momentum(i) = (g / (2 * dx)) * (h_star_R_i^2 - h_star_L_i^2);
    end
    dwdt_source(:, 2) = dwdt_source(:, 2) + S_b_momentum;

    % --- Friction Source Term (S_f) ---
    if isfield(cfg.phys, 'friction_model') && ~isempty(cfg.phys.friction_model)
        H_physical = w(:, 1);  % N x 1
        HU_physical = w(:, 2); % N x 1
        wet_indices_physical = H_physical > dry_tol;
        
        if any(wet_indices_physical)
            friction_term_values = cfg.phys.friction_model(H_physical(wet_indices_physical), HU_physical(wet_indices_physical), g, cfg);
            % The friction term from model is typically S_f for d(hu)/dt = ... + S_f (so sign is included)
            % Ensure friction_term_values is a column vector if not already
            dwdt_source(wet_indices_physical, 2) = dwdt_source(wet_indices_physical, 2) + friction_term_values(:);
        end
    elseif isfield(cfg.phys, 'Cf') && cfg.phys.Cf > 0 % Legacy friction
        H_physical = w(:, 1);
        HU_physical = w(:, 2);
        U_physical = zeros(N,1);
        wet_indices_physical = H_physical > dry_tol;
        U_physical(wet_indices_physical) = HU_physical(wet_indices_physical) ./ H_physical(wet_indices_physical);
        
        friction_legacy = zeros(N,1);
        friction_legacy(wet_indices_physical) = -g * cfg.phys.Cf^2 * abs(U_physical(wet_indices_physical)) .* U_physical(wet_indices_physical) ./ (H_physical(wet_indices_physical).^(1/3) + cfg.numerics.epsilon); % Manning like term (H^(1/3) in denominator)
        % This implies Sf_mom = -g * Cf^2 * U^2 * sign(U) * H / H^(7/3) = -g * Cf^2 * |U|U / H^(4/3)                                                                                                                                                                                   % So the denominator should be H^(4/3)
        friction_legacy(wet_indices_physical) = -g * cfg.phys.Cf^2 .* abs(U_physical(wet_indices_physical)) .* U_physical(wet_indices_physical) ./ (H_physical(wet_indices_physical).^(4/3) + cfg.numerics.epsilon);

        dwdt_source(wet_indices_physical, 2) = dwdt_source(wet_indices_physical, 2) + friction_legacy(wet_indices_physical);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Combine Flux and Source Terms, Flatten Output               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dwdt = dwdt_flux + dwdt_source; % [N, 2] Total time derivative
    dwdt_flat = [dwdt(:,1); dwdt(:,2)]; % [2N x 1] Flattened column vector

end