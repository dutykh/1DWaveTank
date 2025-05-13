function dwdt_flat = rhs_nsw_high_order(t, w_flat, cfg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +core/rhs_nsw_high_order.m
%
% Purpose:
%   Computes the right-hand side (RHS) for the 1D Nonlinear Shallow Water (NSW)
%   equations using a high-order finite volume scheme with reconstruction.
%   Handles boundary conditions, performs reconstruction, calculates numerical
%   fluxes, and computes source terms (friction, bed slope).
%
% Syntax:
%   dwdt_flat = rhs_nsw_high_order(t, w_flat, cfg)
%
% Inputs:
%   t       - [scalar, double] Current simulation time [s].
%   w_flat  - [2N x 1, double] Flattened state vector [H1;...;HN; HU1;...;HUN].
%   cfg     - [struct] Configuration structure. Required fields:
%               cfg.mesh.N:   [integer] Number of spatial cells
%               cfg.mesh.dx:  [double] Cell width [m]
%               cfg.phys.g:   [double] Acceleration due to gravity [m/s^2]
%               cfg.phys.dry_tolerance:  [double] Threshold for dry cells
%               cfg.numFlux:  [function handle] Numerical flux function
%               cfg.bc.left.handle:  [function handle] Left BC
%               cfg.bc.right.handle: [function handle] Right BC
%               cfg.reconstruct.handle: [function handle] Reconstruction method
%             Optional:
%               cfg.bathyHandle: [function handle] Bathymetry (default: flat)
%
% Outputs:
%   dwdt_flat - [2N x 1, double] Flattened time derivative vector [dH/dt; dHU/dt].
%
% Dependencies:
%   Expects correct configuration and function handles for BCs, fluxes, and reconstruction.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   April 24, 2025
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Preparation and Parameter Extraction                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w_flat = w_flat(:); % Ensure w_flat is a column vector

    N = cfg.mesh.N;         % [integer] Number of spatial cells
    dx = cfg.mesh.dx;       % [m] Cell width
    g = cfg.phys.g;         % [m/s^2] Gravity

    % Determine number of ghost cells required by reconstruction method
    if ~isfield(cfg.bc, 'num_ghost_cells')
        if isfield(cfg.reconstruct, 'order')
            ng = max(1, cfg.reconstruct.order);
            cfg.bc.num_ghost_cells = ng;
        else
            ng = 1; % Default: 1st order requires 1 ghost cell
            cfg.bc.num_ghost_cells = ng;
        end
    else
        ng = cfg.bc.num_ghost_cells;
    end

    % Reshape the flattened state vector w_flat into an N x 2 array [H, HU]
    % H = w(:,1), HU = w(:,2)
    w = [w_flat(1:N), w_flat(N+1:2*N)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ghost Cells and Boundary Conditions                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pad the domain with ghost cells (zeros for now)
    w_padded = zeros(N + 2*ng, 2); % [N+2*ng, 2] array
    w_padded(ng+1 : N+ng, :) = w; % Fill interior domain data

    % --- Apply Boundary Conditions ---
    % Check for periodic BCs first
    is_periodic = isequal(cfg.bc.left.handle, @bc.periodic) && isequal(cfg.bc.right.handle, @bc.periodic);

    if is_periodic
        % Apply periodic BCs once for both sides
        w_padded = cfg.bc.left.handle(w_padded, t, 'both', cfg, ng);
    else
        % Apply left BC
        if isfield(cfg.bc.left, 'handle') && ~isempty(cfg.bc.left.handle)
            w_padded = cfg.bc.left.handle(w_padded, t, 'left', cfg, ng);
        else
            warning('core:rhs:NoLeftBC', 'No left boundary condition handle specified.');
        end

        % Apply right BC
        if isfield(cfg.bc.right, 'handle') && ~isempty(cfg.bc.right.handle)
            w_padded = cfg.bc.right.handle(w_padded, t, 'right', cfg, ng);
        else
            warning('core:rhs:NoRightBC', 'No right boundary condition handle specified.');
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare Bathymetry for All Cells (including ghost cells)   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_physical_cells = cfg.mesh.xc; % Physical cell centers
    dx_val = cfg.mesh.dx;
    x_ghost_left = zeros(1, ng);
    for k_ghost = 1:ng
        x_ghost_left(k_ghost) = x_physical_cells(1) - (ng - k_ghost + 1) * dx_val;
    end
    x_ghost_right = zeros(1, ng);
    for k_ghost = 1:ng
        x_ghost_right(k_ghost) = x_physical_cells(end) + k_ghost * dx_val;
    end
    x_all_cell_centers = [x_ghost_left, x_physical_cells, x_ghost_right];
    z_cell_padded = cfg.bathyHandle(cfg, x_all_cell_centers); % Bathymetry at all cell centers
    z_cell_padded = z_cell_padded(:)'; % Ensure it's a row vector for consistent indexing

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reconstruction at Cell Interfaces                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine which reconstruction method to use
    if ~isfield(cfg, 'reconstruct') || ~isfield(cfg.reconstruct, 'handle') || isempty(cfg.reconstruct.handle)
        % Default to no reconstruction (1st order)
        reconstruct_handle = @reconstruct.none;
    else
        reconstruct_handle = cfg.reconstruct.handle;
    end
    
    % Apply the reconstruction method to get left and right states at interfaces
    [wL_interface, wR_interface] = reconstruct_handle(w_padded, cfg);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hydrostatic Reconstruction at Cell Interfaces (Well-balanced)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wL_interface_recon = wL_interface;
    wR_interface_recon = wR_interface;
    dry_tol = cfg.phys.dry_tolerance;
    num_interfaces = N + 1;
    for k_int = 1:num_interfaces
        HL_k = wL_interface(k_int, 1);
        HUL_k = wL_interface(k_int, 2);
        HR_k = wR_interface(k_int, 1);
        HUR_k = wR_interface(k_int, 2);
        z_cell_L_k = z_cell_padded(ng + k_int - 1);
        z_cell_R_k = z_cell_padded(ng + k_int);
        z_interface_eff_k = max(z_cell_L_k, z_cell_R_k);
        HL_k_recon = max(0, HL_k + z_cell_L_k - z_interface_eff_k);
        HR_k_recon = max(0, HR_k + z_cell_R_k - z_interface_eff_k);
        uL_k = 0;
        if HL_k > dry_tol
            uL_k = HUL_k / HL_k;
        end
        HUL_k_recon = HL_k_recon * uL_k;
        uR_k = 0;
        if HR_k > dry_tol
            uR_k = HUR_k / HR_k;
        end
        HUR_k_recon = HR_k_recon * uR_k;
        wL_interface_recon(k_int, 1) = HL_k_recon;
        wL_interface_recon(k_int, 2) = HUL_k_recon;
        wR_interface_recon(k_int, 1) = HR_k_recon;
        wR_interface_recon(k_int, 2) = HUR_k_recon;
    end
    % Use hydrostatically reconstructed states for flux calculation
    F_num = cfg.numFlux(wL_interface_recon, wR_interface_recon, cfg);  % [N+1, 2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Flux Divergence (Spatial Derivative)                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finite volume update: d(w_i)/dt = -(F_{i+1/2} - F_{i-1/2}) / dx
    dwdt_flux = -(F_num(2:N+1,:) - F_num(1:N,:)) / dx; % [N, 2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Source Terms (Friction, Bed Slope)                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dwdt_source = zeros(N, 2); % [N, 2] Initialize source term array

    % --- 1. Bed slope source term (S_b) ---
    % Only needed for non-flat bathymetry
    if isfield(cfg, 'bathyHandle') && ~isequal(func2str(cfg.bathyHandle), 'bathy.flat')
        % Calculate bathymetry at cell centers
        z_centers = cfg.bathyHandle(cfg, cfg.mesh.xc);  % [N, 1]
        % Calculate bathymetry at cell interfaces
        xf = [cfg.mesh.xc(1) - dx/2; (cfg.mesh.xc + dx/2)'];  % Interface positions (column vector)
        z_interfaces = cfg.bathyHandle(cfg, xf);  % [N+1, 1]
        % Calculate bed slope source term for momentum equation
        % S_b = -g * H * dz/dx
        for i = 1:N
            H_i = w(i, 1);  % Water depth at cell i
            if H_i > cfg.phys.dry_tolerance
                % Well-balanced discretization of bed slope term
                % Use central difference of cell-centered bathymetry for well-balanced source term
                z_im1_center = z_cell_padded(ng + i - 1); % Bathymetry at center of cell to the left of current cell i
                z_ip1_center = z_cell_padded(ng + i + 1); % Bathymetry at center of cell to the right of current cell i
                bathymetry_slope_at_i = (z_ip1_center - z_im1_center) / (2 * dx);
                dwdt_source(i, 2) = dwdt_source(i, 2) - g * H_i * bathymetry_slope_at_i;
            end
        end
    end

    % --- 2. Friction source term (S_f) ---
    % Apply friction if a friction model is specified
    if isfield(cfg.phys, 'friction_model') && ~isempty(cfg.phys.friction_model)
        % Call the selected friction model to get the friction term
        H = w(:, 1);  % [m] Water depth
        HU = w(:, 2); % [m^2/s] Discharge
        wet_indices = H > cfg.phys.dry_tolerance; % Indices of wet cells
        
        % Apply to momentum equation (only for wet cells)
        if any(wet_indices)
            friction_term = cfg.phys.friction_model(H(wet_indices), HU(wet_indices), g, cfg);
            dwdt_source(wet_indices, 2) = dwdt_source(wet_indices, 2) + friction_term;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Combine Flux and Source Terms, Flatten Output               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dwdt = dwdt_flux + dwdt_source; % [N, 2] Total time derivative
    dwdt_flat = [dwdt(:,1); dwdt(:,2)]; % [2N x 1] Flattened column vector

end