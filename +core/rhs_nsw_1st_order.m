%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +core/rhs_nsw_1st_order.m
%
% Purpose:
%   Computes the right-hand side (RHS) for the 1D Nonlinear Shallow Water (NSW)
%   equations using a first-order finite volume (FV) scheme. Handles boundary
%   conditions, numerical fluxes, and source terms (friction, bed slope).
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
%               cfg.phys.g:   [double] Acceleration due to gravity [m/s^2]
%               cfg.phys.Cf:  [double] Friction coefficient (Manning's n)
%               cfg.numFlux:  [function handle] Numerical flux function
%               cfg.bc.left.handle:  [function handle] Left BC
%               cfg.bc.right.handle: [function handle] Right BC
%             Optional:
%               bathyHandle: [function handle] Bathymetry (default: flat)
%
% Outputs:
%   dwdt_flat - [2N x 1, double] Flattened time derivative vector [dH/dt; dHU/dt].
%
% Dependencies:
%   Expects correct configuration and function handles for BCs and fluxes.
%
% References:
%   - LeVeque, R.J. (2002). Finite Volume Methods for Hyperbolic Problems.
%   - Toro, E.F. (2001). Shock-Capturing Methods for Free-Surface Shallow Flows.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
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
    Cf = cfg.phys.Cf;       % [unitless] Friction coefficient (Manning's n)

    % Reshape the flattened state vector w_flat into an N x 2 array [H, HU]
    % H = w(:,1), HU = w(:,2)
    w = [w_flat(1:N), w_flat(N+1:2*N)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ghost Cells and Boundary Conditions                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For 1st order FV, only 1 ghost cell on each side is required.
    num_ghost_cells = 1; % Sufficient for 1st order fluxes

    % Pad the domain with ghost cells (zeros for now)
    w_padded = zeros(N + 2*num_ghost_cells, 2); % [N+2, 2] array
    w_padded(num_ghost_cells+1 : N+num_ghost_cells, :) = w; % Fill interior domain data

    % --- Apply Boundary Conditions ---
    Ng = num_ghost_cells;        % Correct - use the locally defined value

    % Check for periodic BCs first
    is_periodic = isequal(cfg.bc.left.handle, @bc.periodic) && isequal(cfg.bc.right.handle, @bc.periodic);

    if is_periodic
        % Apply periodic BCs once for both sides
        w_padded = cfg.bc.left.handle(w_padded, t, 'both', cfg, Ng);
    else
        % Apply left BC
        if isfield(cfg.bc.left, 'handle') && ~isempty(cfg.bc.left.handle)
            w_padded = cfg.bc.left.handle(w_padded, t, 'left', cfg, Ng);
        else
            warning('core:rhs:NoLeftBC', 'No left boundary condition handle specified.');
        end

        % Apply right BC
        if isfield(cfg.bc.right, 'handle') && ~isempty(cfg.bc.right.handle)
            w_padded = cfg.bc.right.handle(w_padded, t, 'right', cfg, Ng);
        else
            warning('core:rhs:NoRightBC', 'No right boundary condition handle specified.');
        end
    end

    % Extract H and HU components after BCs applied
    H_padded = w_padded(:, 1);
    HU_padded = w_padded(:, 2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Numerical Flux Calculation at Cell Interfaces              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For 1st order FV, the state left (wL) and right (wR) of interface i+1/2
    % are simply the cell-centered values from cells i and i+1, respectively.
    % There are N+1 interfaces in total (from 1/2 to N+1/2).
    F_num = zeros(N+1, 2); % Preallocate (though overwritten below)
    
    % Vectorized approach: Prepare all left/right states first
    idxL = (1:(N+1)) + num_ghost_cells - 1; % Indices of cells left of each interface [1 x (N+1)]
    idxR = (1:(N+1)) + num_ghost_cells;     % Indices of cells right of each interface [1 x (N+1)]
    wL = w_padded(idxL, :);                 % [(N+1) x 2] States left of interfaces
    wR = w_padded(idxR, :);                 % [(N+1) x 2] States right of interfaces
    
    % Compute all numerical fluxes at once by calling the (assumed vectorized) flux function
    F_num = cfg.numFlux(wL, wR, cfg);       % [(N+1) x 2] 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Flux Divergence (Spatial Derivative)                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finite volume update: d(w_i)/dt = -(F_{i+1/2} - F_{i-1/2}) / dx
    dwdt_flux = -(F_num(2:N+1,:) - F_num(1:N,:)) / dx; % [N, 2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Source Terms (Friction, Bed Slope)                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dwdt_source = zeros(N, 2); % [N, 2] Initialize source term array

    % --- START: Bathymetric Source Term Calculation ---
    % Implements S_b = g * H * (dh/dx) using a cell-centered, central difference
    % approach for well-balanced shallow water equations.
    % This works for arbitrary bathymetry.

    % Prepare cell center coordinates
    % Prepare bathymetry including ghost cells (z_bc)
    x_cell = cfg.mesh.xc; % [1 x N] cell centers
    dx = cfg.mesh.dx;
    Ng = num_ghost_cells;
    % Ghost cell centers
    x_ghost_left = x_cell(1) - dx;
    x_ghost_right = x_cell(end) + dx;
    bathy_func = cfg.bathyHandle;
    z_ghost_left = bathy_func(cfg, x_ghost_left);
    z_ghost_right = bathy_func(cfg, x_ghost_right);
    z_cell = bathy_func(cfg, x_cell); z_cell = z_cell(:)';
    z_bc = [z_ghost_left, z_cell, z_ghost_right]; % [1 x (N+2)]

    % --- START: Bathymetric Source Term Calculation (Well-Balanced) ---
    % Compute source terms (bathymetry)
    S_src = zeros(2, N); % Initialize S_src, S_src(1,:) remains zero
    % z_bc must be available: bathymetry at all cells including ghost cells
    % w is [N x 2]: w(:,1) = h (water depth), w(:,2) = hu (discharge)
    % cfg.phys.g: gravity, cfg.mesh.dx: grid spacing, Ng: ghost cells, N: # cells
    for i = 1:N
        h_local = w(i, 1);
        z_im1 = z_bc(Ng + i - 1); % bathymetry at i-1
        z_ip1 = z_bc(Ng + i + 1); % bathymetry at i+1
        S_src(2, i) = -cfg.phys.g * h_local * (z_ip1 - z_im1) / (2 * cfg.mesh.dx);
    end
    % Add to source term array used in dwdt update
    dwdt_source(:,2) = dwdt_source(:,2) + S_src(2,:)';
    % --- END: Bathymetric Source Term Calculation ---

    % 2. Friction source term (S_f)
    % Apply friction if a friction model is specified
    if isfield(cfg.phys, 'friction_model') && ~isempty(cfg.phys.friction_model)
        % Call the selected friction model to get the friction term
        H = w(:, 1);  % [m] Water depth
        HU = w(:, 2); % [m^2/s] Discharge
        wet_indices = H > cfg.phys.dry_tolerance; % Indices of wet cells
        
        % Apply to momentum equation (only for wet cells)
        friction_term = cfg.phys.friction_model(H(wet_indices), HU(wet_indices), g, cfg);
        dwdt_source(wet_indices, 2) = dwdt_source(wet_indices, 2) + friction_term;
    end
    
    % Legacy friction support (for backward compatibility)
    % Only use if friction_model is not set and Cf > 0
    if (~isfield(cfg.phys, 'friction_model') || isempty(cfg.phys.friction_model)) && Cf > 0
        H = w(:, 1);  % [m] Water depth
        HU = w(:, 2); % [m^2/s] Discharge
        wet_indices = H > 1e-6; % Indices of wet cells
        U_wet = zeros(size(H)); % [m/s] Velocity, only for wet cells
        U_wet(wet_indices) = HU(wet_indices) ./ H(wet_indices);
        % Calculate Manning friction term for wet cells
        % Add small epsilon to denominator H^(7/3) for numerical stability near H=0
        friction_term = -g * Cf^2 * abs(U_wet(wet_indices)) .* HU(wet_indices) ./ (H(wet_indices).^(7/3) + 1e-10);
        dwdt_source(wet_indices, 2) = dwdt_source(wet_indices, 2) + friction_term;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Combine Flux and Source Terms, Flatten Output               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dwdt = dwdt_flux + dwdt_source; % [N, 2] Total time derivative
    dwdt_flat = [dwdt(:,1); dwdt(:,2)]; % [2N x 1] Flattened column vector

end