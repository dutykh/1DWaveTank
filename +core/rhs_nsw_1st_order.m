function dwdt_flat = rhs_nsw_1st_order(t, w_flat, cfg)

    %RHS_NSW_1ST_ORDER Computes the RHS for the 1D NSW equations using a 1st order FV scheme.
    %   dwdt_flat = RHS_NSW_1ST_ORDER(t, w_flat, cfg) calculates the time
    %   derivative of the state vector for the Nonlinear Shallow Water Equations.
    %   It uses the numerical flux specified in cfg.numFlux (e.g., @flux.FVCF)
    %   and boundary conditions specified by handles in cfg.bc.left.handle and
    %   cfg.bc.right.handle. Assumes 1st order accuracy (no reconstruction).
    %
    %   Inputs:
    %       t         - Current time (scalar).
    %       w_flat    - Flattened state vector [H1;...;HN; HU1;...;HUN] (column vector).
    %       cfg       - Configuration structure, must contain:
    %                     cfg.mesh.N  - Number of cells
    %                     cfg.mesh.dx - Cell width
    %                     cfg.phys.g  - Acceleration due to gravity
    %                     cfg.phys.Cf - Friction coefficient (Manning's n), set to 0 for no friction
    %                     cfg.numFlux - Function handle for numerical flux (e.g., @flux.FORCE)
    %                     cfg.bc.left.handle  - Function handle for left boundary condition
    %                     cfg.bc.right.handle - Function handle for right boundary condition
    %                   Optional:
    %                     cfg.bathyHandle - Function handle for bathymetry (default: flat)
    %
    %   Outputs:
    %       dwdt_flat - Flattened time derivative vector [dH/dt; dHU/dt] (column vector).
    %
    %   Author: Denys Dutykh
    %   Date:   20 April 2025

    w_flat = w_flat(:); % Ensure w_flat is a column vector

    % Extract parameters from config
    N = cfg.mesh.N;
    dx = cfg.mesh.dx;
    g = cfg.phys.g; % Gravity
    Cf = cfg.phys.Cf; % Friction coefficient (Manning's n)

    % Reshape the flattened state vector w_flat into an N x 2 array [H, HU]
    % H = w(:,1), HU = w(:,2)
    w = [w_flat(1:N), w_flat(N+1:2*N)];

    % --- Apply Boundary Conditions using function handles --- 
    % We need ghost cells to calculate fluxes at the domain boundaries (interfaces 1/2 and N+1/2).
    % For 1st order fluxes, only 1 ghost cell on each side is required.
    % Higher order reconstructions would typically require more ghost cells.
    % Determine required number of ghost cells
    % TODO: Adaptively determine num_ghost_cells based on reconstruction order if higher order is implemented.
    num_ghost_cells = 1; % Sufficient for 1st order fluxes

    w_padded = zeros(N + 2*num_ghost_cells, 2); % Array to hold interior cells + ghost cells
    w_padded(num_ghost_cells+1 : N+num_ghost_cells, :) = w; % Fill interior domain data

    % Check if BC handles are provided
    if ~isfield(cfg, 'bc') || ~isfield(cfg.bc, 'left') || ~isfield(cfg.bc.left, 'handle')
        error('Left boundary condition handle (cfg.bc.left.handle) not specified.');
    end
    if ~isfield(cfg.bc, 'right') || ~isfield(cfg.bc.right, 'handle')
        error('Right boundary condition handle (cfg.bc.right.handle) not specified.');
    end

    % Apply left boundary condition using the specified handle via feval.
    % The handle function modifies the first `num_ghost_cells` rows of w_padded.
    w_padded = feval(cfg.bc.left.handle, w_padded, t, 'left', cfg, num_ghost_cells);

    % Apply right boundary condition using the specified handle via feval.
    % The handle function modifies the last `num_ghost_cells` rows of w_padded.
    w_padded = feval(cfg.bc.right.handle, w_padded, t, 'right', cfg, num_ghost_cells);
    % --- End Boundary Conditions --- 

    % --- Calculate Numerical Fluxes at Interfaces --- 
    % Interfaces are indexed i+1/2, running from 1/2 (left boundary) to N+1/2 (right boundary).
    % There are N+1 interfaces in total.
    F_num = zeros(N+1, 2); % Array to store fluxes [F_H, F_HU] at each interface

    % Loop over all N+1 interfaces
    % For a 1st order scheme, the state left (wL) and right (wR) of interface i+1/2
    % are simply the cell-centered values from cells i and i+1, respectively.
    for i = 1:(N+1) % Loop interfaces from i=1 (interface 1/2) to i=N+1 (interface N+1/2)
        % Indices in the padded array corresponding to cells i (left) and i+1 (right)
        idxL = i + num_ghost_cells - 1; % Index of cell i in w_padded
        idxR = i + num_ghost_cells;     % Index of cell i+1 in w_padded

        wL = w_padded(idxL, :); % State [H, HU] in cell i (left of interface i+1/2)
        wR = w_padded(idxR, :); % State [H, HU] in cell i+1 (right of interface i+1/2)

        % Call the numerical flux function (handle provided in cfg.numFlux)
        % The flux function returns the numerical flux vector [F_H, F_HU] across interface i+1/2
        F_num(i,:) = cfg.numFlux(wL, wR, cfg);
    end

    % --- Calculate Spatial Derivative (Flux Divergence) --- 
    % The finite volume update for cell i is: d(w_i)/dt = -(F_{i+1/2} - F_{i-1/2}) / dx
    % This calculates the contribution to dw/dt from the flux differences.
    dwdt_flux = -(F_num(2:N+1,:) - F_num(1:N,:)) / dx; % Result is N x 2 array

    % --- Source Terms --- 
    dwdt_source = zeros(N, 2);

    % 1. Bed slope source term (S_b)
    % IMPORTANT: This implementation assumes a flat bottom (h=constant).
    % For non-flat bathymetry, a well-balanced scheme is crucial to correctly
    % handle steady states (e.g., lake at rest) and accurately compute flows.
    % Adding a naive bed slope term here can lead to spurious oscillations.
    % A proper implementation involves either modifying the numerical flux
    % (flux balancing) or adding carefully constructed source terms that
    % balance the flux gradient (source term balancing).
    % Example (naive, NOT well-balanced, requires cfg.bathyHandle):
    % if isfield(cfg, 'bathyHandle') && ~isequal(func2str(cfg.bathyHandle), 'cfg.bathy.flat')
    %     % Needs careful calculation of h at interfaces and inside cells
    %     h_interfaces = ... % Interpolate/calculate bathymetry at interfaces
    %     bed_slope_term = -g * w(:,1) .* (h_interfaces(2:N+1) - h_interfaces(1:N)) / dx; % Approx -g*H*dh/dx
    %     dwdt_source(:,2) = dwdt_source(:,2) + bed_slope_term;
    % end

    % 2. Friction source term (S_f) using Manning's formula
    % Formula: S_f = -g * n^2 * |U| * U / H^(4/3), where n = Cf
    % This term acts only on the momentum equation (dHU/dt).
    if Cf > 0
        H = w(:, 1);  % Water depth H in each cell
        HU = w(:, 2); % Discharge HU in each cell

        % Avoid division by zero and apply friction only in sufficiently wet cells
        wet_indices = H > 1e-6; % Tolerance for identifying wet cells
        U_wet = zeros(size(H)); % Velocity U = HU/H, only calculated for wet cells
        U_wet(wet_indices) = HU(wet_indices) ./ H(wet_indices);

        % Calculate Manning friction term for wet cells
        % Add small epsilon to denominator H^(7/3) for numerical stability near H=0
        friction_term = -g * Cf^2 * abs(U_wet(wet_indices)) .* HU(wet_indices) ./ (H(wet_indices).^(7/3) + 1e-10);

        % Add friction term to the momentum source term for wet cells
        dwdt_source(wet_indices, 2) = dwdt_source(wet_indices, 2) + friction_term;
    end

    % --- Combine terms --- 
    % Total time derivative is sum of flux divergence and source terms
    dwdt = dwdt_flux + dwdt_source; % Result is N x 2 array

    % --- Flatten output --- 
    % Reshape the N x 2 derivative array back into a 2N x 1 column vector
    dwdt_flat = [dwdt(:,1); dwdt(:,2)];

end