function dwdt_flat = rhs_nsw_1st_order(t, w_flat, cfg)

    %RHS_NSW_1ST_ORDER Computes the RHS for the 1D NSW equations using a 1st order FV scheme.
    %   dwdt_flat = RHS_NSW_1ST_ORDER(t, w_flat, cfg) calculates the time
    %   derivative of the state vector for the Nonlinear Shallow Water Equations.
    %   It uses the numerical flux specified in cfg.numFlux (e.g., @flux.FVCF)
    %   and boundary conditions from cfg.bcL and cfg.bcR. Assumes 1st order
    %   accuracy (no reconstruction).
    %
    %   Inputs:
    %       t         - Current time.
    %       w_flat    - Flattened state vector [H1;...;HN; HU1;...;HUN].
    %       cfg       - Configuration structure.
    %
    %   Outputs:
    %       dwdt_flat - Flattened time derivative vector [dH/dt; dHU/dt].
    
    % Extract parameters from config
    N = cfg.mesh.N;
    dx = cfg.mesh.dx;
    g = cfg.param.g;
    Cf = cfg.param.Cf; % Friction coefficient (will be 0 for this case)
        
    % Reshape the flattened state vector w_flat into N x 2 array [H, HU]
    w = [w_flat(1:N), w_flat(N+1:2*N)];
        
    % --- Apply Boundary Conditions ---
    % We need ghost cells for flux calculation at domain boundaries.
    % For 1st order, we need 1 ghost cell on each side.
    num_ghost_cells = 1;
    w_padded = zeros(N + 2*num_ghost_cells, 2); % Array with ghost cells
    w_padded(num_ghost_cells+1 : N+num_ghost_cells, :) = w; % Fill interior
        
    % Apply left boundary condition (function specified in cfg.bcL)
    % The BC function should fill w_padded(1:num_ghost_cells, :)
    w_padded = cfg.bcL(w_padded, t, 'left', cfg, num_ghost_cells);
        
    % Apply right boundary condition (function specified in cfg.bcR)
    % The BC function should fill w_padded(N+num_ghost_cells+1 : N+2*num_ghost_cells, :)
    w_padded = cfg.bcR(w_padded, t, 'right', cfg, num_ghost_cells);
        
    % --- Calculate Numerical Fluxes at Interfaces ---
    % Interfaces are indexed i+1/2, from 1/2 to N+1/2.
    % There are N+1 interfaces.
    F_num = zeros(N+1, 2); % Array to store fluxes at interfaces
        
    % For 1st order scheme: wL at i+1/2 is w_i, wR at i+1/2 is w_{i+1}
    % We use the padded array indices: w_i corresponds to w_padded(i+num_ghost_cells,:)
    for i = 1:(N+1) % Loop over interfaces i=1/2 to N+1/2
        idxL = i + num_ghost_cells - 1; % Index of cell i (left of interface i+1/2)
        idxR = i + num_ghost_cells;     % Index of cell i+1 (right of interface i+1/2)
            
        wL = w_padded(idxL, :);
        wR = w_padded(idxR, :);
            
        % Call the numerical flux function specified in the config
        F_num(i,:) = cfg.numFlux(wL, wR, cfg);
    end
        
    % --- Calculate Spatial Derivative (Flux Divergence) ---
    % dw_i/dt = (F_{i-1/2} - F_{i+1/2}) / dx
    dwdt_flux = (F_num(1:N,:) - F_num(2:N+1,:)) / dx;
        
    % --- Source Terms ---
    dwdt_source = zeros(N, 2);
        
    % 1. Bed slope source term (zero for flat bottom)
    %    If bathymetry h(x) is not constant, add:
    %    h_padded = cfg.bathyHandle(cfg.mesh.x_edge_padded); % Need padded bathy
    %    h_x = (h_padded(idxR) - h_padded(idxL)) / dx; % Approx slope at cell i
    %    dwdt_source(:,2) = dwdt_source(:,2) - g * w(:,1) .* h_x; % -g*H*dh/dx
        
    % 2. Friction source term (using Manning formula, Cf is Manning's n)
    if Cf > 0
        H = w(:, 1);
        HU = w(:, 2);
        wet_indices = H > 1e-6; % Apply friction only in wet cells
        U_wet = HU(wet_indices) ./ H(wet_indices);
        friction_term = -g * Cf^2 * abs(U_wet) .* HU(wet_indices) ./ (H(wet_indices).^(7/3)); % Manning friction
        dwdt_source(wet_indices, 2) = dwdt_source(wet_indices, 2) + friction_term;
    end
        
    % --- Combine terms ---
    dwdt = dwdt_flux + dwdt_source;
        
    % --- Flatten output ---
    dwdt_flat = [dwdt(:,1); dwdt(:,2)];

end