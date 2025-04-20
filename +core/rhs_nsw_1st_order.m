function dwdt_flat = rhs_nsw_1st_order(t, w_flat, cfg)

    %RHS_NSW_1ST_ORDER Computes the RHS for the 1D NSW equations using a 1st order FV scheme.
    %   dwdt_flat = RHS_NSW_1ST_ORDER(t, w_flat, cfg) calculates the time
    %   derivative of the state vector for the Nonlinear Shallow Water Equations.
    %   It uses the numerical flux specified in cfg.numFlux (e.g., @flux.FVCF)
    %   and boundary conditions specified by handles in cfg.bc.left.handle and
    %   cfg.bc.right.handle. Assumes 1st order accuracy (no reconstruction).
    %
    %   Inputs:
    %       t         - Current time.
    %       w_flat    - Flattened state vector [H1;...;HN; HU1;...;HUN].
    %       cfg       - Configuration structure.
    %
    %   Outputs:
    %       dwdt_flat - Flattened time derivative vector [dH/dt; dHU/dt].

    w_flat = w_flat(:); % Ensure w_flat is a column vector

    % Extract parameters from config
    N = cfg.mesh.N;
    dx = cfg.mesh.dx;
    g = cfg.phys.g; % Corrected path
    Cf = cfg.phys.Cf; % Corrected path

    % Reshape the flattened state vector w_flat into N x 2 array [H, HU]
    w = [w_flat(1:N), w_flat(N+1:2*N)];

    % --- Apply Boundary Conditions using function handles --- 
    % We need ghost cells for flux calculation at domain boundaries.
    % For 1st order, we need 1 ghost cell on each side. Higher order methods
    % might require more ghost cells, handled by the BC function itself.
    % Determine required number of ghost cells (minimum 1 for 1st order flux)
    % TODO: Adaptively determine num_ghost_cells based on scheme/reconstruction order
    num_ghost_cells = 1;

    w_padded = zeros(N + 2*num_ghost_cells, 2); % Array with ghost cells
    w_padded(num_ghost_cells+1 : N+num_ghost_cells, :) = w; % Fill interior

    % --- Apply boundary conditions using handles from cfg ---
    if ~isfield(cfg, 'bc') || ~isfield(cfg.bc, 'left') || ~isfield(cfg.bc.left, 'handle')
        error('Left boundary condition handle (cfg.bc.left.handle) not specified.');
    end
    if ~isfield(cfg.bc, 'right') || ~isfield(cfg.bc.right, 'handle')
        error('Right boundary condition handle (cfg.bc.right.handle) not specified.');
    end

    % Apply left boundary condition using the specified handle
    w_padded = feval(cfg.bc.left.handle, w_padded, t, 'left', cfg, num_ghost_cells);

    % Apply right boundary condition using the specified handle
    w_padded = feval(cfg.bc.right.handle, w_padded, t, 'right', cfg, num_ghost_cells);
    % --- End Boundary Conditions ---

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
    % dw_i/dt = -(F_{i+1/2} - F_{i-1/2}) / dx
    dwdt_flux = -(F_num(2:N+1,:) - F_num(1:N,:)) / dx; % Note the sign convention change

    % --- Source Terms ---
    dwdt_source = zeros(N, 2);

    % 1. Bed slope source term (currently assumes flat bottom, h=const)
    %    If bathymetry h(x) is not constant, need a well-balanced scheme.
    %    For non-flat bathy, this term should be added carefully, potentially
    %    modifying the numerical flux or adding source term balancing.
    %    Example (naive, not well-balanced):
    %    if ~isequal(cfg.bathyHandle, @cfg.bathy.flat)
    %        xc_padded_centers = [cfg.mesh.xc(1)-dx*num_ghost_cells:dx:cfg.mesh.xc(1)-dx, ...
    %                             cfg.mesh.xc, ...
    %                             cfg.mesh.xc(end)+dx:dx:cfg.mesh.xc(end)+dx*num_ghost_cells]'; % Approximate centers of ghost cells
    %        h_padded = cfg.bathyHandle(xc_padded_centers, cfg);
    %        h_interfaces = (h_padded(num_ghost_cells+1:N+num_ghost_cells) + h_padded(num_ghost_cells+2:N+num_ghost_cells+1))/2; % h at i+1/2
    %        bed_slope_term = -g * w(:,1) .* (h_interfaces(2:N+1) - h_interfaces(1:N)) / dx; % Term approx -g*H*dh/dx
    %        dwdt_source(:,2) = dwdt_source(:,2) + bed_slope_term;
    %    end

    % 2. Friction source term (using Manning formula if Cf > 0)
    if isfield(cfg.phys,'Cf') && cfg.phys.Cf > 0
        Cf = cfg.phys.Cf; % Manning's n
        H = w(:, 1);
        HU = w(:, 2);
        wet_indices = H > 1e-6; % Apply friction only in sufficiently wet cells
        U_wet = zeros(size(H));
        U_wet(wet_indices) = HU(wet_indices) ./ H(wet_indices);
        % Manning friction term: -g * n^2 * |U| * U / H^(4/3) = -g * n^2 * |U| * HU / H^(7/3)
        friction_term = -g * Cf^2 * abs(U_wet(wet_indices)) .* HU(wet_indices) ./ (H(wet_indices).^(7/3) + 1e-10); % Added epsilon for stability
        dwdt_source(wet_indices, 2) = dwdt_source(wet_indices, 2) + friction_term;
    end

    % --- Combine terms ---
    dwdt = dwdt_flux + dwdt_source;

    % --- Flatten output ---
    dwdt_flat = [dwdt(:,1); dwdt(:,2)];

end