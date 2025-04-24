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
    % Numerical Flux Calculation at Cell Interfaces              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate numerical fluxes using the reconstructed states
    F_num = cfg.numFlux(wL_interface, wR_interface, cfg);  % [N+1, 2]

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
        h_centers = cfg.bathyHandle(cfg.mesh.xc, cfg);  % [N, 1]
        
        % Calculate bathymetry at cell interfaces
        xf = [cfg.mesh.xc(1) - dx/2; cfg.mesh.xc + dx/2];  % Interface positions
        h_interfaces = cfg.bathyHandle(xf, cfg);  % [N+1, 1]
        
        % Calculate bed slope source term for momentum equation
        % S_b = -g * H * dh/dx
        for i = 1:N
            H_i = w(i, 1);  % Water depth at cell i
            if H_i > cfg.phys.dry_tolerance
                % Well-balanced discretization of bed slope term
                dwdt_source(i, 2) = dwdt_source(i, 2) - g * H_i * (h_interfaces(i+1) - h_interfaces(i)) / dx;
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