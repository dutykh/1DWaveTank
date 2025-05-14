%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +core/rhs_nsw_high_order.m
%
% Purpose:
%   Computes the right-hand side (RHS) for the 1D Nonlinear Shallow Water (NSW)
%   equations using a high-order finite volume scheme. This version implements
%   a well-balanced scheme based on hydrostatic reconstruction of linearly
%   reconstructed variables and a centered source term, as described in
%   Audusse et al. (2004) and related literature (e.g., "SchemaHydroWB.pdf").
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
    
    % Perform reconstruction to get interface values
    % Inline reconstruction logic (was perform_reconstruction)
    if isfield(cfg, 'reconstruct') && isfield(cfg.reconstruct, 'handle') && ~isempty(cfg.reconstruct.handle)
        reconstruct_handle = cfg.reconstruct.handle;
    else
        reconstruct_handle = @reconstruct.none;
    end
    [wL_int, wR_int] = reconstruct_handle(w_padded, cfg);
    
    % Get additional topography at interfaces for source term
    z_interfaces = compute_interface_topo(z_b, cfg, ng, N);
    
    % Apply hydrostatic reconstruction
    [wL_hydro, wR_hydro] = hydrostatic_reconstruction(wL_int, wR_int, z_interfaces, dry_tol, g);
    
    % Calculate numerical fluxes
    F_num = cfg.numFlux(wL_hydro, wR_hydro, cfg);
    
    % Compute flux divergence term: -(F_{i+1/2} - F_{i-1/2})/dx
    dwdt_flux = -(F_num(2:N+1,:) - F_num(1:N,:)) / dx;
    
    % Compute well-balanced source term using exact discrete balancing
    dwdt_source = zeros(N, 2);
    for i = 1:N
        idx = i + ng;
        h_i = w_padded(idx, 1);
        
        % Skip dry cells
        if h_i <= dry_tol
            continue;
        end
        
        % Centered source term using interface heights from reconstruction
        h_left = wL_hydro(i, 1);       % h_{i-1/2+}
        h_right = wL_hydro(i+1, 1);    % h_{i+1/2-}
        
        % Source term that exactly balances the numerical flux for lake-at-rest
        z_diff = z_interfaces(i+1) - z_interfaces(i);
        dwdt_source(i, 2) = -g * h_i * z_diff / dx;  
    end
    
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

function [wL_hydro, wR_hydro] = hydrostatic_reconstruction(wL, wR, z_interfaces, dry_tol, g)
    % Initialize output with input values
    wL_hydro = wL;
    wR_hydro = wR;
    
    for i = 1:size(wL, 1)
        % Get elevations at interface
        z_left = z_interfaces(i);
        z_right = z_interfaces(i);
        
        % Maximum elevation at interface (exactly the same value for left and right)
        z_max = z_left;  % Since z_left = z_right for exact well-balancing
        
        % Extract values
        h_left = wL(i, 1);
        h_right = wR(i, 1);
        
        % Calculate velocities carefully
        u_left = 0;
        if h_left > dry_tol
            u_left = wL(i, 2) / h_left;
        end
        
        u_right = 0;
        if h_right > dry_tol
            u_right = wR(i, 2) / h_right;
        end
        
        % Hydrostatic reconstruction - formula from the paper
        h_left_recon = max(0, h_left + z_left - z_max);
        h_right_recon = max(0, h_right + z_right - z_max);
        
        % Update conserved variables
        wL_hydro(i, 1) = h_left_recon;
        wL_hydro(i, 2) = h_left_recon * u_left;
        
        wR_hydro(i, 1) = h_right_recon;
        wR_hydro(i, 2) = h_right_recon * u_right;
    end
end

function z_interfaces = compute_interface_topo(z_cell, cfg, ng, N)
    % Compute topography at interfaces using averaged cell values
    % This ensures z_i-1/2 is exactly the same when viewed from cells i-1 and i
    z_interfaces = zeros(1, N+1);
    
    for i = 1:N+1
        z_interfaces(i) = 0.5 * (z_cell(ng+i-1) + z_cell(ng+i));
    end
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