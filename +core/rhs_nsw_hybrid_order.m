function dwdt_flat = rhs_nsw_hybrid_order(t, w_flat, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % +core/rhs_nsw_hybrid_order.m
    %
    % Purpose:
    %   Computes the right-hand side (RHS) for the 1D Nonlinear Shallow Water (NSW)
    %   equations using a hybrid approach:
    %   - Higher-order reconstruction for flux computation
    %   - First-order hydrostatic reconstruction for source terms
    %   This ensures well-balanced properties while maintaining higher accuracy.
    %
    % Syntax:
    %   dwdt_flat = rhs_nsw_hybrid_order(t, w_flat, cfg)
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
    % Date:   May 15, 2025
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Input Preparation and Parameter Extraction                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        w_flat = w_flat(:); % Ensure w_flat is a column vector
    
        N = cfg.mesh.N;             % Number of spatial cells
        dx = cfg.mesh.dx;           % Cell width [m]
        g = cfg.phys.g;             % Gravity [m/s^2]
        dry_tol = cfg.phys.dry_tolerance; % Threshold for dry cells
    
        % Reshape the flattened state vector
        w = [w_flat(1:N), w_flat(N+1:2*N)]; % [H, HU]
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Ghost Cells and Boundary Conditions                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ng = cfg.bc.num_ghost_cells; % Number of ghost cells (for higher-order)
    
        w_padded = zeros(N + 2*ng, 2);
        w_padded(ng+1 : N+ng, :) = w;
    
        % Apply boundary conditions
        is_periodic = isequal(cfg.bc.left.handle, @bc.periodic) && isequal(cfg.bc.right.handle, @bc.periodic);
        if is_periodic
            w_padded = cfg.bc.left.handle(w_padded, t, 'both', cfg, ng);
        else
            if isfield(cfg.bc.left, 'handle') && ~isempty(cfg.bc.left.handle)
                w_padded = cfg.bc.left.handle(w_padded, t, 'left', cfg, ng);
            else
                warning('core:rhs:NoLeftBC', 'No left boundary condition handle specified.');
            end
            if isfield(cfg.bc.right, 'handle') && ~isempty(cfg.bc.right.handle)
                w_padded = cfg.bc.right.handle(w_padded, t, 'right', cfg, ng);
            else
                warning('core:rhs:NoRightBC', 'No right boundary condition handle specified.');
            end
        end
    
        % Extract H and HU, calculate U for wet cells
        H_padded = w_padded(:,1);
        HU_padded = w_padded(:,2);
        U_padded = zeros(size(H_padded));
        wet_idx_padded = H_padded > dry_tol;
        U_padded(wet_idx_padded) = HU_padded(wet_idx_padded) ./ H_padded(wet_idx_padded);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Bathymetry at Cell Centers (including ghost cells)         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x_physical_centers = cfg.mesh.xc; % 1xN
        x_ghost_left_centers = zeros(1, ng);
        for k_ghost = 1:ng
            x_ghost_left_centers(k_ghost) = x_physical_centers(1) - (ng - k_ghost + 1) * dx;
        end
        x_ghost_right_centers = zeros(1, ng);
        for k_ghost = 1:ng
            x_ghost_right_centers(k_ghost) = x_physical_centers(end) + k_ghost * dx;
        end
        x_all_cell_centers = [x_ghost_left_centers, x_physical_centers, x_ghost_right_centers];
        z_bc = cfg.bathyHandle(cfg, x_all_cell_centers); % Bathymetry at all cell centers
        z_bc = z_bc(:)'; % Ensure row vector
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FLUX COMPUTATION: Higher-Order Approach                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use higher-order reconstruction to get interface values
        [wL_ho, wR_ho] = cfg.reconstruct.handle(w_padded, cfg);
        
        % Apply hydrostatic reconstruction to higher-order values
        wL_recon = zeros(N+1, 2);
        wR_recon = zeros(N+1, 2);
    
        for k_int = 1:N+1 % Loop over interfaces
            % Get cell indices adjacent to interface
            idx_i = ng + k_int - 1;  % Left cell
            idx_ip1 = idx_i + 1;     % Right cell
    
            % Extract bathymetry values
            z_i = z_bc(idx_i);
            z_ip1 = z_bc(idx_ip1);
            max_z_interface = max(z_i, z_ip1);
    
            % Get reconstructed states
            h_L = wL_ho(k_int, 1);
            hu_L = wL_ho(k_int, 2);
            h_R = wR_ho(k_int, 1);
            hu_R = wR_ho(k_int, 2);
    
            % Calculate velocities
            u_L = 0;
            if h_L > dry_tol
                u_L = hu_L / h_L;
            end
            
            u_R = 0;
            if h_R > dry_tol
                u_R = hu_R / h_R;
            end
    
            % Apply hydrostatic reconstruction to higher-order values
            h_G_k = max(0, h_L + z_i - max_z_interface);
            hu_G_k = h_G_k * u_L;
            wL_recon(k_int, :) = [h_G_k, hu_G_k];
    
            h_D_k = max(0, h_R + z_ip1 - max_z_interface);
            hu_D_k = h_D_k * u_R;
            wR_recon(k_int, :) = [h_D_k, hu_D_k];
        end
        
        % Calculate numerical fluxes
        F_num = cfg.numFlux(wL_recon, wR_recon, cfg);
    
        % Compute flux divergence
        dwdt_flux = -(F_num(2:N+1,:) - F_num(1:N,:)) / dx;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SOURCE TERMS: First-Order Approach                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dwdt_source = zeros(N, 2);
    
        % --- Well-Balanced Bed Slope Source Term (First-Order) ---
        S_b_momentum = zeros(N,1);
        for i = 1:N % Loop over physical cells
            idx_cell_i = ng + i; % Padded index for current cell i
    
            h_i_cell = H_padded(idx_cell_i);
            z_i_cell = z_bc(idx_cell_i);
            
            % Skip dry cells
            if h_i_cell <= dry_tol
                continue;
            end
            
            % Bathymetry of neighboring cells
            z_im1_cell = z_bc(idx_cell_i - 1);
            z_ip1_cell = z_bc(idx_cell_i + 1);
    
            % First-order hydrostatic reconstruction at right interface (i+1/2)
            h_star_R_i = max(0, h_i_cell + z_i_cell - max(z_i_cell, z_ip1_cell));
            
            % First-order hydrostatic reconstruction at left interface (i-1/2)
            h_star_L_i = max(0, h_i_cell + z_i_cell - max(z_im1_cell, z_i_cell));
            
            % Source term formula from first-order scheme
            S_b_momentum(i) = (g / (2 * dx)) * (h_star_R_i^2 - h_star_L_i^2);
        end
        dwdt_source(:, 2) = dwdt_source(:, 2) + S_b_momentum;
    
        % --- Friction Source Term (S_f) ---
        if isfield(cfg.phys, 'friction_model') && ~isempty(cfg.phys.friction_model)
            H_physical = w(:, 1);
            HU_physical = w(:, 2);
            wet_indices_physical = H_physical > dry_tol;
            
            if any(wet_indices_physical)
                friction_term_values = cfg.phys.friction_model(H_physical(wet_indices_physical), HU_physical(wet_indices_physical), g, cfg);
                dwdt_source(wet_indices_physical, 2) = dwdt_source(wet_indices_physical, 2) + friction_term_values(:);
            end
        elseif isfield(cfg.phys, 'Cf') && cfg.phys.Cf > 0 % Legacy friction
            H_physical = w(:, 1);
            HU_physical = w(:, 2);
            U_physical = zeros(N,1);
            wet_indices_physical = H_physical > dry_tol;
            U_physical(wet_indices_physical) = HU_physical(wet_indices_physical) ./ H_physical(wet_indices_physical);
            
            friction_legacy = zeros(N,1);
            friction_legacy(wet_indices_physical) = -g * cfg.phys.Cf^2 .* abs(U_physical(wet_indices_physical)) .* U_physical(wet_indices_physical) ./ (H_physical(wet_indices_physical).^(4/3) + cfg.numerics.epsilon);
    
            dwdt_source(wet_indices_physical, 2) = dwdt_source(wet_indices_physical, 2) + friction_legacy(wet_indices_physical);
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Combine Flux and Source Terms, Flatten Output               %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dwdt = dwdt_flux + dwdt_source;
        dwdt_flat = [dwdt(:,1); dwdt(:,2)];
    
    end