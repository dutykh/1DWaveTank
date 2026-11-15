function dwdt_flat = rhs_nsw_hybrid_order(t, w_flat, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % +core/rhs_nsw_hybrid_order.m
    %
    % Purpose:
    %   Computes the right-hand side (RHS) for the 1D Nonlinear Shallow Water (NSW)
    %   equations by the second-order well-balanced finite volume scheme of
    %   Audusse, Bouchut, Bristeau, Klein and Perthame (2004, Section 3), as
    %   rhs_nsw_high_order does: the depth h, the free surface h + z and the
    %   discharge hu are reconstructed at the edges of the cells by the
    %   reconstruction of the configuration, the hydrostatic reconstruction
    %   gives the two states passed to the numerical flux at each interface,
    %   and the bed slope is the one that matches these states, the sum of the
    %   hydrostatic terms at the interfaces and of the centred term
    %   g (h_l + h_r)/2 (z_l - z_r) of each cell
    %   (core.hydrostatic_reconstruction_2nd_order). A lake at rest is an exact
    %   steady state over any bottom, at a shoreline, next to a wall and at an
    %   open or wave-making end (bc.ghost_bathymetry and bc.open give the ghost
    %   cells the free surface of the end cell), the depth stays nonnegative with a limited
    %   reconstruction (MUSCL), and the scheme is second-order accurate on
    %   smooth solutions. This right-hand side also accepts the friction
    %   coefficient cfg.phys.Cf when no friction law is given.
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
        % The ghost cells take the bottom that matches the state of their boundary
        % condition (mirror image at a wall, copy of the end cell at a wave-making
        % end, periodic image, bottom function at an open end, whose ghost depth
        % bc.open sets from the free surface), so that their free surface is that
        % of the water they stand for
        z_bc = bc.ghost_bathymetry(z_bc, cfg, ng);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Second-Order Hydrostatic Reconstruction and Fluxes         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Depth, free surface and discharge reconstructed at the edges of the
        % cells, states U_{i+1/2-}, U_{i+1/2+} of the hydrostatic reconstruction,
        % and the bed slope that balances their fluxes at rest
        [wL_recon, wR_recon, S_bed] = core.hydrostatic_reconstruction_2nd_order(w_padded, z_bc, cfg);

        % Calculate numerical fluxes
        F_num = cfg.numFlux(wL_recon, wR_recon, cfg);

        % Compute flux divergence
        dwdt_flux = -(F_num(2:N+1,:) - F_num(1:N,:)) / dx;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SOURCE TERMS                                                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dwdt_source = zeros(N, 2);

        % --- Bed slope of the second-order hydrostatic reconstruction ---
        dwdt_source(:, 2) = dwdt_source(:, 2) + S_bed;
    
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