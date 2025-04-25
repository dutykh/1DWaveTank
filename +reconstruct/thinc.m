function [wL_interface, wR_interface] = thinc(w_padded, cfg)
% THINC - Tangent of Hyperbola for Interface Capturing reconstruction
%
% Purpose:
%   Performs reconstruction using the THINC (Tangent of Hyperbola for 
%   Interface Capturing) method. This method is particularly effective for
%   problems with sharp discontinuities, as it maintains sharp interfaces
%   without oscillations.
%
% Syntax:
%   [wL_interface, wR_interface] = thinc(w_padded, cfg)
%
% Inputs:
%   w_padded - [(N+2*ng) x 2, double] Padded state array including ghost cells.
%              w_padded(:,1) = water depth [m],
%              w_padded(:,2) = discharge [m^2/s]
%   cfg      - [struct] Configuration structure containing:
%              cfg.mesh.N: Number of cells
%              cfg.bc.num_ghost_cells: Number of ghost cells
%              cfg.phys.dry_tolerance: Threshold for dry cells
%              cfg.reconstruct.thinc_beta: (Optional) Steepness parameter β (default: 1.5)
%
% Outputs:
%   wL_interface - [(N+1) x 2, double] Left states at all interfaces
%   wR_interface - [(N+1) x 2, double] Right states at all interfaces
%
% References:
%   - Xiao, F., Honma, Y., & Kono, T. (2005). A simple algebraic interface
%     capturing scheme using hyperbolic tangent function. International
%     Journal for Numerical Methods in Fluids, 48(9), 1023-1040.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: April 26, 2025

    % Extract parameters
    ng = cfg.bc.num_ghost_cells;  % Number of ghost cells
    N = cfg.mesh.N;               % Number of physical cells
    g = cfg.phys.g;               % Gravitational acceleration
    dry_tol = cfg.phys.dry_tolerance; % Threshold for dry cells
    
    % THINC parameter β (controls steepness of interface)
    if isfield(cfg, 'reconstruct') && isfield(cfg.reconstruct, 'thinc_beta')
        beta = cfg.reconstruct.thinc_beta;
    else
        beta = 1.5; % Default value (typical range: 1.0-2.0)
    end
    
    % Check for characteristic reconstruction flag
    use_characteristic = isfield(cfg, 'reconstruct') && isfield(cfg.reconstruct, 'characteristic') && cfg.reconstruct.characteristic;
    
    % Initialize arrays for interface values
    wL_interface = zeros(N+1, 2);
    wR_interface = zeros(N+1, 2);
    
    if use_characteristic
        % --- Characteristic-based Reconstruction --- 
        for i = 1:N+1
            idx_L = i + ng - 1;  % Left cell index (i)
            idx_R = idx_L + 1;   % Right cell index (i+1)
            
            % Check for valid cells with sufficient stencil
            if idx_L < 2 || idx_R > size(w_padded, 1)-1
                % Near boundary - use first-order
                wL_interface(i, :) = w_padded(idx_L, :);
                wR_interface(i, :) = w_padded(idx_R, :);
                continue;
            end
            
            % Use Roe average state for transformation matrix
            wL_local = w_padded(idx_L, :);   % State in cell i
            wR_local = w_padded(idx_R, :); % State in cell i+1
            
            % Check for dry cells
            if wL_local(1) < dry_tol || wR_local(1) < dry_tol
                % Fallback to first-order if cells are dry
                wL_interface(i, :) = w_padded(idx_L, :);
                wR_interface(i, :) = w_padded(idx_R, :);
                continue;
            end
            
            % Calculate Roe-averaged state
            H_roe = 0.5 * (wL_local(1) + wR_local(1)); % Simple average for depth
            
            % Compute Roe-averaged velocity
            HL_sqrt = sqrt(wL_local(1));
            HR_sqrt = sqrt(wR_local(1));
            UL = wL_local(2) / wL_local(1);
            UR = wR_local(2) / wR_local(1);
            U_roe = (HL_sqrt * UL + HR_sqrt * UR) / (HL_sqrt + HR_sqrt);
            
            % Wave speed at Roe state
            C_roe = sqrt(g * H_roe);
            
            % Eigenvectors at Roe state (for shallow water equations)
            % R = [1, 1; U_roe-C_roe, U_roe+C_roe]  (right eigenvectors as columns)
            % L = 1/(2*C_roe) * [U_roe+C_roe, -1; -U_roe+C_roe, 1]  (left eigenvectors as rows)
            R = [1, 1; U_roe-C_roe, U_roe+C_roe];
            inv_2C = 1.0 / (2.0 * C_roe + 1e-12); % Add epsilon to avoid division by zero
            L = inv_2C * [U_roe+C_roe, -1; -U_roe+C_roe, 1];
            
            % --- Left cell THINC reconstruction ---
            % Get stencil in primitive variables {i-2, i-1, i, i+1, i+2}
            N = size(w_padded, 1);
            stencil_L = [idx_L-2, idx_L-1, idx_L, idx_L+1, idx_L+2];
            stencil_L = max(min(stencil_L, N), 1); % Clamp indices to [1, N]
            prim_stencil_L = zeros(5, 2);
            
            % Convert to primitive variables (H, U)
            for j = 1:5
                prim_stencil_L(j, 1) = w_padded(stencil_L(j), 1); % H
                if prim_stencil_L(j, 1) > dry_tol
                    prim_stencil_L(j, 2) = w_padded(stencil_L(j), 2) / prim_stencil_L(j, 1); % U
                else
                    prim_stencil_L(j, 2) = 0;
                end
            end
            
            % Project onto characteristic variables
            char_stencil_L = zeros(5, 2);
            for j = 1:5
                char_stencil_L(j, :) = (L * prim_stencil_L(j, :)')';
            end
            
            % Apply THINC to each characteristic variable
            char_wL = zeros(1, 2);
            for var_idx = 1:2
                q = char_stencil_L(:, var_idx);
                
                % Use center three points for min/max
                idx_center = 3; % Center cell (i) in stencil
                q_min = min([q(idx_center-1), q(idx_center), q(idx_center+1)]);
                q_max = max([q(idx_center-1), q(idx_center), q(idx_center+1)]);
                
                % Skip THINC if the variation is too small
                if abs(q_max - q_min) < 1e-10
                    char_wL(var_idx) = q(idx_center);
                else
                    % Calculate parameter θ
                    delta = q(idx_center+1) - q(idx_center-1);
                    if abs(delta) < 1e-10
                        char_wL(var_idx) = q(idx_center);
                    else
                        theta = (2 * q(idx_center) - q_min - q_max) / (q_max - q_min);
                        B = theta * beta;
                        tanh_term = tanh(B + 0.5 * beta);
                        char_wL(var_idx) = q_min + (q_max - q_min) * (0.5 + 0.5 * tanh_term);
                    end
                end
            end
            
            % Project back to primitive variables
            prim_wL = (R * char_wL')';
            
            % Convert to conservative variables
            wL_interface(i, 1) = prim_wL(1);
            wL_interface(i, 2) = prim_wL(1) * prim_wL(2);
            
            % --- Right cell THINC reconstruction ---
            % Get stencil in primitive variables {i-1, i, i+1, i+2, i+3}
            N = size(w_padded, 1);
            stencil_R = [idx_R-2, idx_R-1, idx_R, idx_R+1, idx_R+2];
            stencil_R = max(min(stencil_R, N), 1); % Clamp indices to [1, N]
            prim_stencil_R = zeros(5, 2);
            
            % Convert to primitive variables
            for j = 1:5
                prim_stencil_R(j, 1) = w_padded(stencil_R(j), 1); % H
                if prim_stencil_R(j, 1) > dry_tol
                    prim_stencil_R(j, 2) = w_padded(stencil_R(j), 2) / prim_stencil_R(j, 1); % U
                else
                    prim_stencil_R(j, 2) = 0;
                end
            end
            
            % Project onto characteristic variables
            char_stencil_R = zeros(5, 2);
            for j = 1:5
                char_stencil_R(j, :) = (L * prim_stencil_R(j, :)')';
            end
            
            % Apply THINC to each characteristic variable
            char_wR = zeros(1, 2);
            for var_idx = 1:2
                q = char_stencil_R(:, var_idx);
                
                % Use center three points for min/max
                idx_center = 3; % Center cell (i+1) in stencil
                q_min = min([q(idx_center-1), q(idx_center), q(idx_center+1)]);
                q_max = max([q(idx_center-1), q(idx_center), q(idx_center+1)]);
                
                % Skip THINC if the variation is too small
                if abs(q_max - q_min) < 1e-10
                    char_wR(var_idx) = q(idx_center);
                else
                    % Calculate parameter θ
                    delta = q(idx_center+1) - q(idx_center-1);
                    if abs(delta) < 1e-10
                        char_wR(var_idx) = q(idx_center);
                    else
                        theta = (2 * q(idx_center) - q_min - q_max) / (q_max - q_min);
                        B = theta * beta;
                        tanh_term = tanh(B - 0.5 * beta);
                        char_wR(var_idx) = q_min + (q_max - q_min) * (0.5 + 0.5 * tanh_term);
                    end
                end
            end
            
            % Project back to primitive variables
            prim_wR = (R * char_wR')';
            
            % Convert to conservative variables
            wR_interface(i, 1) = prim_wR(1);
            wR_interface(i, 2) = prim_wR(1) * prim_wR(2);
        end
    else
        % --- Component-wise Reconstruction --- 
        % Loop over state variables (H and HU)
        for var_idx = 1:2
            % Extract the component
            q = w_padded(:, var_idx);
            
            % THINC reconstruction for this variable
            for i = 1:N+1
                % Left and right cells adjacent to the interface
                idx_L = i + ng - 1;  % Left cell index (i)
                idx_R = idx_L + 1;   % Right cell index (i+1)
                
                % Check if we have enough cells for valid reconstruction
                if idx_L < 2 || idx_R > length(q)-1
                    % Near boundary - use first-order (piecewise constant)
                    wL_interface(i, var_idx) = q(idx_L);
                    wR_interface(i, var_idx) = q(idx_R);
                    continue;
                end
                
                % --- Left cell reconstruction (for right face) ---
                % Determine local min and max for bounds
                q_min_L = min([q(idx_L-1), q(idx_L), q(idx_L+1)]);
                q_max_L = max([q(idx_L-1), q(idx_L), q(idx_L+1)]);
                
                % Skip THINC if the variation is too small (nearly constant)
                if abs(q_max_L - q_min_L) < 1e-10
                    wL_interface(i, var_idx) = q(idx_L);
                else
                    % Calculate the parameter θ which determines interface location
                    delta_L = q(idx_L+1) - q(idx_L-1);
                    if abs(delta_L) < 1e-10
                        % No clear gradient direction, use cell center value
                        wL_interface(i, var_idx) = q(idx_L);
                    else
                        % Determine interface location based on upwind/downwind values
                        theta_L = (2 * q(idx_L) - q_min_L - q_max_L) / (q_max_L - q_min_L);
                        
                        % Convert theta to interface location parameter
                        B_L = theta_L * beta;
                        
                        % THINC reconstruction at the right face (x_{i+1/2})
                        % q(x) = q_min + (q_max - q_min) * (0.5 + 0.5*tanh(beta*(x-x0)/dx))
                        % At x_{i+1/2}, x-x0 = 0.5*dx
                        tanh_term = tanh(B_L + 0.5 * beta);
                        wL_interface(i, var_idx) = q_min_L + (q_max_L - q_min_L) * (0.5 + 0.5 * tanh_term);
                    end
                end
                
                % --- Right cell reconstruction (for left face) ---
                % Determine local min and max for bounds
                q_min_R = min([q(idx_R-1), q(idx_R), q(idx_R+1)]);
                q_max_R = max([q(idx_R-1), q(idx_R), q(idx_R+1)]);
                
                % Skip THINC if the variation is too small (nearly constant)
                if abs(q_max_R - q_min_R) < 1e-10
                    wR_interface(i, var_idx) = q(idx_R);
                else
                    % Calculate the parameter θ which determines interface location
                    delta_R = q(idx_R+1) - q(idx_R-1);
                    if abs(delta_R) < 1e-10
                        % No clear gradient direction, use cell center value
                        wR_interface(i, var_idx) = q(idx_R);
                    else
                        % Determine interface location based on upwind/downwind values
                        theta_R = (2 * q(idx_R) - q_min_R - q_max_R) / (q_max_R - q_min_R);
                        
                        % Convert theta to interface location parameter
                        B_R = theta_R * beta;
                        
                        % THINC reconstruction at the left face (x_{i-1/2})
                        % q(x) = q_min + (q_max - q_min) * (0.5 + 0.5*tanh(beta*(x-x0)/dx))
                        % At x_{i-1/2}, x-x0 = -0.5*dx
                        tanh_term = tanh(B_R - 0.5 * beta);
                        wR_interface(i, var_idx) = q_min_R + (q_max_R - q_min_R) * (0.5 + 0.5 * tanh_term);
                    end
                end
            end
        end
    end
    
    % Ensure non-negative water depths (H)
    wL_interface(:,1) = max(wL_interface(:,1), 0);
    wR_interface(:,1) = max(wR_interface(:,1), 0);
    
    % Handle dry states - set momentum to zero for dry cells
    dry_states_L = wL_interface(:,1) < dry_tol;
    dry_states_R = wR_interface(:,1) < dry_tol;
    
    wL_interface(dry_states_L, 2) = 0;
    wR_interface(dry_states_R, 2) = 0;

end