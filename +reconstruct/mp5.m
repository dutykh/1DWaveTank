function [wL_interface, wR_interface] = mp5(w_padded, cfg)

% MP5 - Monotonicity Preserving 5th order reconstruction
%
% Purpose:
%   Performs 5th-order MP5 (Monotonicity Preserving) reconstruction for the
%   finite volume method. MP5 achieves high-order accuracy in smooth regions
%   while preserving monotonicity near discontinuities using specialized limiting.
%   This implementation supports both component-wise and characteristic-based
%   reconstruction modes.
%
% Syntax:
%   [wL_interface, wR_interface] = mp5(w_padded, cfg)
%
% Inputs:
%   w_padded - [(N+2*ng) x 2, double] Padded state array including ghost cells.
%             w_padded(:,1) = water depth [m],
%             w_padded(:,2) = discharge [m^2/s]
%   cfg      - [struct] Configuration structure containing:
%              cfg.mesh.N: Number of cells
%              cfg.bc.num_ghost_cells: Number of ghost cells
%              cfg.phys.g: Gravitational acceleration
%              cfg.phys.dry_tolerance: Threshold for dry cells
%              cfg.reconstruct.mp5_mode: (optional) 'component' or 'characteristic'
%
% Outputs:
%   wL_interface - [(N+1) x 2, double] Left states at all interfaces
%   wR_interface - [(N+1) x 2, double] Right states at all interfaces
%
% References:
%   - Suresh, A., & Huynh, H. T. (1997). Accurate monotonicity-preserving
%     schemes with Runge–Kutta time stepping. Journal of Computational
%     Physics, 136(1), 83-99.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: April 24, 2025

    % Extract parameters
    ng = cfg.bc.num_ghost_cells;  % Number of ghost cells
    N = cfg.mesh.N;               % Number of physical cells
    g = cfg.phys.g;               % Gravitational acceleration
    dry_tol = cfg.phys.dry_tolerance;  % Threshold for dry cells
    
    % Check if we have enough ghost cells for MP5 (need at least 3)
    if ng < 3
        error('MP5 reconstruction requires at least 3 ghost cells. Current: %d', ng);
    end
    
    % Determine reconstruction mode (characteristic-based is default)
    mp5_mode = 'characteristic';
    if isfield(cfg, 'reconstruct') && isfield(cfg.reconstruct, 'mp5_mode')
        mp5_mode = cfg.reconstruct.mp5_mode;
    end
    
    % Initialize arrays for interface values
    wL_interface = zeros(N+1, 2);
    wR_interface = zeros(N+1, 2);
    
    % MP5 parameters
    alpha_mp = 4.0;  % MP parameter controlling monotonicity
    
    % Process according to the selected mode
    if strcmpi(mp5_mode, 'component')
        % Apply MP5 directly to the conservative variables
        [wL_interface, wR_interface] = mp5_component_wise(w_padded, ng, N, alpha_mp);
    else
        % Apply MP5 to characteristic variables
        [wL_interface, wR_interface] = mp5_characteristic(w_padded, ng, N, g, alpha_mp, cfg);
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

%% Component-wise MP5 reconstruction implementation
function [wL_interface, wR_interface] = mp5_component_wise(w_padded, ng, N, alpha_mp)
    % Applies MP5 reconstruction directly to conservative variables
    
    % Initialize arrays for interface values
    wL_interface = zeros(N+1, 2);
    wR_interface = zeros(N+1, 2);
    
    % Loop over state variables (H and HU)
    for var_idx = 1:2
        % Extract the component
        q = w_padded(:, var_idx);
        
        % Loop over interfaces
        for i = 1:N+1
            % Cell indices adjacent to interface i+1/2
            idx_L = i + ng - 1;  % Left cell index (i)
            idx_R = idx_L + 1;   % Right cell index (i+1)
            
            % Check if we have enough cells for a proper stencil
            if idx_L < 3 || idx_R > length(q)-2
                % Near boundary - use first-order
                wL_interface(i, var_idx) = q(idx_L);
                wR_interface(i, var_idx) = q(idx_R);
                continue;
            end
            
            % --- Left cell reconstruction (for right face) ---
            % For 5th order MP method, we need a 5-point stencil centered at i
            % {i-2, i-1, i, i+1, i+2}
            
            % Step 1: Compute the unconstrained 5th order polynomial interpolation
            q_im2 = q(idx_L-2);
            q_im1 = q(idx_L-1);
            q_i   = q(idx_L);
            q_ip1 = q(idx_L+1);
            q_ip2 = q(idx_L+2);
            
            % 5th order polynomial interpolation at x_{i+1/2}
            q_L = (2*q_im2 - 13*q_im1 + 47*q_i + 27*q_ip1 - 3*q_ip2)/60;
            
            % Step 2: Monotonicity-preserving constraint
            % Compute bounds for monotonicity preservation
            dq_min = min([q_i, q_ip1]) - min([q_im1, q_i, q_ip1]);
            dq_max = max([q_i, q_ip1]) - max([q_im1, q_i, q_ip1]);
            
            % Add buffer parameter alpha
            dq_min = q_i + alpha_mp * dq_min;
            dq_max = q_i + alpha_mp * dq_max;
            
            % Apply bounds to ensure monotonicity
            q_MP = min(max(q_L, dq_min), dq_max);
            
            % Step 3: Additional constraints to ensure MP property
            % Compute median-based constraints
            q_UL = q_i + alpha_mp * (q_i - q_im1);
            q_AV = 0.5 * (q_i + q_ip1);
            q_MD = q_AV - 0.5 * (q_ip1 - q_i);
            q_LC = q_i + 0.5 * (q_i - q_im1);
            
            % Monotonicity constraint from Suresh & Huynh
            if (q_i - q_im1)*(q_ip1 - q_i) <= 0
                % Local extremum - use center value
                q_final = q_i;
            else
                % Apply MP constraints
                % First compute the MP bounds
                q_min = min([q_i, q_ip1, q_MD]);
                q_max = max([q_i, q_ip1, q_MD]);
                
                % Then apply the bounds to the preliminary value
                q_MP = min(max(q_MP, q_min), q_max);
                
                % Final condition based on the local curvature
                if (q_MP - q_i)*(q_MP - q_UL) <= 0
                    q_final = q_MP;
                else if (q_i - q_im1)*(q_i - q_im2) <= 0
                        q_final = q_LC;
                    else
                        q_final = q_UL;
                    end
                end
            end
            
            % Store the final interpolated value
            wL_interface(i, var_idx) = q_final;
            
            % --- Right cell reconstruction (for left face) ---
            % Similar process but mirrored for the right cell
            
            % For right cell, we need stencil {i-1, i, i+1, i+2, i+3}
            q_im1 = q(idx_R-1);
            q_i   = q(idx_R);
            q_ip1 = q(idx_R+1);
            q_ip2 = q(idx_R+2);
            
            % For the left face of right cell, we use the mirrored formula
            % The stencil shifts by one compared to the left cell
            q_R = (2*q_ip2 - 13*q_ip1 + 47*q_i + 27*q_im1 - 3*q(idx_R-2))/60;
            
            % Apply monotonicity constraints (similar to above but mirrored)
            dq_min = min([q_i, q_im1]) - min([q_ip1, q_i, q_im1]);
            dq_max = max([q_i, q_im1]) - max([q_ip1, q_i, q_im1]);
            
            dq_min = q_i + alpha_mp * dq_min;
            dq_max = q_i + alpha_mp * dq_max;
            
            q_MP = min(max(q_R, dq_min), dq_max);
            
            q_UL = q_i + alpha_mp * (q_i - q_ip1);
            q_AV = 0.5 * (q_i + q_im1);
            q_MD = q_AV - 0.5 * (q_im1 - q_i);
            q_LC = q_i + 0.5 * (q_i - q_ip1);
            
            if (q_i - q_ip1)*(q_im1 - q_i) <= 0
                q_final = q_i;
            else
                q_min = min([q_i, q_im1, q_MD]);
                q_max = max([q_i, q_im1, q_MD]);
                
                q_MP = min(max(q_MP, q_min), q_max);
                
                if (q_MP - q_i)*(q_MP - q_UL) <= 0
                    q_final = q_MP;
                else if (q_i - q_ip1)*(q_i - q_ip2) <= 0
                        q_final = q_LC;
                    else
                        q_final = q_UL;
                    end
                end
            end
            
            wR_interface(i, var_idx) = q_final;
        end
    end
end

%% Characteristic-based MP5 reconstruction implementation
function [wL_interface, wR_interface] = mp5_characteristic(w_padded, ng, N, g, alpha_mp, cfg)
    % Applies MP5 reconstruction in characteristic variables
    
    % Initialize arrays for interface values
    wL_interface = zeros(N+1, 2);
    wR_interface = zeros(N+1, 2);
    
    % Loop over interfaces
    for i = 1:N+1
        % Cell indices adjacent to interface i+1/2
        idx_L = i + ng - 1;  % Left cell index (i)
        idx_R = idx_L + 1;   % Right cell index (i+1)
        
        % Check if we have enough cells for a proper stencil
        if idx_L < 3 || idx_R > size(w_padded, 1)-2
            % Near boundary - use first-order
            wL_interface(i, :) = w_padded(idx_L, :);
            wR_interface(i, :) = w_padded(idx_R, :);
            continue;
        end
        
        % Use Roe average state for the transformation matrix
        HL = w_padded(idx_L, 1);
        HR = w_padded(idx_R, 1);
        if HL < 1e-6 || HR < 1e-6
            % One of the cells is dry - revert to component-wise
            wL_interface(i, :) = w_padded(idx_L, :);
            wR_interface(i, :) = w_padded(idx_R, :);
            continue;
        end
        
        % Calculate Roe-averaged state at the interface
        HL_sqrt = sqrt(HL);
        HR_sqrt = sqrt(HR);
        H_roe = 0.5 * (HL + HR); % Average depth
        
        % Compute Roe-averaged velocity
        UL = w_padded(idx_L, 2) / HL;
        UR = w_padded(idx_R, 2) / HR;
        U_roe = (HL_sqrt * UL + HR_sqrt * UR) / (HL_sqrt + HR_sqrt);
        
        % Wave speed at Roe state
        c_roe = sqrt(g * H_roe);
        
        % --- Transform to Characteristic Variables ---
        % Characteristic variables:
        % W1 = U - 2*sqrt(g*H)  (left-going)
        % W2 = U + 2*sqrt(g*H)  (right-going)
        
        % Create transformation matrix and its inverse
        % For shallow water, we use the eigenvectors of the Jacobian
        R = [1, 1; U_roe-c_roe, U_roe+c_roe];  % Right eigenvectors as columns
        
        % --- Corrected R_inv Calculation ---
        num_eps = cfg.numerics.epsilon; % Use config epsilon for stability
        inv_2C = 1.0 / (2.0 * c_roe + num_eps); % Add epsilon to denominator
        R_inv = inv_2C * [ U_roe+c_roe,  -1; 
                          -U_roe+c_roe,   1 ];
        % -----------------------------------
        
        % --- Left Cell MP5 Reconstruction ---
        % We need 5 cells for MP5: {i-2, i-1, i, i+1, i+2}
        stencil_L = [idx_L-2, idx_L-1, idx_L, idx_L+1, idx_L+2];
        
        % Extract and convert to primitive variables
        prim_stencil_L = zeros(5, 2);
        for j = 1:5
            prim_stencil_L(j, 1) = w_padded(stencil_L(j), 1);  % H
            if prim_stencil_L(j, 1) > cfg.phys.dry_tolerance
                prim_stencil_L(j, 2) = w_padded(stencil_L(j), 2) / prim_stencil_L(j, 1); % U
            else
                prim_stencil_L(j, 2) = 0;
            end
        end
        
        % Project onto characteristic variables
        char_stencil_L = zeros(5, 2);
        for j = 1:5
            char_stencil_L(j, :) = (R_inv * prim_stencil_L(j, :)')';
        end
        
        % Apply MP5 reconstruction to each characteristic variable
        char_wL = zeros(1, 2);
        for var_idx = 1:2
            % Extract the characteristic variable from the stencil
            q_im2 = char_stencil_L(1, var_idx);
            q_im1 = char_stencil_L(2, var_idx);
            q_i   = char_stencil_L(3, var_idx);
            q_ip1 = char_stencil_L(4, var_idx);
            q_ip2 = char_stencil_L(5, var_idx);
            
            % 5th order polynomial interpolation at x_{i+1/2}
            q_L = (2*q_im2 - 13*q_im1 + 47*q_i + 27*q_ip1 - 3*q_ip2)/60;
            
            % Apply monotonicity constraints
            dq_min = min([q_i, q_ip1]) - min([q_im1, q_i, q_ip1]);
            dq_max = max([q_i, q_ip1]) - max([q_im1, q_i, q_ip1]);
            
            dq_min = q_i + alpha_mp * dq_min;
            dq_max = q_i + alpha_mp * dq_max;
            
            q_MP = min(max(q_L, dq_min), dq_max);
            
            q_UL = q_i + alpha_mp * (q_i - q_im1);
            q_AV = 0.5 * (q_i + q_ip1);
            q_MD = q_AV - 0.5 * (q_ip1 - q_i);
            q_LC = q_i + 0.5 * (q_i - q_im1);
            
            if (q_i - q_im1)*(q_ip1 - q_i) <= 0
                q_final = q_i;
            else
                q_min = min([q_i, q_ip1, q_MD]);
                q_max = max([q_i, q_ip1, q_MD]);
                
                q_MP = min(max(q_MP, q_min), q_max);
                
                if (q_MP - q_i)*(q_MP - q_UL) <= 0
                    q_final = q_MP;
                else if (q_i - q_im1)*(q_i - q_im2) <= 0
                        q_final = q_LC;
                    else
                        q_final = q_UL;
                    end
                end
            end
            
            char_wL(var_idx) = q_final;
        end
        
        % Project back to primitive variables
        prim_wL = (R * char_wL')';
        
        % Convert to conservative variables
        cons_wL = [prim_wL(1), prim_wL(1)*prim_wL(2)];
        wL_interface(i, :) = cons_wL;
        
        % --- Right Cell MP5 Reconstruction ---
        % We need 5 cells for MP5: {i-1, i, i+1, i+2, i+3}
        stencil_R = [idx_R-2, idx_R-1, idx_R, idx_R+1, idx_R+2];
        
        % Extract and convert to primitive variables
        prim_stencil_R = zeros(5, 2);
        for j = 1:5
            prim_stencil_R(j, 1) = w_padded(stencil_R(j), 1);  % H
            if prim_stencil_R(j, 1) > cfg.phys.dry_tolerance
                prim_stencil_R(j, 2) = w_padded(stencil_R(j), 2) / prim_stencil_R(j, 1); % U
            else
                prim_stencil_R(j, 2) = 0;
            end
        end
        
        % Project onto characteristic variables
        char_stencil_R = zeros(5, 2);
        for j = 1:5
            char_stencil_R(j, :) = (R_inv * prim_stencil_R(j, :)')';
        end
        
        % Apply MP5 reconstruction to each characteristic variable
        char_wR = zeros(1, 2);
        for var_idx = 1:2
            % Extract the characteristic variable from the stencil
            q_im2 = char_stencil_R(1, var_idx);
            q_im1 = char_stencil_R(2, var_idx);
            q_i   = char_stencil_R(3, var_idx);
            q_ip1 = char_stencil_R(4, var_idx);
            q_ip2 = char_stencil_R(5, var_idx);
            
            % For the left face, we use the mirrored formula
            q_R = (2*q_ip2 - 13*q_ip1 + 47*q_i + 27*q_im1 - 3*q_im2)/60;
            
            % Apply mirrored monotonicity constraints 
            dq_min = min([q_i, q_im1]) - min([q_ip1, q_i, q_im1]);
            dq_max = max([q_i, q_im1]) - max([q_ip1, q_i, q_im1]);
            
            dq_min = q_i + alpha_mp * dq_min;
            dq_max = q_i + alpha_mp * dq_max;
            
            q_MP = min(max(q_R, dq_min), dq_max);
            
            q_UL = q_i + alpha_mp * (q_i - q_ip1);
            q_AV = 0.5 * (q_i + q_im1);
            q_MD = q_AV - 0.5 * (q_im1 - q_i);
            q_LC = q_i + 0.5 * (q_i - q_ip1);
            
            if (q_i - q_ip1)*(q_im1 - q_i) <= 0
                q_final = q_i;
            else
                q_min = min([q_i, q_im1, q_MD]);
                q_max = max([q_i, q_im1, q_MD]);
                
                q_MP = min(max(q_MP, q_min), q_max);
                
                if (q_MP - q_i)*(q_MP - q_UL) <= 0
                    q_final = q_MP;
                else if (q_i - q_ip1)*(q_i - q_ip2) <= 0
                        q_final = q_LC;
                    else
                        q_final = q_UL;
                    end
                end
            end
            
            char_wR(var_idx) = q_final;
        end
        
        % Project back to primitive variables
        prim_wR = (R * char_wR')';
        
        % Convert to conservative variables
        cons_wR = [prim_wR(1), prim_wR(1)*prim_wR(2)];
        wR_interface(i, :) = cons_wR;
    end

end