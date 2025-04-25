function [wL_interface, wR_interface] = ppm(w_padded, cfg)
% PPM - Piecewise Parabolic Method reconstruction
%
% Purpose:
%   Performs 3rd-order accurate PPM (Piecewise Parabolic Method) reconstruction
%   for the finite volume method. PPM uses parabolic profiles within each cell
%   to achieve higher-order accuracy while maintaining non-oscillatory behavior.
%   This implementation supports both component-wise and characteristic-based
%   reconstruction modes.
%
% Syntax:
%   [wL_interface, wR_interface] = ppm(w_padded, cfg)
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
%              cfg.reconstruct.ppm_mode: (optional) 'component' or 'characteristic'
%
% Outputs:
%   wL_interface - [(N+1) x 2, double] Left states at all interfaces
%   wR_interface - [(N+1) x 2, double] Right states at all interfaces
%
% References:
%   - Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM)
%     for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201.
%   - Carpenter, R. L., et al. (1990). Application of the Piecewise Parabolic Method
%     (PPM) to meteorological modeling. Monthly Weather Review, 118(3), 586-612.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: April 24, 2025

    % Extract parameters
    ng = cfg.bc.num_ghost_cells;  % Number of ghost cells
    N = cfg.mesh.N;               % Number of physical cells
    g = cfg.phys.g;               % Gravitational acceleration
    dry_tol = cfg.phys.dry_tolerance;  % Threshold for dry cells
    
    % Check if we have enough ghost cells for PPM (need at least 2)
    if ng < 2
        error('PPM reconstruction requires at least 2 ghost cells. Current: %d', ng);
    end
    
    % Determine reconstruction mode (characteristic-based is default)
    ppm_mode = 'characteristic';
    if isfield(cfg, 'reconstruct') && isfield(cfg.reconstruct, 'ppm_mode')
        ppm_mode = cfg.reconstruct.ppm_mode;
    end
    
    % Initialize arrays for interface values
    wL_interface = zeros(N+1, 2);
    wR_interface = zeros(N+1, 2);
    
    % Process according to the selected mode
    if strcmpi(ppm_mode, 'component')
        % Apply PPM directly to the conservative variables
        [wL_interface, wR_interface] = ppm_component_wise(w_padded, ng, N);
    else
        % Apply PPM to characteristic variables
        [wL_interface, wR_interface] = ppm_characteristic(w_padded, ng, N, g);
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

%% Component-wise PPM reconstruction implementation
function [wL_interface, wR_interface] = ppm_component_wise(w_padded, ng, N)
    % Applies PPM reconstruction directly to conservative variables
    
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
            if idx_L < 2 || idx_R > length(q)-1
                % Near boundary - use first-order
                wL_interface(i, var_idx) = q(idx_L);
                wR_interface(i, var_idx) = q(idx_R);
                continue;
            end
            
            % --- Left cell reconstruction (for right face) ---
            % Compute limited interface values for cell boundaries
            % using parabolic interpolation with monotonicity constraints
            
            % Initial interface value estimates
            qL_minus = 0.5*(q(idx_L-1) + q(idx_L)) - (1/6)*(q(idx_L) - q(idx_L-1));  % q_{i-1/2}
            qL_plus = 0.5*(q(idx_L) + q(idx_L+1)) - (1/6)*(q(idx_L+1) - q(idx_L));   % q_{i+1/2}
            
            % Apply monotonicity constraint - no new extrema
            if (qL_plus - q(idx_L))*(q(idx_L) - qL_minus) <= 0
                qL_minus = q(idx_L);
                qL_plus = q(idx_L);
            elseif 6*(qL_plus - qL_minus)*(q(idx_L) - 0.5*(qL_plus + qL_minus)) > (qL_plus - qL_minus)^2
                qL_minus = 3*q(idx_L) - 2*qL_plus;
            elseif 6*(qL_plus - qL_minus)*(q(idx_L) - 0.5*(qL_plus + qL_minus)) < -(qL_plus - qL_minus)^2
                qL_plus = 3*q(idx_L) - 2*qL_minus;
            end
            
            % Compute parabola coefficients
            qL_bar = q(idx_L);  % Cell average
            dqL = qL_plus - qL_minus;  % Difference across cell
            q6L = 6*(qL_bar - 0.5*(qL_plus + qL_minus));  % Parabolic term
            
            % Evaluate parabola at right face of left cell (x_{i+1/2})
            wL_interface(i, var_idx) = qL_plus;
            
            % --- Right cell reconstruction (for left face) ---
            % Similar process for the right cell
            
            % Initial interface value estimates
            qR_minus = 0.5*(q(idx_R-1) + q(idx_R)) - (1/6)*(q(idx_R) - q(idx_R-1));  % q_{i+1/2}
            qR_plus = 0.5*(q(idx_R) + q(idx_R+1)) - (1/6)*(q(idx_R+1) - q(idx_R));   % q_{i+3/2}
            
            % Apply monotonicity constraint
            if (qR_plus - q(idx_R))*(q(idx_R) - qR_minus) <= 0
                qR_minus = q(idx_R);
                qR_plus = q(idx_R);
            elseif 6*(qR_plus - qR_minus)*(q(idx_R) - 0.5*(qR_plus + qR_minus)) > (qR_plus - qR_minus)^2
                qR_minus = 3*q(idx_R) - 2*qR_plus;
            elseif 6*(qR_plus - qR_minus)*(q(idx_R) - 0.5*(qR_plus + qR_minus)) < -(qR_plus - qR_minus)^2
                qR_plus = 3*q(idx_R) - 2*qR_minus;
            end
            
            % Evaluate parabola at left face of right cell (x_{i+1/2})
            wR_interface(i, var_idx) = qR_minus;
        end
    end
end

%% Characteristic-based PPM reconstruction implementation
function [wL_interface, wR_interface] = ppm_characteristic(w_padded, ng, N, g)
    % Applies PPM reconstruction in characteristic variables
    
    % Initialize arrays for interface values
    wL_interface = zeros(N+1, 2);
    wR_interface = zeros(N+1, 2);
    
    % Loop over interfaces
    for i = 1:N+1
        % Cell indices adjacent to interface i+1/2
        idx_L = i + ng - 1;  % Left cell index (i)
        idx_R = idx_L + 1;   % Right cell index (i+1)
        
        % Check if we have enough cells for a proper stencil
        if idx_L < 2 || idx_R > size(w_padded, 1)-1
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
        
        % Compute inverse of R (left eigenvectors for shallow water equations)
        eps_val = eps; % Use MATLAB's built-in epsilon for stability
        inv_2C = 1.0 / (2.0 * c_roe + eps_val); % Add small epsilon for numerical stability
        R_inv = inv_2C * [ U_roe + c_roe,  -1; 
                          -U_roe + c_roe,   1 ];
        
        % Create stencil for characteristic decomposition
        stencil_L = [idx_L-1, idx_L, idx_L+1];
        stencil_R = [idx_R-1, idx_R, idx_R+1];
        
        % Extract full stencil for left cell
        stencil_data_L = w_padded(stencil_L, :);
        
        % Convert to primitive variables [H, U]
        prim_stencil_L = stencil_data_L;
        for j = 1:3
            if stencil_data_L(j, 1) > 1e-6
                prim_stencil_L(j, 2) = stencil_data_L(j, 2) / stencil_data_L(j, 1);
            else
                prim_stencil_L(j, 2) = 0;
            end
        end
        
        % Project onto characteristic variables
        char_stencil_L = zeros(3, 2);
        for j = 1:3
            char_stencil_L(j, :) = (R_inv * prim_stencil_L(j, :)')';
        end
        
        % Apply PPM to each characteristic variable
        char_wL = zeros(1, 2);
        for var_idx = 1:2
            q = char_stencil_L(:, var_idx);
            
            % Initial interface value estimates
            q_minus = 0.5*(q(1) + q(2)) - (1/6)*(q(2) - q(1));  % q_{i-1/2}
            q_plus = 0.5*(q(2) + q(3)) - (1/6)*(q(3) - q(2));   % q_{i+1/2}
            
            % Apply monotonicity constraint
            if (q_plus - q(2))*(q(2) - q_minus) <= 0
                q_minus = q(2);
                q_plus = q(2);
            elseif 6*(q_plus - q_minus)*(q(2) - 0.5*(q_plus + q_minus)) > (q_plus - q_minus)^2
                q_minus = 3*q(2) - 2*q_plus;
            elseif 6*(q_plus - q_minus)*(q(2) - 0.5*(q_plus + q_minus)) < -(q_plus - q_minus)^2
                q_plus = 3*q(2) - 2*q_minus;
            end
            
            % Evaluate characteristic at right interface
            char_wL(var_idx) = q_plus;
        end
        
        % Project back to primitive variables
        prim_wL = (R * char_wL')';
        
        % Convert to conservative variables
        cons_wL = [prim_wL(1), prim_wL(1)*prim_wL(2)];
        wL_interface(i, :) = cons_wL;
        
        % --- Right Cell (Similar Process) ---
        % Extract stencil for right cell
        stencil_data_R = w_padded(stencil_R, :);
        
        % Convert to primitive variables
        prim_stencil_R = stencil_data_R;
        for j = 1:3
            if stencil_data_R(j, 1) > 1e-6
                prim_stencil_R(j, 2) = stencil_data_R(j, 2) / stencil_data_R(j, 1);
            else
                prim_stencil_R(j, 2) = 0;
            end
        end
        
        % Project onto characteristic variables
        char_stencil_R = zeros(3, 2);
        for j = 1:3
            char_stencil_R(j, :) = (R_inv * prim_stencil_R(j, :)')';
        end
        
        % Apply PPM to each characteristic variable
        char_wR = zeros(1, 2);
        for var_idx = 1:2
            q = char_stencil_R(:, var_idx);
            
            % Initial interface value estimates
            q_minus = 0.5*(q(1) + q(2)) - (1/6)*(q(2) - q(1));  % q_{i+1/2}
            q_plus = 0.5*(q(2) + q(3)) - (1/6)*(q(3) - q(2));   % q_{i+3/2}
            
            % Apply monotonicity constraint
            if (q_plus - q(2))*(q(2) - q_minus) <= 0
                q_minus = q(2);
                q_plus = q(2);
            elseif 6*(q_plus - q_minus)*(q(2) - 0.5*(q_plus + q_minus)) > (q_plus - q_minus)^2
                q_minus = 3*q(2) - 2*q_plus;
            elseif 6*(q_plus - q_minus)*(q(2) - 0.5*(q_plus + q_minus)) < -(q_plus - q_minus)^2
                q_plus = 3*q(2) - 2*q_minus;
            end
            
            % Evaluate characteristic at left interface
            char_wR(var_idx) = q_minus;
        end
        
        % Project back to primitive variables
        prim_wR = (R * char_wR')';
        
        % Convert to conservative variables
        cons_wR = [prim_wR(1), prim_wR(1)*prim_wR(2)];
        wR_interface(i, :) = cons_wR;
    end
end
