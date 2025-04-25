function [wL_interface, wR_interface] = weno5(w_padded, cfg)
% WENO5 - Weighted Essentially Non-Oscillatory 5th order reconstruction
%
% Purpose:
%   Performs 5th-order WENO reconstruction, either component-wise (default)
%   or characteristic-wise (if cfg.reconstruct.characteristic = true).
%
% Syntax:
%   [wL_interface, wR_interface] = weno5(w_padded, cfg)
%
% Inputs:
%   w_padded - [(N+2*ng) x 2, double] Padded state array [H, HU].
%   cfg      - [struct] Configuration structure containing:
%              cfg.mesh.N: Number of cells
%              cfg.bc.num_ghost_cells: Number of ghost cells (must be >= 3)
%              cfg.phys.dry_tolerance: Threshold for dry cells
%              cfg.phys.g: Gravity
%              cfg.reconstruct.characteristic: (Optional) Boolean, true for characteristic-wise.
%
% Outputs:
%   wL_interface - [(N+1) x 2, double] Left states [H, HU] at interfaces.
%   wR_interface - [(N+1) x 2, double] Right states [H, HU] at interfaces.
%
% References:
%   - Jiang, G.-S., & Shu, C.-W. (1996).
%   - Shu, C.-W. (1998).
%   - Based on component-wise version by Dr. Denys Dutykh.
%   - Characteristic decomposition adapted from standard methods (e.g., LeVeque 2002).
%
% Author: Dr. Denys Dutykh (Original component-wise)
%         Cascade (Characteristic modification)
% Date: April 25, 2025

    % --- Extract parameters ---
    ng = cfg.bc.num_ghost_cells;  % Number of ghost cells
    N = cfg.mesh.N;               % Number of physical cells
    dry_tol = cfg.phys.dry_tolerance; % Threshold for dry cells
    g = cfg.phys.g;               % Gravity

    % --- Configuration --- 
    % Default to component-wise if not specified
    use_characteristic = isfield(cfg.reconstruct, 'characteristic') && cfg.reconstruct.characteristic;

    % Ensure we have enough ghost cells
    if ng < 3
        error('WENO5 reconstruction requires at least 3 ghost cells. Current: %d', ng);
    end

    % --- Initialize --- 
    wL_interface = zeros(N+1, 2);
    wR_interface = zeros(N+1, 2);
    num_vars = size(w_padded, 2); % Should be 2 (H, HU)

    % --- Loop over interfaces --- 
    for i = 1:N+1
        idx = i + ng - 1; % Center index in padded array for cell i (left of interface i+1/2)

        % --- Boundary Handling Check --- 
        % Check if indices for the full 5-point stencils are valid
        % Stencil for wL_{i+1/2}: {idx-2, idx-1, idx, idx+1, idx+2}
        use_full_wL = (idx >= 3 && idx <= size(w_padded, 1) - 2);
        % Stencil for wR_{i+1/2}: {idx-1, idx, idx+1, idx+2, idx+3}
        use_full_wR = (idx >= 2 && idx <= size(w_padded, 1) - 3);

        if use_characteristic && use_full_wL && use_full_wR
            % --- Characteristic-wise Reconstruction --- 
            
            % 1. Roe Average at interface i+1/2
            wL_local = w_padded(idx, :);   % State in cell i
            wR_local = w_padded(idx+1, :); % State in cell i+1
            [H_roe, U_roe, C_roe] = utils.roe_average(wL_local, wR_local, cfg);
            
            % Avoid issues if Roe state is dry/problematic
            if C_roe < dry_tol
                % Fallback to first-order component-wise if Roe avg is bad
                 wL_interface(i, :) = w_padded(idx, :);   
                 wR_interface(i, :) = w_padded(idx+1, :); 
                 continue; % Skip to next interface
            end

            % 2. Eigenvectors at Roe average state
            [L, R, ~] = utils.sw_eigenvectors(H_roe, U_roe, C_roe, cfg);

            % --- Reconstruct Left State (wL) at i+1/2 --- 
            stencil_indices_L = idx-2 : idx+2;
            W_stencil_L = w_padded(stencil_indices_L, :); % [5 x 2] conservative vars
            W_char_stencil_L = (L * W_stencil_L')'; % [5 x 2] characteristic vars

            reconstructed_char_L = zeros(1, num_vars);
            for k = 1:num_vars % Loop over characteristic fields
                q_char_k = W_char_stencil_L(:, k); % [5 x 1] stencil for k-th char var
                reconstructed_char_L(k) = local_weno5(q_char_k, 'L', cfg); % Pass cfg for epsilon?
            end
            wL_interface(i, :) = (R * reconstructed_char_L')'; % Project back

            % --- Reconstruct Right State (wR) at i+1/2 --- 
            stencil_indices_R = idx-1 : idx+3;
            W_stencil_R = w_padded(stencil_indices_R, :); % [5 x 2] conservative vars
            W_char_stencil_R = (L * W_stencil_R')'; % [5 x 2] characteristic vars

            reconstructed_char_R = zeros(1, num_vars);
            for k = 1:num_vars % Loop over characteristic fields
                q_char_k = W_char_stencil_R(:, k); % [5 x 1] stencil for k-th char var
                reconstructed_char_R(k) = local_weno5(q_char_k, 'R', cfg); % Pass cfg for epsilon?
            end
            wR_interface(i, :) = (R * reconstructed_char_R')'; % Project back

        else
            % --- Component-wise Reconstruction (or 1st Order Boundary Fallback) --- 
            for var_idx = 1:num_vars
                q = w_padded(:, var_idx);
                
                % Left State (wL)
                if use_full_wL
                    q_stencil_L = q(idx-2 : idx+2); % Stencil {i-2, ..., i+2}
                    wL_interface(i, var_idx) = local_weno5(q_stencil_L, 'L', cfg);
                else
                    wL_interface(i, var_idx) = q(idx); % Fallback to 1st order
                end
                
                % Right State (wR)
                if use_full_wR
                    q_stencil_R = q(idx-1 : idx+3); % Stencil {i-1, ..., i+3}
                    wR_interface(i, var_idx) = local_weno5(q_stencil_R, 'R', cfg);
                else
                    wR_interface(i, var_idx) = q(idx+1); % Fallback to 1st order
                end
            end
        end % End characteristic vs component-wise logic

    end % End loop over interfaces

    % --- Ensure non-negative water depths (H) --- 
    wL_interface(:,1) = max(wL_interface(:,1), 0);
    wR_interface(:,1) = max(wR_interface(:,1), 0);

    % --- Handle dry states - set momentum to zero --- 
    dry_states_L = wL_interface(:,1) < dry_tol;
    dry_states_R = wR_interface(:,1) < dry_tol;

    wL_interface(dry_states_L, 2) = 0;
    wR_interface(dry_states_R, 2) = 0;

end

% =========================================================================
% Local Helper Function for Core WENO5 Calculation
% =========================================================================
function q_rec = local_weno5(q_stencil, side, cfg)
    % Performs the core WENO5 calculation on a 5-point stencil.
    % Inputs:
    %   q_stencil - [5x1] vector of values from the stencil.
    %   side      - ['L'|'R'] Indicates whether to compute left or right state.
    %   cfg       - Configuration struct (used for epsilon).
    % Output:
    %   q_rec     - Reconstructed value.

    % WENO parameters
    if isfield(cfg.reconstruct, 'epsilon')
        epsilon = cfg.reconstruct.epsilon;
    else
        epsilon = 1e-6;  % Default small number
    end
    
    % Linear weights (optimal for smooth solutions) - use same for L/R
    d0 = 0.1;  % Weight for stencil shifted left-most relative to interface
    d1 = 0.6;  % Weight for centered stencil
    d2 = 0.3;  % Weight for stencil shifted right-most relative to interface
    
    if side == 'L'
        % For q_{i+1/2}^- using stencil {i-2, i-1, i, i+1, i+2}
        % q_stencil = [q(i-2), q(i-1), q(i), q(i+1), q(i+2)]'
        qim2 = q_stencil(1); qim1 = q_stencil(2); qi = q_stencil(3); qip1 = q_stencil(4); qip2 = q_stencil(5);
        
        % Candidate polynomials (reconstruction values from each substencil)
        v0 = (1/3)*qim2 - (7/6)*qim1 + (11/6)*qi;   % Stencil {i-2, i-1, i}
        v1 = -(1/6)*qim1 + (5/6)*qi   + (1/3)*qip1;  % Stencil {i-1, i,   i+1}
        v2 = (1/3)*qi   + (5/6)*qip1 - (1/6)*qip2;  % Stencil {i,   i+1, i+2}
        
        % Smoothness indicators (IS)
        beta0 = (13/12)*(qim2 - 2*qim1 + qi)^2 + (1/4)*(qim2 - 4*qim1 + 3*qi)^2;
        beta1 = (13/12)*(qim1 - 2*qi + qip1)^2 + (1/4)*(qim1 - qip1)^2;
        beta2 = (13/12)*(qi - 2*qip1 + qip2)^2 + (1/4)*(3*qi - 4*qip1 + qip2)^2;
        
        % Non-linear weights (alpha)
        alpha0 = d0 / ((epsilon + beta0)^2);
        alpha1 = d1 / ((epsilon + beta1)^2);
        alpha2 = d2 / ((epsilon + beta2)^2);
        
    elseif side == 'R'
        % For q_{i+1/2}^+ using stencil {i-1, i, i+1, i+2, i+3}
        % q_stencil = [q(i-1), q(i), q(i+1), q(i+2), q(i+3)]'
        qim1 = q_stencil(1); qi = q_stencil(2); qip1 = q_stencil(3); qip2 = q_stencil(4); qip3 = q_stencil(5);
        
        % Candidate polynomials (reconstruction values from each substencil)
        v0 = (1/3)*qip3 - (7/6)*qip2 + (11/6)*qip1; % Stencil {i+1, i+2, i+3}
        v1 = -(1/6)*qip2 + (5/6)*qip1 + (1/3)*qi;   % Stencil {i,   i+1, i+2}
        v2 = (1/3)*qip1 + (5/6)*qi   - (1/6)*qim1; % Stencil {i-1, i,   i+1}
        
        % Smoothness indicators (IS)
        beta0 = (13/12)*(qip3 - 2*qip2 + qip1)^2 + (1/4)*(qip3 - 4*qip2 + 3*qip1)^2;
        beta1 = (13/12)*(qip2 - 2*qip1 + qi)^2   + (1/4)*(qip2 - qi)^2;
        beta2 = (13/12)*(qip1 - 2*qi + qim1)^2   + (1/4)*(3*qip1 - 4*qi + qim1)^2;
        
        % Non-linear weights (alpha) - Use same linear weights d0, d1, d2
        alpha0 = d0 / ((epsilon + beta0)^2);
        alpha1 = d1 / ((epsilon + beta1)^2);
        alpha2 = d2 / ((epsilon + beta2)^2);
        
    else
        error('Invalid side specified for local_weno5. Use ''L'' or ''R''.');
    end
    
    % Normalize weights (omega)
    alpha_sum = alpha0 + alpha1 + alpha2;
    
    % Avoid division by zero if alpha_sum is extremely small (e.g., all betas large)
    if alpha_sum < eps
        % Fallback to linear weights (or equal weights)
        omega0 = d0; omega1 = d1; omega2 = d2;
    else
        omega0 = alpha0 / alpha_sum;
        omega1 = alpha1 / alpha_sum;
        omega2 = alpha2 / alpha_sum;
    end
    
    % Compute the reconstructed value
    q_rec = omega0*v0 + omega1*v1 + omega2*v2;

end