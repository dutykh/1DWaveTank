function [wL_interface, wR_interface] = muscl(w_padded, cfg)
% MUSCL Second-order MUSCL reconstruction for shallow water equations.
%
% Purpose:
%   Performs second-order spatial reconstruction of cell-averaged data to
%   interface values using the MUSCL approach with a slope limiter.
%   Supports both component-wise and characteristic-based reconstruction.
%
% Syntax:
%   [wL_interface, wR_interface] = muscl(w_padded, cfg)
%
% Inputs:
%   w_padded - [(N+2*ng) x 2, double] Padded state array including ghost cells.
%             w_padded(:,1) = water depth H [m],
%             w_padded(:,2) = discharge HU [m^2/s]
%   cfg      - [struct] Configuration structure containing:
%              cfg.mesh.N: Number of cells
%              cfg.bc.num_ghost_cells: Number of ghost cells (ng >= 2 needed)
%              cfg.reconstruct.limiter: Limiter function handle (e.g., @minmod)
%              cfg.reconstruct.characteristic: (Optional) Boolean, true to use
%                                              characteristic reconstruction (default: false).
%              cfg.phys.g: Gravitational acceleration
%              cfg.phys.dry_tolerance: Threshold for dry cells
%
% Outputs:
%   wL_interface - [(N+1) x 2, double] Left conservative states at all interfaces [H; HU]
%   wR_interface - [(N+1) x 2, double] Right conservative states at all interfaces [H; HU]
%
% Description:
%   If cfg.reconstruct.characteristic is true, reconstruction is performed on
%   Riemann invariants (W1 = u - 2c, W2 = u + 2c) to improve stability near shocks.
%   Otherwise, reconstruction is performed component-wise on H and HU.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   25 April 2025

    % Extract parameters
    ng = cfg.bc.num_ghost_cells;
    N = cfg.mesh.N;
    g = cfg.phys.g;
    dry_tolerance = cfg.phys.dry_tolerance;
    total_cells_padded = N + 2*ng;

    % Check for characteristic reconstruction flag
    use_characteristic = isfield(cfg.reconstruct, 'characteristic') && cfg.reconstruct.characteristic;

    % Get limiter function handle
    if isfield(cfg.reconstruct, 'limiter_handle') && ~isempty(cfg.reconstruct.limiter_handle)
        limiter_handle = cfg.reconstruct.limiter_handle;
    else
        limiter_handle = @reconstruct.limiters.minmod;
        warning('MUSCL:NoLimiter', 'No slope limiter specified, using minmod by default.');
    end

    % Initialize output arrays
    wL_interface = zeros(N+1, 2);
    wR_interface = zeros(N+1, 2);

    % Extract H and HU
    H_padded = w_padded(:, 1);
    HU_padded = w_padded(:, 2);

    if use_characteristic
        % --- Characteristic Reconstruction --- 
        
        % Convert to primitive variables (H, U)
        U_padded = zeros(size(H_padded));
        idx_wet = H_padded > dry_tolerance;
        U_padded(idx_wet) = HU_padded(idx_wet) ./ H_padded(idx_wet);
        
        % Compute Riemann invariants (W1=u-2c, W2=u+2c)
        W1 = zeros(size(H_padded));
        W2 = zeros(size(H_padded));
        c_padded = sqrt(g * H_padded(idx_wet)); % Wave speed
        W1(idx_wet) = U_padded(idx_wet) - 2 * c_padded;
        W2(idx_wet) = U_padded(idx_wet) + 2 * c_padded;
        
        % Calculate limited slopes for characteristic variables
        slopes_W1 = zeros(size(W1));
        slopes_W2 = zeros(size(W2));
        for i = 2:(total_cells_padded-1)
            % Only compute slopes where all cells involved are wet
            if H_padded(i-1) > dry_tolerance && H_padded(i) > dry_tolerance && H_padded(i+1) > dry_tolerance
                % W1 (left-going)
                delta_minus_W1 = W1(i) - W1(i-1);
                delta_plus_W1 = W1(i+1) - W1(i);
                slopes_W1(i) = limiter_handle(delta_minus_W1, delta_plus_W1);
                
                % W2 (right-going)
                delta_minus_W2 = W2(i) - W2(i-1);
                delta_plus_W2 = W2(i+1) - W2(i);
                slopes_W2(i) = limiter_handle(delta_minus_W2, delta_plus_W2);
            end % else slopes remain zero
        end
        
        % Reconstruct characteristic variables at interfaces
        for i = 1:N+1
            idx_L = i + ng - 1; % Left cell index
            idx_R = i + ng;     % Right cell index
            
            % If either neighbor cell is dry, fallback to 1st order for safety?
            % Or rely on the slope being zero from the calculation above?
            % Current muscl_characteristic uses 1st order if *either* cell is dry.
            if H_padded(idx_L) > dry_tolerance && H_padded(idx_R) > dry_tolerance
                % Reconstruct characteristics at interface i+1/2
                W1_L = W1(idx_L) + 0.5 * slopes_W1(idx_L);
                W1_R = W1(idx_R) - 0.5 * slopes_W1(idx_R);
                W2_L = W2(idx_L) + 0.5 * slopes_W2(idx_L);
                W2_R = W2(idx_R) - 0.5 * slopes_W2(idx_R);
                
                % Transform back to primitive variables (H, U)
                U_L = 0.5 * (W1_L + W2_L);
                c_L = 0.25 * (W2_L - W1_L);
                H_L = max((c_L * c_L) / g, 0); % Ensure non-negative H
                
                U_R = 0.5 * (W1_R + W2_R);
                c_R = 0.25 * (W2_R - W1_R);
                H_R = max((c_R * c_R) / g, 0); % Ensure non-negative H
                
                % Convert back to conservative variables [H; HU]
                wL_interface(i, 1) = H_L;
                wL_interface(i, 2) = H_L * U_L;
                wR_interface(i, 1) = H_R;
                wR_interface(i, 2) = H_R * U_R;
                
            else % One or both cells near interface are dry: Use 1st order
                wL_interface(i, 1) = H_padded(idx_L);
                wL_interface(i, 2) = HU_padded(idx_L);
                wR_interface(i, 1) = H_padded(idx_R);
                wR_interface(i, 2) = HU_padded(idx_R);
            end
        end
        
    else
        % --- Component-wise Reconstruction --- 
        slopes = zeros(total_cells_padded, 2); % Slopes for H, HU
        
        % Step 1: Calculate Limited Slopes for Conservative Variables (H, HU)
        for var_idx = 1:2 % Loop over H (1) and HU (2)
            q = w_padded(:, var_idx);
            current_var_slopes = zeros(total_cells_padded, 1);
            for i = 2:(total_cells_padded-1)
                delta_minus = q(i) - q(i-1);
                delta_plus = q(i+1) - q(i);
                current_var_slopes(i) = limiter_handle(delta_minus, delta_plus);
            end
            slopes(:, var_idx) = current_var_slopes;
        end
        
        % Step 2: Reconstruct Conservative Variables at Interfaces
        for i = 1:N+1
            idx_left = i + ng - 1;
            idx_right = i + ng;
            for var_idx = 1:2
                q = w_padded(:, var_idx);
                current_slopes = slopes(:, var_idx);
                % Left state at interface i+1/2: q_L + 0.5*slope_L
                wL_interface(i, var_idx) = q(idx_left) + 0.5 * current_slopes(idx_left);
                % Right state at interface i+1/2: q_R - 0.5*slope_R
                wR_interface(i, var_idx) = q(idx_right) - 0.5 * current_slopes(idx_right);
            end
        end
    end % End characteristic vs component-wise

    % --- Final Step: Ensure Positivity and Handle Dry States --- 
    % This applies regardless of reconstruction method
    for i = 1:N+1
        % Check Left State
        if wL_interface(i,1) < dry_tolerance
            wL_interface(i,1) = 0; % Enforce non-negative depth
            wL_interface(i,2) = 0; % Set momentum to zero if dry
        end
        
        % Check Right State
        if wR_interface(i,1) < dry_tolerance
            wR_interface(i,1) = 0; % Enforce non-negative depth
            wR_interface(i,2) = 0; % Set momentum to zero if dry
        end
    end
end
