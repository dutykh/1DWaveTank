function [wL_interface, wR_interface] = eno2(w_padded, cfg)

% ENO2 - Essentially Non-Oscillatory 2nd order reconstruction
%
% Purpose:
%   Performs 2nd-order ENO (Essentially Non-Oscillatory) reconstruction for
%   the finite volume method. ENO schemes adaptively choose the smoothest
%   stencil by comparing undivided differences, maintaining high-order
%   accuracy while avoiding oscillations near discontinuities.
%
% Syntax:
%   [wL_interface, wR_interface] = eno2(w_padded, cfg)
%
% Inputs:
%   w_padded - [(N+2*ng) x 2, double] Padded state array including ghost cells.
%             w_padded(:,1) = water depth [m],
%             w_padded(:,2) = discharge [m^2/s]
%   cfg      - [struct] Configuration structure containing:
%              cfg.mesh.N: Number of cells
%              cfg.bc.num_ghost_cells: Number of ghost cells
%              cfg.phys.dry_tolerance: Threshold for dry cells
%
% Outputs:
%   wL_interface - [(N+1) x 2, double] Left states at all interfaces
%   wR_interface - [(N+1) x 2, double] Right states at all interfaces
%
% References:
%   - Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987).
%     Uniformly high order accurate essentially non-oscillatory schemes, III.
%     Journal of Computational Physics, 71(2), 231-303.
%   - Shu, C.-W. (1998). Essentially non-oscillatory and weighted essentially
%     non-oscillatory schemes for hyperbolic conservation laws.
%     In Advanced Numerical Approximation of Nonlinear Hyperbolic Equations.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: April 24, 2025

    % Extract parameters
    ng = cfg.bc.num_ghost_cells;  % Number of ghost cells
    N = cfg.mesh.N;               % Number of physical cells
    dry_tol = cfg.phys.dry_tolerance; % Threshold for dry cells

    % Initialize arrays for interface values
    wL_interface = zeros(N+1, 2);
    wR_interface = zeros(N+1, 2);
    
    % Loop over state variables (H and HU)
    for var_idx = 1:2
        % Extract the component
        q = w_padded(:, var_idx);
        
        % ENO2 reconstruction for this variable
        for i = 1:N+1
            % Left and right cells adjacent to the interface
            idx_L = i + ng - 1;  % Left cell index (i)
            idx_R = idx_L + 1;   % Right cell index (i+1)
            
            % --- Left Interface State (from left cell) ---
            if idx_L > 1 && idx_L < N+2*ng
                % Calculate undivided differences (slopes)
                diff_L = q(idx_L) - q(idx_L-1);  % q_i - q_{i-1}
                diff_R = q(idx_L+1) - q(idx_L);  % q_{i+1} - q_i
                
                % Choose smoother stencil by comparing absolute differences
                if abs(diff_L) <= abs(diff_R)
                    % Use left-biased stencil {i-1, i}
                    slope = diff_L;
                else
                    % Use right-biased stencil {i, i+1}
                    slope = diff_R;
                end
                
                % Reconstruct to right face of cell idx_L (i.e., x_{i+1/2})
                wL_interface(i, var_idx) = q(idx_L) + 0.5 * slope;
            else
                % Near boundary - use first-order (piecewise constant)
                wL_interface(i, var_idx) = q(idx_L);
            end
            
            % --- Right Interface State (from right cell) ---
            if idx_R > 1 && idx_R < N+2*ng
                % Calculate undivided differences (slopes)
                diff_L = q(idx_R) - q(idx_R-1);  % q_{i+1} - q_i
                diff_R = q(idx_R+1) - q(idx_R);  % q_{i+2} - q_{i+1}
                
                % Choose smoother stencil by comparing absolute differences
                if abs(diff_L) <= abs(diff_R)
                    % Use left-biased stencil {i, i+1}
                    slope = diff_L;
                else
                    % Use right-biased stencil {i+1, i+2}
                    slope = diff_R;
                end
                
                % Reconstruct to left face of cell idx_R (i.e., x_{i+1/2})
                wR_interface(i, var_idx) = q(idx_R) - 0.5 * slope;
            else
                % Near boundary - use first-order (piecewise constant)
                wR_interface(i, var_idx) = q(idx_R);
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