function [wL_interface, wR_interface] = uno2(w_padded, cfg)

    % UNO2 - Uniformly Non-Oscillatory 2nd order reconstruction
%
% Purpose:
%   Performs 2nd-order UNO (Uniformly Non-Oscillatory) reconstruction for
%   the finite volume method. UNO is a variant of ENO that uses a modified
%   limiter approach to achieve higher order accuracy while maintaining
%   non-oscillatory properties.
%
% Syntax:
%   [wL_interface, wR_interface] = uno2(w_padded, cfg)
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
%   - Harten, A., & Osher, S. (1987). Uniformly high-order accurate 
%     nonoscillatory schemes, I. SIAM Journal on Numerical Analysis, 24(2), 279-309.
%   - Marquina, A. (1994). Local piecewise hyperbolic reconstruction of 
%     numerical fluxes for nonlinear scalar conservation laws. 
%     SIAM Journal on Scientific Computing, 15(4), 892-915.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: April 24, 2025

    % Extract parameters
    ng = cfg.bc.num_ghost_cells;  % Number of ghost cells
    N = cfg.mesh.N;              % Number of physical cells
    dry_tol = cfg.phys.dry_tolerance; % Threshold for dry cells

    % Initialize arrays for interface values
    wL_interface = zeros(N+1, 2);
    wR_interface = zeros(N+1, 2);
    
    % Loop over state variables (H and HU)
    for var_idx = 1:2
        % Extract the component
        q = w_padded(:, var_idx);
        
        % UNO2 reconstruction for this variable
        for i = 1:N+1
            % Left and right cells adjacent to the interface
            idx_L = i + ng - 1;
            idx_R = idx_L + 1;
            
            % We need values from a wider stencil for UNO2
            idx_LL = idx_L - 1; % Two cells to the left
            idx_RR = idx_R + 1; % Two cells to the right
            
            % Check if we have valid cell indices (not at boundaries)
            if idx_LL >= 1 && idx_RR <= N+2*ng
                % UNO2 slope limiter calculation
                
                % Forward differences
                dq_L = q(idx_L) - q(idx_LL);  % q_i - q_{i-1}
                dq_C = q(idx_R) - q(idx_L);   % q_{i+1} - q_i
                dq_R = q(idx_RR) - q(idx_R);  % q_{i+2} - q_{i+1}
                
                % UNO2 limiter - left interface
                % The key idea of UNO2 is to modify the minmod limiter to allow more
                % accuracy while still preventing oscillations
                
                % Left slope estimation using UNO2 approach
                % First get the standard minmod of adjacent differences
                if dq_L * dq_C <= 0
                    slope_L = 0; % Opposite signs - set to zero
                else
                    % Same sign - relaxed limiter
                    % UNO calculates the bounded slope that prevents creating new extrema
                    % while allowing for higher accuracy than standard minmod
                    dq_UNO_L = minmod2(2*dq_L - dq_C, 2*dq_C - dq_L, dq_L, dq_C);
                    slope_L = dq_UNO_L / 2; % Divide by 2 for the half-cell extrapolation
                end
                
                % Right slope estimation using UNO2 approach
                if dq_C * dq_R <= 0
                    slope_R = 0; % Opposite signs - set to zero
                else
                    % Same sign - relaxed limiter
                    dq_UNO_R = minmod2(2*dq_C - dq_R, 2*dq_R - dq_C, dq_C, dq_R);
                    slope_R = dq_UNO_R / 2; % Divide by 2 for the half-cell extrapolation
                end
                
                % Extrapolate to interfaces using the limited slopes
                wL_interface(i, var_idx) = q(idx_L) + slope_L;
                wR_interface(i, var_idx) = q(idx_R) - slope_R;
            else
                % Near boundary - revert to first-order (piecewise constant)
                wL_interface(i, var_idx) = q(idx_L);
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

% Helper function: modified minmod limiter for UNO2
function result = minmod2(a, b, c, d)

    % UNO-specific minmod function that takes 4 arguments
    % This allows more accurate reconstruction than standard minmod
    % while still maintaining the non-oscillatory property
    
    if (a*b <= 0) || (a*c <= 0) || (a*d <= 0)
        result = 0;
    else
        % All values have same sign
        result = sign(a) * min([abs(a), abs(b), abs(c), abs(d)]);
    end

end