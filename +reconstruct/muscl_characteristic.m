function [wL_interface, wR_interface] = muscl_characteristic(w_padded, cfg)
% MUSCL_CHARACTERISTIC MUSCL scheme with characteristic decomposition
%
% Purpose:
%   Performs second-order reconstruction in characteristic variables
%   for the shallow water equations, providing better oscillation control
%   at discontinuities than component-wise reconstruction.
%
% Syntax:
%   [wL_interface, wR_interface] = muscl_characteristic(w_padded, cfg)
%
% Inputs:
%   w_padded - [(N+2*ng) x 2, double] Padded state array including ghost cells.
%             w_padded(:,1) = water depth [m],
%             w_padded(:,2) = discharge [m^2/s]
%   cfg      - [struct] Configuration structure containing:
%              cfg.mesh.N: Number of cells
%              cfg.bc.num_ghost_cells: Number of ghost cells
%              cfg.reconstruct.limiter: Limiter function handle
%              cfg.phys.g: Gravitational acceleration
%              cfg.phys.dry_tolerance: Threshold for dry cells
%
% Outputs:
%   wL_interface - [(N+1) x 2, double] Left conservative states at all interfaces
%   wR_interface - [(N+1) x 2, double] Right conservative states at all interfaces
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: April 24, 2025

% Extract parameters
ng = cfg.bc.num_ghost_cells;  % Number of ghost cells
N = cfg.mesh.N;               % Number of physical cells
g = cfg.phys.g;               % Gravitational acceleration
dry_tol = cfg.phys.dry_tolerance; % Threshold for dry cells

% Check if a limiter function is specified
if ~isfield(cfg.reconstruct, 'limiter') || isempty(cfg.reconstruct.limiter)
    limiter_handle = @reconstruct.limiters.minmod; % Default limiter
else
    limiter_handle = cfg.reconstruct.limiter;
end

% Preallocate arrays for interface values
wL_interface = zeros(N+1, 2);
wR_interface = zeros(N+1, 2);

% Extract primitive variables (H, U) from conservative variables (H, HU)
H_padded = w_padded(:, 1);
HU_padded = w_padded(:, 2);
U_padded = zeros(size(H_padded));
idx_wet = H_padded > dry_tol;
U_padded(idx_wet) = HU_padded(idx_wet) ./ H_padded(idx_wet);

% Compute characteristic variables (Riemann invariants)
% W1 = u - 2√(gh) (left-going characteristic)
% W2 = u + 2√(gh) (right-going characteristic)
W1 = zeros(size(H_padded));
W2 = zeros(size(H_padded));
idx_wet = H_padded > dry_tol;
c_padded = sqrt(g * H_padded(idx_wet));  % Wave speed
W1(idx_wet) = U_padded(idx_wet) - 2 * c_padded;
W2(idx_wet) = U_padded(idx_wet) + 2 * c_padded;

% Calculate slopes for characteristic variables
slopes_W1 = zeros(size(W1));
slopes_W2 = zeros(size(W2));

for i = 2:(N+2*ng-1)
    % Only compute slopes where all cells involved are wet
    if H_padded(i-1) > dry_tol && H_padded(i) > dry_tol && H_padded(i+1) > dry_tol
        % Characteristic variable 1 (left-going)
        delta_minus = W1(i) - W1(i-1);
        delta_plus = W1(i+1) - W1(i);
        slopes_W1(i) = limiter_handle(delta_minus, delta_plus);
        
        % Characteristic variable 2 (right-going)
        delta_minus = W2(i) - W2(i-1);
        delta_plus = W2(i+1) - W2(i);
        slopes_W2(i) = limiter_handle(delta_minus, delta_plus);
    end
end

% Reconstruct characteristic variables at interfaces
for i = 1:N+1
    % Index in padded array
    idx_L = i + ng - 1;  % Left cell index
    idx_R = i + ng;      % Right cell index
    
    % Only process if both cells are wet
    if H_padded(idx_L) > dry_tol && H_padded(idx_R) > dry_tol
        % Reconstruct characteristic variables at interface
        W1_L = W1(idx_L) + 0.5 * slopes_W1(idx_L);
        W1_R = W1(idx_R) - 0.5 * slopes_W1(idx_R);
        
        W2_L = W2(idx_L) + 0.5 * slopes_W2(idx_L);
        W2_R = W2(idx_R) - 0.5 * slopes_W2(idx_R);
        
        % Transform back to primitive variables (H, U)
        % U = (W1 + W2)/2, c = (W2 - W1)/4
        U_L = 0.5 * (W1_L + W2_L);
        c_L = 0.25 * (W2_L - W1_L);
        H_L = (c_L * c_L) / g;
        
        U_R = 0.5 * (W1_R + W2_R);
        c_R = 0.25 * (W2_R - W1_R);
        H_R = (c_R * c_R) / g;
        
        % Ensure non-negative water depth
        H_L = max(H_L, 0);
        H_R = max(H_R, 0);
        
        % Convert to conservative variables
        wL_interface(i, 1) = H_L;
        wL_interface(i, 2) = H_L * U_L;
        
        wR_interface(i, 1) = H_R;
        wR_interface(i, 2) = H_R * U_R;
    else
        % One or both cells are dry - use first-order
        wL_interface(i, 1) = H_padded(idx_L);
        wL_interface(i, 2) = HU_padded(idx_L);
        
        wR_interface(i, 1) = H_padded(idx_R);
        wR_interface(i, 2) = HU_padded(idx_R);
    end
    
    % Safety check for dry states in the result
    if wL_interface(i, 1) < dry_tol
        wL_interface(i, 1) = 0;
        wL_interface(i, 2) = 0;
    end
    
    if wR_interface(i, 1) < dry_tol
        wR_interface(i, 1) = 0;
        wR_interface(i, 2) = 0;
    end
end
end
