function [wL_interface, wR_interface] = none(w_padded, cfg)
% NONE First-order method (no reconstruction)
%
% Purpose:
%   Implements a first-order method with no reconstruction.
%   Simply copies the cell-centered values to the interfaces.
%   This provides compatibility with the high-order framework while
%   maintaining the original first-order behavior.
%
% Syntax:
%   [wL_interface, wR_interface] = none(w_padded, cfg)
%
% Inputs:
%   w_padded - [(N+2*ng) x 2, double] Padded state array including ghost cells.
%             w_padded(:,1) = water depth [m],
%             w_padded(:,2) = discharge [m^2/s]
%   cfg      - [struct] Configuration structure containing:
%              cfg.mesh.N: Number of cells
%              cfg.bc.num_ghost_cells: Number of ghost cells
%
% Outputs:
%   wL_interface - [(N+1) x 2, double] Left states at all interfaces
%   wR_interface - [(N+1) x 2, double] Right states at all interfaces
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: April 24, 2025

% Extract parameters
ng = cfg.bc.num_ghost_cells;  % Number of ghost cells
N = cfg.mesh.N;               % Number of physical cells

% Preallocate arrays for interface values
wL_interface = zeros(N+1, 2);
wR_interface = zeros(N+1, 2);

% For the first-order method, the interface values are simply the 
% cell-centered values from adjacent cells:
% wL at interface i+1/2 = value from cell i
% wR at interface i+1/2 = value from cell i+1

% Loop over state variables (H and HU)
for var_idx = 1:2
    % Extract the component
    q = w_padded(:, var_idx);
    
    % Compute interface values for all interfaces
    for i = 1:N+1
        idx_padded = i + ng - 1; % Adjust for ghost cells
        
        % Left state at interface i+1/2 is cell-centered value from cell i
        wL_interface(i, var_idx) = q(idx_padded);
        
        % Right state at interface i+1/2 is cell-centered value from cell i+1
        wR_interface(i, var_idx) = q(idx_padded+1);
    end
end
end
