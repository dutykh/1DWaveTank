%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +bc/periodic.m
%
% Purpose:
%   Applies periodic boundary conditions to the padded state vector.
%   Copies data from the interior domain at one end to the ghost cells
%   at the opposite end.
%
% Syntax:
%   w_padded = periodic(w_padded, t, side, cfg, num_ghost_cells)
%
% Inputs:
%   w_padded        - [N + 2*num_ghost_cells, 2, double] Padded state vector [H, HU].
%   t               - [double] Current time (unused for periodic BC).
%   side            - [char] Specifies which side ('left', 'right', or 'both').
%                     Typically called once with 'both' or internally handled.
%   cfg             - [struct] Configuration structure (contains N).
%   num_ghost_cells - [integer] Number of ghost cells on each side.
%
% Outputs:
%   w_padded        - [N + 2*num_ghost_cells, 2, double] State vector with
%                     periodic boundary conditions applied.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   22 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w_padded = periodic(w_padded, t, side, cfg, num_ghost_cells) %#ok<INUSD> t, side unused
    % Assumes w_padded is [N + 2*Ng, 2] where column 1 is H, column 2 is HU
    N = cfg.mesh.N;
    Ng = num_ghost_cells;
    N_total_padded_rows = N + 2*Ng;

    % Indices for ghost and corresponding interior cells
    left_ghost_idx = 1:Ng;
    right_interior_idx = N + 1 : N + Ng;

    right_ghost_idx = N + Ng + 1 : N_total_padded_rows;
    left_interior_idx = Ng + 1 : Ng + Ng; 

    % Apply periodic BCs to H (column 1)
    w_padded(left_ghost_idx, 1) = w_padded(right_interior_idx, 1);
    w_padded(right_ghost_idx, 1) = w_padded(left_interior_idx, 1);

    % Apply periodic BCs to HU (column 2)
    w_padded(left_ghost_idx, 2) = w_padded(right_interior_idx, 2);
    w_padded(right_ghost_idx, 2) = w_padded(left_interior_idx, 2);

    % Return the modified w_padded [N+2*Ng, 2]
end
