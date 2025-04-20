% +bc/open.m
function w_padded = open(w_padded, t, side, cfg, num_ghost_cells)

    %OPEN Implements an open (zeroth-order extrapolation) boundary condition.
    %   Copies the state from the nearest interior cell(s) to the ghost cell(s).
    %   w_padded = OPEN(w_padded, t, side, cfg, num_ghost_cells)
    %
    %   Inputs:
    %       w_padded        - State vector array including space for ghost cells.
    %       t               - Current time (unused).
    %       side            - String, either 'left' or 'right'.
    %       cfg             - Configuration structure (unused).
    %       num_ghost_cells - Number of ghost cells to fill.
    %
    %   Outputs:
    %       w_padded        - State vector array with ghost cells filled.

    N = cfg.mesh.N; % Number of actual cells

    if strcmp(side, 'left')
        % Left boundary (indices 1:num_ghost_cells)
        interior_cell_idx = num_ghost_cells + 1; % Index of the first interior cell
        for i = 1:num_ghost_cells
            ghost_idx = i;
            w_padded(ghost_idx, :) = w_padded(interior_cell_idx, :);
        end
    elseif strcmp(side, 'right')
        % Right boundary (indices N+num_ghost_cells+1 : N+2*num_ghost_cells)
        interior_cell_idx = N + num_ghost_cells; % Index of the last interior cell
        for i = 1:num_ghost_cells
            ghost_idx = interior_cell_idx + i;
             w_padded(ghost_idx, :) = w_padded(interior_cell_idx, :);
        end
    else
        error('Invalid side specified. Use ''left'' or ''right''.');
    end

end