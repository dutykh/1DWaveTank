% +bc/wall.m
function w_padded = wall(w_padded, t, side, cfg, num_ghost_cells)

    %WALL Implements a solid wall boundary condition using ghost cells.
    %   Sets ghost cell values such that velocity at the boundary is zero.
    %   w_padded = WALL(w_padded, t, side, cfg, num_ghost_cells)
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
        interior_idx_start = num_ghost_cells + 1;
        for i = 1:num_ghost_cells
            ghost_idx = num_ghost_cells - i + 1;
            interior_idx = interior_idx_start + i - 1;
            w_padded(ghost_idx, 1) =  w_padded(interior_idx, 1);  % H_ghost = H_interior
            w_padded(ghost_idx, 2) = -w_padded(interior_idx, 2);  % HU_ghost = -HU_interior
        end
    elseif strcmp(side, 'right')
        % Right boundary (indices N+num_ghost_cells+1 : N+2*num_ghost_cells)
        interior_idx_end = N + num_ghost_cells;
        for i = 1:num_ghost_cells
            ghost_idx = interior_idx_end + i;
            interior_idx = interior_idx_end - i + 1;
            w_padded(ghost_idx, 1) =  w_padded(interior_idx, 1);  % H_ghost = H_interior
            w_padded(ghost_idx, 2) = -w_padded(interior_idx, 2);  % HU_ghost = -HU_interior
        end
    else
        error('Invalid side specified. Use ''left'' or ''right''.');
    end

end