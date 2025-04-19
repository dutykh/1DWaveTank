function w_padded = placeholder_bc_1st_order(w_padded, t, side, cfg, num_ghost_cells)

    %PLACEHOLDER_BC_1ST_ORDER Dummy BC for 1st order, adds 1 ghost cell.
    %   w_padded = PLACEHOLDER_BC_1ST_ORDER(w_padded, t, side, cfg, num_ghost_cells)
    %   fills the single ghost cell needed for the 1st order RHS function by
    %   copying the value from the nearest interior cell (zeroth-order extrapolation).
    %
    %   Inputs:
    %       w_padded        - State vector array including space for ghost cells.
    %       t               - Current time (unused).
    %       side            - String, either 'left' or 'right'.
    %       cfg             - Configuration structure (unused).
    %       num_ghost_cells - Number of ghost cells (should be 1 for this).
    %
    %   Outputs:
    %       w_padded        - State vector array with the specified ghost cell filled.
    
    N = cfg.mesh.N; % Number of actual cells
    
    if num_ghost_cells ~= 1
        error('This placeholder BC is designed for exactly 1 ghost cell.');
    end
    
    warning('Using placeholder boundary condition on %s side (zeroth-order extrapolation).', side);
    
    if strcmp(side, 'left')
        % Fill ghost cell at index 1 by copying from cell 1 (index 2 in padded array)
        w_padded(1, :) = w_padded(2, :);
    elseif strcmp(side, 'right')
        % Fill ghost cell at index N+2 by copying from cell N (index N+1 in padded array)
        w_padded(N+2, :) = w_padded(N+1, :);
    else
        error('Invalid side specified. Use ''left'' or ''right''.');
    end

end