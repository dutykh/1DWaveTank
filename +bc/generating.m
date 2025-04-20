% +bc/generating.m
function w_padded = generating(w_padded, t, side, cfg, num_ghost_cells)
    %GENERATING Implements a wave generating boundary condition.
    %   Sets H in the ghost cell based on a time function (e.g., sine wave)
    %   and computes HU based on characteristic theory (simplified here).
    %   w_padded = GENERATING(w_padded, t, side, cfg, num_ghost_cells)
    %
    %   Requires parameters in cfg:
    %       cfg.bc.[side].param.a (amplitude, default 0.1)
    %       cfg.bc.[side].param.T (period, default 2*pi)
    %       cfg.phys.g (gravity)
    %       cfg.param.H0 (still water depth, needed for c0)
    %
    %   Inputs:
    %       w_padded        - State vector array including space for ghost cells.
    %       t               - Current time.
    %       side            - String, either 'left' or 'right'.
    %       cfg             - Configuration structure.
    %       num_ghost_cells - Number of ghost cells to fill (typically 1 for this BC).
    %
    %   Outputs:
    %       w_padded        - State vector array with ghost cells filled.

    if num_ghost_cells ~= 1
        warning('Generating BC typically implemented for 1 ghost cell. Using only the outermost.');
    end

    g = cfg.phys.g;
    N = cfg.mesh.N;

    % --- Get parameters for BoundaryValue function ---
    if strcmp(side,'left')
        params = cfg.bc.left.param;
        interior_cell_idx = num_ghost_cells + 1;
        ghost_cell_idx = num_ghost_cells; % Outermost ghost cell
    elseif strcmp(side,'right')
         params = cfg.bc.right.param;
         interior_cell_idx = N + num_ghost_cells;
         ghost_cell_idx = N + num_ghost_cells + 1; % Outermost ghost cell
    else
        error('Invalid side specified. Use ''left'' or ''right''.');
    end

    a = params.a; % Amplitude
    T = params.T; % Period
    H0 = cfg.param.H0; % Reference depth for characteristics

    % --- BoundaryValue function (Nested) ---
    function H_boundary = BoundaryValue(time)
        % Calculates target water depth at the boundary at a given time.
         omega = 2*pi/T;
         H_boundary = H0 + a * sin(omega * time);
         H_boundary = max(H_boundary, 1e-6); % Ensure non-negative depth
    end
    % --- End of BoundaryValue ---

    % Set H in the ghost cell
    H_ghost = BoundaryValue(t);
    w_padded(ghost_cell_idx, 1) = H_ghost;

    % --- Calculate HU in the ghost cell using Riemann invariants (simplified) ---
    % This assumes subcritical flow entering the domain from the boundary.
    % We use the outward propagating characteristic from the interior cell
    % and the inward propagating characteristic determined by H_ghost.
    % C- = u - 2*c
    % C+ = u + 2*c
    % For left boundary (inflow): C- from interior, C+ from boundary H_ghost
    % For right boundary (inflow): C+ from interior, C- from boundary H_ghost

    H_interior = w_padded(interior_cell_idx, 1);
    HU_interior = w_padded(interior_cell_idx, 2);
    if H_interior > 1e-6
         U_interior = HU_interior / H_interior;
         c_interior = sqrt(g * H_interior);
    else
         U_interior = 0;
         c_interior = 0;
    end

    c_ghost = sqrt(g * H_ghost);

    if strcmp(side, 'left')
        % Riemann invariant C- from interior cell
        Riemann_minus = U_interior - 2 * c_interior;
        % Calculate U_ghost using C- and c_ghost (derived from H_ghost)
        % U_ghost - 2 * c_ghost = Riemann_minus  => U_ghost = Riemann_minus + 2*c_ghost
        U_ghost = Riemann_minus + 2 * c_ghost;
        % Ensure inflow condition (U_ghost > 0 for left boundary)
        % Note: This simple approach might not be robust for all cases.
        % More sophisticated methods exist (e.g., solving the Riemann problem).
         U_ghost = max(U_ghost, 0); % Simplified: prevent outflow if generating

    elseif strcmp(side, 'right')
        % Riemann invariant C+ from interior cell
        Riemann_plus = U_interior + 2 * c_interior;
        % Calculate U_ghost using C+ and c_ghost (derived from H_ghost)
        % U_ghost + 2 * c_ghost = Riemann_plus => U_ghost = Riemann_plus - 2*c_ghost
        U_ghost = Riemann_plus - 2 * c_ghost;
        % Ensure inflow condition (U_ghost < 0 for right boundary)
         U_ghost = min(U_ghost, 0); % Simplified: prevent outflow if generating
    end

    % Set HU in the ghost cell
    w_padded(ghost_cell_idx, 2) = H_ghost * U_ghost;

    % --- Extrapolate to other ghost cells if num_ghost_cells > 1 ---
    % Simplest: copy the state from the outermost computed ghost cell
    if num_ghost_cells > 1
        if strcmp(side,'left')
             for i = 1:(num_ghost_cells-1)
                  w_padded(i,:) = w_padded(ghost_cell_idx,:);
             end
        elseif strcmp(side,'right')
             for i = (ghost_cell_idx+1):(N + 2*num_ghost_cells)
                  w_padded(i,:) = w_padded(ghost_cell_idx,:);
             end
        end
    end

end