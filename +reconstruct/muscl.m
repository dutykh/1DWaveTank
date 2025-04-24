function [wL_interface, wR_interface] = muscl(w_padded, cfg)
    % Extract parameters
    ng = cfg.bc.num_ghost_cells;
    N = cfg.mesh.N;
    dry_tolerance = cfg.phys.dry_tolerance;
    total_cells_padded = N + 2*ng;

    % Get theta parameter - Note: theta is typically used in specific slope
    % calculation formulas (like ULTIMATE), not directly with standard limiters.
    % Standard MUSCL usually just computes limited slopes from fwd/bwd diffs.
    if ~isfield(cfg.reconstruct, 'theta') || isempty(cfg.reconstruct.theta)
        theta = 0.0; % Defaulting to 0, but not used in slope calc below
    else
        theta = cfg.reconstruct.theta;
    end

    % Get limiter
    if ~isfield(cfg.reconstruct, 'limiter') || isempty(cfg.reconstruct.limiter)
        limiter_handle = @reconstruct.limiters.minmod;
        warning('MUSCL:NoLimiter', 'No slope limiter specified, using minmod by default.');
    else
        limiter_handle = cfg.reconstruct.limiter;
    end

    % Initialize arrays
    wL_interface = zeros(N+1, 2);
    wR_interface = zeros(N+1, 2);
    slopes = zeros(total_cells_padded, 2); % Store slopes for H and HU

    % --- Step 1: Calculate Limited Slopes for Conservative Variables (H, HU) --- 
    for var_idx = 1:2 % Loop over H (1) and HU (2)
        q = w_padded(:, var_idx);
        current_var_slopes = zeros(total_cells_padded, 1);

        % Loop over all cells where slopes can be computed (indices 2 to total_cells-1)
        for i = 2:(total_cells_padded-1)
            % Forward difference (q_{i+1} - q_i)
            delta_plus = q(i+1) - q(i);

            % Backward difference (q_i - q_{i-1})
            delta_minus = q(i) - q(i-1);

            % Apply slope limiter to forward/backward differences
            current_var_slopes(i) = limiter_handle(delta_minus, delta_plus);
        end
        slopes(:, var_idx) = current_var_slopes;
    end

    % --- Step 2: Reconstruct Conservative Variables at Interfaces --- 
    % Loop over all interior interfaces (1 to N+1)
    for i = 1:N+1
        % Index of the cell to the LEFT of interface i+1/2 in the padded array
        idx_left = i + ng - 1;
        % Index of the cell to the RIGHT of interface i+1/2 in the padded array
        idx_right = i + ng; % same as idx_left + 1

        for var_idx = 1:2 % Reconstruct H and HU
            q = w_padded(:, var_idx);
            current_slopes = slopes(:, var_idx);

            % Left state at interface i+1/2: q_i + 0.5*slope_i
            wL_interface(i, var_idx) = q(idx_left) + 0.5 * current_slopes(idx_left);

            % Right state at interface i+1/2: q_{i+1} - 0.5*slope_{i+1}
            wR_interface(i, var_idx) = q(idx_right) - 0.5 * current_slopes(idx_right);
        end
    end

    % --- Step 3: Ensure Positivity and Handle Dry States --- 
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
