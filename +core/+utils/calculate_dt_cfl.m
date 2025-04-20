function dt = calculate_dt_cfl(w, cfg)
    % CALCULATE_DT_CFL Calculates the stable time step based on CFL condition.
    %   dt = CALCULATE_DT_CFL(w, cfg) calculates the maximum allowed time step
    %   for the 1D Nonlinear Shallow Water Equations given the current state 
    %   vector w = [H; HU] and the configuration structure cfg.
    %
    %   Inputs:
    %       w   - Current state vector (flattened: [H; HU], size 2N x 1).
    %       cfg - Configuration structure containing:
    %             cfg.mesh.N    : Number of grid cells.
    %             cfg.mesh.dx   : Grid cell size.
    %             cfg.param.g   : Acceleration due to gravity.
    %             cfg.time.CFL  : Courant number.
    %
    %   Output:
    %       dt  - Maximum stable time step based on the CFL condition.

    % --- Input Checks and Parameter Extraction ---
    if ~isfield(cfg.time, 'CFL') || cfg.time.CFL <= 0
        error('cfg.time.CFL must be a positive value.');
    end
    if ~isfield(cfg.mesh, 'N') || cfg.mesh.N <= 0
        error('cfg.mesh.N must be a positive integer.');
    end
     if ~isfield(cfg.mesh, 'dx') || cfg.mesh.dx <= 0
        error('cfg.mesh.dx must be a positive value.');
    end
    if ~isfield(cfg.param, 'g') || cfg.param.g <= 0
        error('cfg.param.g must be a positive value.');
    end

    N = cfg.mesh.N;
    dx = cfg.mesh.dx;
    g = cfg.param.g;
    CFL = cfg.time.CFL;

    if length(w) ~= 2 * N
        error('Input state vector w must have length 2*N.');
    end

    % --- Calculate Maximum Wave Speed --- 
    H = w(1:N);         % Extract water depth
    HU = w(N+1:2*N);    % Extract discharge

    % Ensure H is non-negative for sqrt
    H = max(H, 0);

    % Calculate velocity U = HU / H only where H is significantly > 0 (wet cells)
    U = zeros(N, 1);
    wet_tol = 1e-6; % Tolerance for determining wet cells (relative to typical depth or absolute? Let's use absolute for now)
    % Alternatively, could use: wet_tol = 1e-6 * cfg.param.H0 if H0 is always available
    wet_indices = (H > wet_tol);
    
    if any(wet_indices)
        U(wet_indices) = HU(wet_indices) ./ H(wet_indices);
    end

    % Calculate wave speed S = |U| + sqrt(g*H)
    wave_speed = abs(U) + sqrt(g * H);

    % Find the maximum wave speed across the domain
    max_S = max(wave_speed);

    % --- Calculate Time Step --- 
    if max_S < 1e-10 % Handle static or completely dry case
        % If the max speed is effectively zero, any dt is stable theoretically.
        % Return a large value, but perhaps limited by dt_plot or tf/10?
        % Let's return a large number for now, the integrator should cap it.
        dt = inf; 
        % Alternatively: dt = (cfg.time.tf - t_current) / 10; % or some other reasonable large step
    else
        % Calculate dt based on CFL condition
        dt = CFL * dx / max_S;
    end
    
    % Ensure dt is positive (should be unless max_S was negative, which is unlikely)
    dt = max(dt, 0);

end
