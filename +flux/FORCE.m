function F_num = FORCE(wL, wR, cfg)
    % FORCE Calculate numerical flux using the FORCE scheme.
    %
    %   F_num = FORCE(wL, wR, cfg)
    %
    %   Inputs:
    %       wL   : State vector [hL; qL] at the left side of the interface.
    %       wR   : State vector [hR; qR] at the right side of the interface.
    %       cfg  : Configuration structure (must contain cfg.param.g and cfg.time.CFL).
    %
    %   Outputs:
    %       F_num: Numerical flux vector [Fh; Fq] across the interface.
    %
    %   Method:
    %       Uses the FORCE (First Order Centred) scheme, which is the arithmetic mean
    %       of the Lax-Friedrichs (LF) and Richtmyer (RI) fluxes:
    %       F_FORCE = 0.5 * (F_LF + F_RI)
    %       F_LF = 0.5 * [F(wL) + F(wR) - S_max*(wR - wL)]
    %       F_RI = F(w_half)
    %       w_half = 0.5*(wL + wR) - 0.5*(dt/dx)*(F(wR) - F(wL))
    %       where dt/dx is approximated using CFL as cfg.time.CFL / S_max.
    %
    %   Reference:
    %       Toro, E. F. (2009). Riemann Solvers and Numerical Methods for 
    %       Fluid Dynamics: A Practical Introduction. Springer.
    %       (Chapter 6)
    %
    %   Author: Denys Dutykh
    %   Date:   20 April 2025

    persistent g cfl % Store gravity and CFL locally for efficiency

    % Initialize persistent gravity 'g' and CFL number 'cfl' if they are empty
    if isempty(g) || isempty(cfl)
        % Get gravity 'g'
        if isfield(cfg, 'param') && isfield(cfg.param, 'g')
            g = cfg.param.g;
        else
            g = 9.81; % Default gravity
            warning('FORCE: Using default g = 9.81 m/s^2');
        end
        % Get CFL number
        if isfield(cfg, 'time') && isfield(cfg.time, 'CFL')
            cfl = cfg.time.CFL;
        else
            cfl = 0.9; % Default CFL
            warning('FORCE: Using default CFL = 0.9');
        end
    end

    h_eps = 1e-6; % Tolerance for dry state

    % --- Physical Flux Function --- 
    function F = physical_flux(w_vec, g_val)
        h = w_vec(1);
        q = w_vec(2);
        if h <= h_eps
            u = 0;
            F = [0; 0];
        else
            u = q / h;
            F = [q; q*u + 0.5*g_val*h^2];
        end
    end

    % --- Calculate Left/Right States and Speeds --- 
    wL = wL(:); % [HL; qL]
    wR = wR(:); % [HR; qR]
    hL = wL(1); qL = wL(2);
    hR = wR(1); qR = wR(2);

    uL = 0; cL = 0;
    if hL > h_eps
        uL = qL / hL;
        cL = sqrt(g * hL);
    end

    uR = 0; cR = 0;
    if hR > h_eps
        uR = qR / hR;
        cR = sqrt(g * hR);
    end

    % Max wave speed estimate
    S_max = max(abs(uL) + cL, abs(uR) + cR);
    if S_max < 1e-9
        S_max = 1e-9; % Avoid division by zero if flow is stagnant
    end

    % --- Calculate Fluxes --- 
    F_L = physical_flux(wL, g);
    F_R = physical_flux(wR, g);

    % 1. Lax-Friedrichs Flux
    F_LF = 0.5 * (F_L + F_R - S_max * (wR(:) - wL(:)));

    % 2. Richtmyer Flux
    % Intermediate state w_half
    dtdx_approx = cfl / S_max;
    w_half = 0.5 * (wL(:) + wR(:)) - 0.5 * dtdx_approx * (F_R - F_L);
    
    % Ensure positivity of h_half (simple clipping)
    w_half(1) = max(w_half(1), h_eps);
    
    % Physical flux at intermediate state
    F_RI = physical_flux(w_half, g);

    % 3. FORCE Flux (Average)
    F_num = 0.5 * (F_LF + F_RI);

    % Ensure output is ROW vector
    F_num = F_num(:).'; 

end
