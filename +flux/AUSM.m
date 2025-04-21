function F = AUSM(wL, wR, cfg)

    %AUSM+ Advection Upstream Splitting Method Plus (AUSM+) numerical flux for SWE.
    %
    %   F = AUSM(WL, WR, CFG) calculates the numerical flux vector across an
    %   interface using the AUSM+ scheme for the 1D Nonlinear Shallow Water
    %   Equations (NSW).
    %
    %   Inputs:
    %       wL  - State vector [H, HU] on the left side of the interface.
    %       wR  - State vector [H, HU] on the right side of the interface.
    %       cfg - Configuration structure containing physical parameters (cfg.phys.g).
    %
    %   Outputs:
    %       F   - Numerical flux vector [F_H, F_HU] across the interface (1x2).
    %
    %   References:
    %     - Liou, M.-S. (1996). A sequel to AUSM: AUSM+.
    %       Journal of Computational Physics, 129(2), 364-382.
    %     - Adapted for Shallow Water Equations.
    %
    %   Author: Denys Dutykh
    %   Date:   21 April 2025

    g = cfg.phys.g; % Acceleration due to gravity
    tol = 1e-10;    % Tolerance for dry state check

    % AUSM+ parameters (common values)
    alpha_param = 3/16;
    beta_param = 1/8;

    % --- Left State --- 
    hL = wL(1);
    huL = wL(2);

    % Check for dry state
    if hL < tol
        hL = 0;
        huL = 0;
        uL = 0;
        aL = 0;
        ML = 0;
    else
        uL = huL / hL;
        aL = sqrt(g * hL);
        if aL < tol
            ML = 0; % Avoid division by zero if aL is effectively zero
        else
            ML = uL / aL;
        end
    end
    pL = 0.5 * g * hL^2; % Physical pressure component

    % --- Right State --- 
    hR = wR(1);
    huR = wR(2);

    % Check for dry state
    if hR < tol
        hR = 0;
        huR = 0;
        uR = 0;
        aR = 0;
        MR = 0;
    else
        uR = huR / hR;
        aR = sqrt(g * hR);
         if aR < tol
            MR = 0; % Avoid division by zero if aR is effectively zero
        else
            MR = uR / aR;
         end
    end
    pR = 0.5 * g * hR^2; % Physical pressure component

    % --- Interface Speed of Sound (using Roe average) --- 
    h_avg = 0.5 * (hL + hR);
    if h_avg < tol
        a_half = 0; % Avoid sqrt of negative/small number if average depth is zero
    else
        a_half = sqrt(g * h_avg); % Roe-averaged speed of sound
    end
    
    % Prevent division by zero if a_half is zero
    if a_half < tol
        a_half_inv = 0;
    else
        a_half_inv = 1.0 / a_half;
    end

    % --- AUSM+ Mach Number Splitting --- 
    if abs(ML) <= 1.0
        M_plus_L = 0.25 * (ML + 1)^2 + beta_param * (ML^2 - 1)^2;
        P_plus_L = 0.25 * (ML + 1)^2 * (2 - ML) + alpha_param * ML * (ML^2 - 1)^2;
    else
        M_plus_L = 0.5 * (ML + abs(ML));
        P_plus_L = 0.5 * (1 + sign(ML)); % Equivalent to M_plus_L / ML for |M| >= 1
    end

    if abs(MR) <= 1.0
        M_minus_R = -0.25 * (MR - 1)^2 - beta_param * (MR^2 - 1)^2;
        P_minus_R = 0.25 * (MR - 1)^2 * (2 + MR) - alpha_param * MR * (MR^2 - 1)^2;
    else
        M_minus_R = 0.5 * (MR - abs(MR));
        P_minus_R = 0.5 * (1 - sign(MR)); % Equivalent to M_minus_R / MR for |M| >= 1
    end

    % --- Interface Convective Mach Number & Velocity --- 
    M_c = M_plus_L + M_minus_R; % Common convective Mach number
    u_c = a_half * M_c;         % Common convective velocity

    % --- Mass Flux (F_H) --- 
    F_H = a_half * (M_plus_L * hL + M_minus_R * hR);
    % Alternatively, using common velocity u_c:
    % F_H = u_c * 0.5 * (hL + hR) - 0.5 * abs(u_c) * (hR - hL) % (Needs verification for AUSM+)
    % The first form F_H = a_half * (...) is more standard for AUSM+ mass flux

    % --- Interface Pressure --- 
    p_half = pL * P_plus_L + pR * P_minus_R;

    % --- Momentum Flux (F_HU) --- 
    % F_HU = a_half * (M_plus_L * huL + M_minus_R * huR) + p_half
    % F_HU = 0.5*(F_H(hL uL+hR uR) + |F_H|(hL uL-hR uR))+ p_half % This is not the correct way?
    F_HU = a_half * (M_plus_L * huL + M_minus_R * huR) + p_half;
    % Alternative using u_c:
    % F_HU = u_c * 0.5 * (huL + huR) - 0.5 * abs(u_c) * (huR - huL) + p_half % (Needs verification)
    % Stick to the standard AUSM+ form

    % --- Final Flux Vector --- 
    F = [F_H, F_HU];

end