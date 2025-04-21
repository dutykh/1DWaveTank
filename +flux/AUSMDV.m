function F = AUSMDV(wL, wR, cfg)

    %AUSMDV AUSM-Derivative Variant numerical flux for SWE.
    %
    %   F = AUSMDV(WL, WR, CFG) calculates the numerical flux vector across an
    %   interface using the AUSMDV scheme for the 1D Nonlinear Shallow Water
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
    %     - Wada, Y., & Liou, M. S. (1997). A flux splitting scheme with high-resolution
    %       and robustness for discontinuities. AIAA paper, 94-0083.
    %     - Adapted for Shallow Water Equations.
    %
    %   Author: Denys Dutykh
    %   Date:   21 April 2025

    g = cfg.phys.g; % Acceleration due to gravity
    tol = 1e-10;    % Tolerance for dry state check

    % AUSM+ parameters used by AUSMDV splitting functions
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
    
    % --- AUSM+ Style Mach/Pressure Splitting Functions (used by AUSMDV) --- 
    if abs(ML) <= 1.0
        M_plus_L = 0.25 * (ML + 1)^2 + beta_param * (ML^2 - 1)^2;
        P_plus_L = 0.25 * (ML + 1)^2 * (2 - ML) + alpha_param * ML * (ML^2 - 1)^2;
    else
        M_plus_L = 0.5 * (ML + abs(ML));
        P_plus_L = 0.5 * (1 + sign(ML)); 
    end

    if abs(MR) <= 1.0
        M_minus_R = -0.25 * (MR - 1)^2 - beta_param * (MR^2 - 1)^2;
        P_minus_R = 0.25 * (MR - 1)^2 * (2 + MR) - alpha_param * MR * (MR^2 - 1)^2;
    else
        M_minus_R = 0.5 * (MR - abs(MR));
        P_minus_R = 0.5 * (1 - sign(MR)); 
    end

    % --- AUSMDV Interface Velocity (u_p) --- 
    % Note: This is a key difference from AUSM+
    u_p = 0.5*(uL + uR) + 0.5*( M_plus_L - M_minus_R ) * a_half;
    
    % --- Interface Pressure --- 
    % Calculated using the AUSM+ pressure splitting functions
    p_half = pL * P_plus_L + pR * P_minus_R;

    % --- Mass Flux (F_H) using convective form with u_p --- 
    F_H = 0.5 * ( u_p*(hL + hR) - abs(u_p)*(hR - hL) );
    
    % --- Momentum Flux (F_HU) using convective form with u_p --- 
    F_HU = 0.5 * ( u_p*(huL + huR) - abs(u_p)*(huR - huL) ) + p_half;
    
    % --- Final Flux Vector --- 
    F = [F_H, F_HU];

end
