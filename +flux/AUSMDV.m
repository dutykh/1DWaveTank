%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/AUSMDV.m
%
% Purpose:
%   Computes the AUSMDV (AUSM-Derivative Variant) numerical flux for the
%   1D Nonlinear Shallow Water Equations (NSW). AUSMDV is a variant of AUSM+
%   designed for improved shock and contact resolution, using a derivative-based
%   interface velocity and upwind flux forms. This implementation is adapted
%   for shallow water equations.
%
% Syntax:
%   F = AUSMDV(wL, wR, cfg)
%
% Inputs:
%   wL  - [1 x 2, double] State vector [H, HU] on the left side of the interface.
%   wR  - [1 x 2, double] State vector [H, HU] on the right side of the interface.
%   cfg - [struct] Configuration structure. Required field: cfg.phys.g (gravity).
%
% Outputs:
%   F   - [1 x 2, double] Numerical flux vector [F_H, F_HU] across the interface.
%
% Dependencies:
%   None (standalone flux function, but expects correct cfg.phys.g).
%
% References:
%   - Wada, Y., & Liou, M. S. (1997). A flux splitting scheme with high-resolution
%     and robustness for discontinuities. AIAA paper, 94-0083.
%   - Adapted for Shallow Water Equations (see e.g. Toro, "Riemann Solvers...").
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = AUSMDV(wL, wR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameter and Tolerance Setup                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = cfg.phys.g; % [m/s^2] Acceleration due to gravity
    tol = 1e-10;    % Tolerance for dry state check

    % AUSM+ parameters (used in AUSMDV splitting functions)
    alpha_param = 3/16;
    beta_param = 1/8;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Left State (wL) Extraction and Mach Calculation             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hL = wL(1);       % [m] Water depth (left)
    huL = wL(2);      % [m^2/s] Discharge (left)

    % Handle dry state: set all variables to zero if hL < tol
    if hL < tol
        hL = 0; huL = 0; uL = 0; aL = 0; ML = 0;
    else
        uL = huL / hL;                 % [m/s] Velocity
        aL = sqrt(g * hL);             % [m/s] Speed of sound
        if aL < tol
            ML = 0; % Avoid division by zero
        else
            ML = uL / aL;              % [unitless] Mach number
        end
    end
    pL = 0.5 * g * hL^2;              % [Pa] Physical pressure component

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Right State (wR) Extraction and Mach Calculation            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hR = wR(1);       % [m] Water depth (right)
    huR = wR(2);      % [m^2/s] Discharge (right)

    % Handle dry state: set all variables to zero if hR < tol
    if hR < tol
        hR = 0; huR = 0; uR = 0; aR = 0; MR = 0;
    else
        uR = huR / hR;                 % [m/s] Velocity
        aR = sqrt(g * hR);             % [m/s] Speed of sound
        if aR < tol
            MR = 0; % Avoid division by zero
        else
            MR = uR / aR;              % [unitless] Mach number
        end
    end
    pR = 0.5 * g * hR^2;              % [Pa] Physical pressure component

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interface Speed of Sound (Roe Average)                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Roe-averaged speed of sound for stability and well-balancing
    h_avg = 0.5 * (hL + hR);
    if h_avg < tol
        a_half = 0; % Avoid sqrt of negative/small number if average depth is zero
    else
        a_half = sqrt(g * h_avg); % [m/s] Roe-averaged speed of sound
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUSM+ Style Mach/Pressure Splitting Functions (AUSMDV uses) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUSMDV Interface Velocity (u_p)                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Key difference from AUSM+: interface velocity uses both average and
    % difference of split Mach numbers for improved resolution of contacts.
    u_p = 0.5*(uL + uR) + 0.5*( M_plus_L - M_minus_R ) * a_half; % [m/s]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interface Pressure                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculated using the AUSM+ pressure splitting functions
    p_half = pL * P_plus_L + pR * P_minus_R; % [Pa]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mass Flux (F_H) using upwinded convective form with u_p    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Classical upwind form: combines average and jump in H
    F_H = 0.5 * ( u_p*(hL + hR) - abs(u_p)*(hR - hL) ); % [m^2/s]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Momentum Flux (F_HU) using upwinded convective form        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F_HU = 0.5 * ( u_p*(huL + huR) - abs(u_p)*(huR - huL) ) + p_half; % [m^3/s^2]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Combine and Return Flux Vector                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F = [F_H, F_HU];

end