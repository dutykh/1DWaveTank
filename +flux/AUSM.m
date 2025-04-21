%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/AUSM.m
%
% Purpose:
%   Computes the AUSM+ (Advection Upstream Splitting Method Plus) numerical
%   flux for the 1D Nonlinear Shallow Water Equations (NSW). This flux is
%   designed to accurately and robustly capture shocks and contact discontinuities
%   by splitting the flux into convective (Mach number) and pressure parts.
%   This implementation is adapted for the shallow water equations.
%
% Syntax:
%   F = AUSM(wL, wR, cfg)
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
%   - Liou, M.-S. (1996). A sequel to AUSM: AUSM+.
%     Journal of Computational Physics, 129(2), 364-382.
%   - Adapted for Shallow Water Equations (see e.g. Toro, "Riemann Solvers...").
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = AUSM(wL, wR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameter and Tolerance Setup                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = cfg.phys.g; % [m/s^2] Acceleration due to gravity
    tol = 1e-10;    % Tolerance for dry state check

    % AUSM+ parameters (tuning for stability and accuracy)
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
    
    % Prevent division by zero if a_half is zero (not used, but kept for clarity)
    if a_half < tol
        a_half_inv = 0;
    else
        a_half_inv = 1.0 / a_half;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUSM+ Mach and Pressure Splitting Polynomials              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % These polynomials are designed for sharp shock/contact capturing
    % and to avoid spurious oscillations. See Liou (1996) for details.
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interface Convective Mach Number and Velocity               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M_c = M_plus_L + M_minus_R; % [unitless] Convective Mach number at interface
    u_c = a_half * M_c;         % [m/s] Convective velocity at interface

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mass Flux (F_H)                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Standard AUSM+ mass flux: upwinded by Mach splitting
    F_H = a_half * (M_plus_L * hL + M_minus_R * hR); % [m^2/s]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interface Pressure                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pressure splitting polynomials (see Liou, 1996)
    p_half = pL * P_plus_L + pR * P_minus_R; % [Pa]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Momentum Flux (F_HU)                                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convective part (Mach splitting) + pressure part
    F_HU = a_half * (M_plus_L * huL + M_minus_R * huR) + p_half; % [m^3/s^2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Combine and Return Flux Vector                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F = [F_H, F_HU];

end