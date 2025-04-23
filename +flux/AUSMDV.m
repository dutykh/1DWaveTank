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
%   wL  - [N x 2, double] State matrix [H, HU] on the left side of the interface.
%   wR  - [N x 2, double] State matrix [H, HU] on the right side of the interface.
%   cfg - [struct] Configuration structure. Required fields: cfg.phys.g, cfg.phys.dry_tolerance.
%
% Outputs:
%   F   - [N x 2, double] Numerical flux matrix [F_H, F_HU] across the interface.
%
% Dependencies:
%   None (standalone flux function, but expects correct cfg.phys.g and cfg.phys.dry_tolerance).
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
    g = cfg.phys.g;    % Acceleration due to gravity
    tol = cfg.phys.dry_tolerance;    % Tolerance for dry state check

    % AUSM+ parameters (used in AUSMDV splitting functions)
    alpha_param = 3/16;
    beta_param = 1/8;

    N = size(wL, 1); % Number of interfaces

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Left State (wL) Extraction and Primitive Variable Calculation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hL = wL(:, 1);       % [m] Water depth (left) [N x 1]
    huL = wL(:, 2);      % [m^2/s] Discharge (left) [N x 1]

    % Initialize primitive variables
    uL = zeros(N, 1);
    aL = zeros(N, 1);
    ML = zeros(N, 1);

    % Calculate for wet cells
    wetL = hL >= tol;
    uL(wetL) = huL(wetL) ./ hL(wetL);                % [m/s] Velocity
    aL(wetL) = sqrt(g * hL(wetL));            % [m/s] Speed of sound
    
    % Calculate Mach number where sound speed is non-zero
    calc_ML = wetL & (aL >= tol);
    ML(calc_ML) = uL(calc_ML) ./ aL(calc_ML);             % [unitless] Mach number

    pL = 0.5 * g * hL.^2;                             % [Pa] Physical pressure component [N x 1]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Right State (wR) Extraction and Primitive Variable Calculation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hR = wR(:, 1);       % [m] Water depth (right) [N x 1]
    huR = wR(:, 2);      % [m^2/s] Discharge (right) [N x 1]

    % Initialize primitive variables
    uR = zeros(N, 1);
    aR = zeros(N, 1);
    MR = zeros(N, 1);

    % Calculate for wet cells
    wetR = hR >= tol;
    uR(wetR) = huR(wetR) ./ hR(wetR);                % [m/s] Velocity
    aR(wetR) = sqrt(g * hR(wetR));            % [m/s] Speed of sound

    % Calculate Mach number where sound speed is non-zero
    calc_MR = wetR & (aR >= tol);
    MR(calc_MR) = uR(calc_MR) ./ aR(calc_MR);             % [unitless] Mach number

    pR = 0.5 * g * hR.^2;                             % [Pa] Physical pressure component [N x 1]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interface Speed of Sound (Roe Average)                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h_avg = 0.5 * (hL + hR); % [N x 1]
    a_half = zeros(N, 1);
    calc_a_half = h_avg >= tol;
    a_half(calc_a_half) = sqrt(g * h_avg(calc_a_half)); % [m/s] Roe-averaged speed of sound [N x 1]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUSM+ Style Mach/Pressure Splitting Functions (AUSMDV uses) % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize splitting functions
    M_plus_L = zeros(N, 1);
    P_plus_L = zeros(N, 1);
    M_minus_R = zeros(N, 1);
    P_minus_R = zeros(N, 1);

    % --- Left State Polynomials --- 
    idx_subsonic_L = abs(ML) <= 1.0;
    idx_supersonic_L = ~idx_subsonic_L;

    % Subsonic cases
    M_plus_L(idx_subsonic_L) = 0.25 .* (ML(idx_subsonic_L) + 1).^2 + beta_param .* (ML(idx_subsonic_L).^2 - 1).^2;
    P_plus_L(idx_subsonic_L) = 0.25 .* (ML(idx_subsonic_L) + 1).^2 .* (2 - ML(idx_subsonic_L)) + alpha_param .* ML(idx_subsonic_L) .* (ML(idx_subsonic_L).^2 - 1).^2;
    
    % Supersonic cases
    M_plus_L(idx_supersonic_L) = 0.5 .* (ML(idx_supersonic_L) + abs(ML(idx_supersonic_L)));
    P_plus_L(idx_supersonic_L) = 0.5 .* (1 + sign(ML(idx_supersonic_L)));

    % --- Right State Polynomials --- 
    idx_subsonic_R = abs(MR) <= 1.0;
    idx_supersonic_R = ~idx_subsonic_R;

    % Subsonic cases
    M_minus_R(idx_subsonic_R) = -0.25 .* (MR(idx_subsonic_R) - 1).^2 - beta_param .* (MR(idx_subsonic_R).^2 - 1).^2;
    P_minus_R(idx_subsonic_R) = 0.25 .* (MR(idx_subsonic_R) - 1).^2 .* (2 + MR(idx_subsonic_R)) - alpha_param .* MR(idx_subsonic_R) .* (MR(idx_subsonic_R).^2 - 1).^2;
    
    % Supersonic cases
    M_minus_R(idx_supersonic_R) = 0.5 .* (MR(idx_supersonic_R) - abs(MR(idx_supersonic_R)));
    P_minus_R(idx_supersonic_R) = 0.5 .* (1 - sign(MR(idx_supersonic_R)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUSMDV Interface Velocity (u_p)                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Key difference from AUSM+: interface velocity uses both average and
    % difference of split Mach numbers for improved resolution of contacts.
    u_p = 0.5.*(uL + uR) + 0.5.*( M_plus_L - M_minus_R ) .* a_half; % [m/s] [N x 1]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interface Pressure                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculated using the AUSM+ pressure splitting functions
    p_half = pL .* P_plus_L + pR .* P_minus_R; % [Pa] [N x 1]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mass Flux (F_H) using upwinded convective form with u_p    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Classical upwind form: combines average and jump in H
    F_H = 0.5 .* ( u_p.*(hL + hR) - abs(u_p).*(hR - hL) ); % [m^2/s] [N x 1]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Momentum Flux (F_HU) using upwinded convective form        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F_HU = 0.5 .* ( u_p.*(huL + huR) - abs(u_p).*(huR - huL) ) + p_half; % [m^3/s^2] [N x 1]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Combine and Return Flux Vector                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F = [F_H, F_HU]; % [N x 2]

end