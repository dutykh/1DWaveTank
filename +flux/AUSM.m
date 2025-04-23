%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/AUSM.m
%
% Purpose:
%   Computes the AUSM+ (Advection Upstream Splitting Method Plus) numerical
%   flux for the 1D Nonlinear Shallow Water Equations (NSW). This flux is
%   designed to accurately and robustly capture shocks and contact discontinuities
%   by splitting the flux into convective (Mach number) and pressure parts.
%   This implementation is adapted for the shallow water equations and vectorized
%   to accept N x 2 state matrices wL and wR and return an N x 2 flux matrix F.
%
% Syntax:
%   F = AUSM(wL, wR, cfg)
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
    g = cfg.phys.g;    % Acceleration due to gravity
    tol = cfg.phys.dry_tolerance;    % Tolerance for dry state check

    % AUSM+ parameters (tuning for stability and accuracy)
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
    % AUSM+ Mach and Pressure Splitting Polynomials              %
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
    % Interface Convective Mach Number and Velocity               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M_c = M_plus_L + M_minus_R;         % [unitless] Convective Mach number at interface [N x 1]
    % u_c = a_half .* M_c;              % [m/s] Convective velocity at interface [N x 1] (Not explicitly needed for flux)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mass Flux (F_H)                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Standard AUSM+ mass flux: upwinded by Mach splitting
    F_H = a_half .* (M_plus_L .* hL + M_minus_R .* hR); % [m^2/s] [N x 1]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interface Pressure                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pressure splitting polynomials (see Liou, 1996)
    p_half = pL .* P_plus_L + pR .* P_minus_R; % [Pa] [N x 1]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Momentum Flux (F_HU)                                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convective part (Mach splitting) + pressure part
    F_HU = a_half .* (M_plus_L .* huL + M_minus_R .* huR) + p_half; % [m^3/s^2] [N x 1]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Combine and Return Flux Vector                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F = [F_H, F_HU]; % [N x 2]

end