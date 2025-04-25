%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/Kinetic.m
%
% Purpose:
%   Computes the Kinetic numerical flux for the 1D Non-Linear Shallow Water
%   Equations (NSW). This implementation is based on the kinetic scheme with
%   indicator function from the original Fortran code in schema2.f.
%
% Syntax:
%   F = Kinetic(wL, wR, cfg)
%
% Inputs:
%   wL  - [N x 2, double] State vector [H, HU] on left side of interface.
%   wR  - [N x 2, double] State vector [H, HU] on right side of interface.
%   cfg - [struct] Configuration structure. Required fields:
%          cfg.phys.g: [double] Gravitational acceleration [m/s^2].
%          cfg.phys.dry_tolerance: [double] Threshold for dry cells [m].
%
% Outputs:
%   F   - [N x 2, double] Kinetic numerical flux vector [F_H, F_HU].
%
% References:
%   - Original Fortran implementation in schema2.f (calcul_s2, F1, F2M, F2P)
%
% Acknowledgment:
%   The kinetic flux algorithm implemented here is based on original Fortran code
%   contributed by Professor Mehmet ERSOY (SEATECH - École d'Ingénieurs de l'Université de Toulon,
%   IMATH - Institut de Mathématiques de Toulon, France).
%   Contact:
%     Prof. Mehmet ERSOY
%     Avenue de l’université BP 20132, 83957 La Garde Cedex, France
%     Tél.: +33 483 166 665 | Port.: +33 672 533 633
%     Web: http://ersoy.univ-tln.fr/
%
%   We gratefully acknowledge Prof. ERSOY's contribution of the original kinetic flux routine.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: 23 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = Kinetic(wL, wR, cfg)

    % Extract parameters
    g = cfg.phys.g;        % Gravitational acceleration
    dry_tol = cfg.phys.dry_tolerance; % Threshold for dry cells
    
    % Extract states on the left and right
    hL = wL(:,1); % Water depth on the left
    qL = wL(:,2); % Discharge on the left
    hR = wR(:,1); % Water depth on the right
    qR = wR(:,2); % Discharge on the right
    
    % Ensure strictly positive water depths for numerical stability
    hL = max(hL, dry_tol);
    hR = max(hR, dry_tol);
    
    % Initialize the flux vector
    F = zeros(size(wL));
    
    % Compute bathymetry difference
    % In the 1DWaveTank code, bathymetry isn't directly accessible in flux functions
    % For now, assume flat bathymetry (dZ = 0)
    dZ = zeros(size(hL));
    
    % First component (mass conservation)
    F(:,1) = kinetic_F1(hL, qL, hR, qR, dZ, g, cfg, dry_tol);
    
    % Second component (momentum conservation)
    % In the original Fortran code, this used F2M for the interface flux
    F(:,2) = kinetic_F2M(hL, qL, hR, qR, dZ, g, cfg, dry_tol);

end

function f1 = kinetic_F1(hL, qL, hR, qR, dZ, g, cfg, dry_tol)
    epsilon = cfg.numerics.epsilon;

    % Compute the first component of the kinetic flux
    % Vectorized implementation
    
    % Compute left contribution
    sonmL = sqrt(g * hL / 2.0);
    uL = qL ./ hL;
    borneL = sqrt(2.0 * g * pos(dZ));
    bgL = max(borneL, uL - sonmL .* sqrt(3.0));
    bdL = max(borneL, uL + sonmL .* sqrt(3.0));
    
    % Avoid division by zero
    denomL = 4.0 * sonmL .* sqrt(3.0);
    idx_validL = denomL > epsilon;
    tmp1 = zeros(size(hL));
    tmp1(idx_validL) = hL(idx_validL) .* (bdL(idx_validL).^2 - bgL(idx_validL).^2) ./ denomL(idx_validL);
    
    % Compute right contribution
    sonmR = sqrt(g * hR / 2.0);
    uR = qR ./ hR;
    borneR = -sqrt(2.0 * g * neg(dZ));
    bgR = min(borneR, uR - sonmR .* sqrt(3.0));
    bdR = min(borneR, uR + sonmR .* sqrt(3.0));
    
    % Avoid division by zero
    denomR = 4.0 * sonmR .* sqrt(3.0);
    idx_validR = denomR > epsilon;
    tmp2 = zeros(size(hR));
    tmp2(idx_validR) = hR(idx_validR) .* (bdR(idx_validR).^2 - bgR(idx_validR).^2) ./ denomR(idx_validR);
    
    % Combine
    f1 = tmp1 + tmp2;

end

function f2m = kinetic_F2M(hL, qL, hR, qR, dZ, g, cfg, dry_tol)
    epsilon = cfg.numerics.epsilon;

    % Compute the second component of the kinetic flux (from left to right)
    % Vectorized implementation
    
    % Left contribution 1
    sonmL = sqrt(g * hL / 2.0);
    uL = qL ./ hL;
    bgL1 = pos(uL - sonmL .* sqrt(3.0));
    bdL1 = pos(uL + sonmL .* sqrt(3.0));
    
    % Avoid division by zero
    denomL = 6.0 * sonmL .* sqrt(3.0);
    idx_validL = denomL > epsilon;
    tmp1 = zeros(size(hL));
    tmp1(idx_validL) = hL(idx_validL) .* (bdL1(idx_validL).^3 - bgL1(idx_validL).^3) ./ denomL(idx_validL);
    
    % Left contribution 2
    borneL = sqrt(2.0 * g * pos(dZ));
    bgL2 = min(borneL, pos(uL - sonmL .* sqrt(3.0)));
    bdL2 = min(borneL, pos(uL + sonmL .* sqrt(3.0)));
    
    % Avoid division by zero
    tmp2 = zeros(size(hL));
    tmp2(idx_validL) = hL(idx_validL) .* (bdL2(idx_validL).^3 - bgL2(idx_validL).^3) ./ denomL(idx_validL);
    
    % Right contribution
    sonmR = sqrt(g * hR / 2.0);
    uR = qR ./ hR;
    borneR = -sqrt(2.0 * g * neg(dZ));
    bgR = min(borneR, uR - sonmR .* sqrt(3.0));
    bdR = min(borneR, uR + sonmR .* sqrt(3.0));
    
    % Avoid division by zero
    denomR = 6.0 * sonmR .* sqrt(3.0);
    idx_validR = denomR > epsilon;
    tmp3 = zeros(size(hR));
    tmp3(idx_validR) = -hR(idx_validR) .* (pos(bdR(idx_validR).^2 + 2.0 * g * dZ(idx_validR)).^1.5 - ...
                                         pos(bgR(idx_validR).^2 + 2.0 * g * dZ(idx_validR)).^1.5) ./ ...
                                         denomR(idx_validR);
    
    % Combine
    f2m = tmp1 + tmp2 + tmp3;

end

function p = pos(x)

    % Compute the positive part of x
    p = max(0.0, x);

end

function n = neg(x)

    % Compute the negative part of x
    n = max(0.0, -x);

end