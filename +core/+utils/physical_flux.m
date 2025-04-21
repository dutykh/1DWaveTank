%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +core/+utils/physical_flux.m
%
% Purpose:
%   Calculates the physical flux vector F(w) for the 1D Non-Linear Shallow
%   Water Equations (NSW). The NSW equations are given by:
%       dH/dt + d(HU)/dx = 0
%       d(HU)/dt + d(HU^2 + 0.5*g*H^2)/dx = -g*H*dh/dx + S_f
%   This function calculates the flux terms appearing in the spatial derivatives:
%       F(w) = [F_H; F_HU] = [HU; HU^2 + 0.5*g*H^2]
%
% Syntax:
%   F = physical_flux(w, cfg)
%
% Inputs:
%   w   - [N x 2, double] State vector array for N cells, where w(:,1) = H
%         (water depth) [m] and w(:,2) = HU (discharge) [m^2/s].
%   cfg - [struct] Configuration structure. Required fields:
%         cfg.phys.g: [double] Acceleration due to gravity [m/s^2].
%
% Outputs:
%   F   - [N x 2, double] Physical flux vector array, where F(:,1) = F_H [m^2/s]
%         and F(:,2) = F_HU [m^3/s^2].
%
% Dependencies:
%   None (utility function, but expects correct cfg structure).
%
% References:
%   - Vreugdenhil, C. B. (1994). Numerical Methods for Shallow-Water Flow.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = physical_flux(w, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract parameters and state variables                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = cfg.phys.g;     % [m/s^2] Acceleration due to gravity

    H = w(:, 1);        % [m] Water depth
    HU = w(:, 2);       % [m^2/s] Discharge (Momentum)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate primitive variable U (velocity)                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Need velocity U = HU/H to compute the momentum flux component.
    % Handle potential division by zero in dry cells (H -> 0).
    U = zeros(size(H)); % [m/s] Initialize velocity
    wet_indices = H > 1e-10; % Indices of cells with non-negligible depth
    U(wet_indices) = HU(wet_indices) ./ H(wet_indices);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate physical flux components                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % F_H: Flux for the continuity equation (mass conservation)
    F_H = HU; % [m^2/s]

    % F_HU: Flux for the momentum equation
    % Includes convective term (rho * u^2 * A -> HU*U) and pressure term (p*A -> 0.5*g*H^2)
    F_HU = HU .* U + 0.5 * g * H.^2; % [m^3/s^2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Combine into flux vector F = [F_H, F_HU]                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F = [F_H, F_HU];

end