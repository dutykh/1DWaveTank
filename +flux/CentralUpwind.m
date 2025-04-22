%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/CentralUpwind.m
%
% Purpose:
%   Computes the Central-Upwind numerical flux for the 1D Non-Linear Shallow
%   Water Equations (NSW). This flux is based on the Kurganov-Tadmor/Kurganov-
%   Noelle-Petrova semi-discrete scheme that combines advantages of central
%   and upwind approaches, resulting in less numerical diffusion than first-order
%   methods while maintaining good stability properties.
%
% Syntax:
%   F = CentralUpwind(wL, wR, cfg)
%
% Inputs:
%   wL  - [1 x 2, double] State vector [H, HU] on the left side of the interface.
%   wR  - [1 x 2, double] State vector [H, HU] on the right side of the interface.
%   cfg - [struct] Configuration structure. Required fields: cfg.phys.g, cfg.phys.dry_tolerance.
%
% Outputs:
%   F   - [1 x 2, double] Central-Upwind numerical flux vector [F_H, F_HU] across the interface.
%
% Dependencies:
%   None (standalone flux function, but expects correct cfg.phys.g and cfg.phys.dry_tolerance).
%
% References:
%   - Kurganov, A., & Tadmor, E. (2000). New high-resolution central schemes 
%     for nonlinear conservation laws and convection-diffusion equations. 
%     Journal of Computational Physics, 160(1), 241-282.
%   - Kurganov, A., Noelle, S., & Petrova, G. (2001). Semidiscrete central-upwind 
%     schemes for hyperbolic conservation laws and Hamilton-Jacobi equations. 
%     SIAM Journal on Scientific Computing, 23(3), 707-740.
%   - Kurganov, A., & Petrova, G. (2007). A second-order well-balanced 
%     positivity preserving central-upwind scheme for the Saint-Venant system. 
%     Communications in Mathematical Sciences, 5(1), 133-160.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
% Date:   22 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = CentralUpwind(wL, wR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract Parameters and State Variables                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = cfg.phys.g;        % [m/s^2] Acceleration due to gravity
    eps_flux = cfg.phys.dry_tolerance;     % Tolerance for numerical stability & dry state

    % Ensure inputs are properly formatted
    if isvector(wL) && length(wL) == 2; wL = wL(:)'; end
    if isvector(wR) && length(wR) == 2; wR = wR(:)'; end

    % Extract states H and HU from left and right vectors
    hL = wL(:,1); huL = wL(:,2); % [m], [m^2/s]
    hR = wR(:,1); huR = wR(:,2); % [m], [m^2/s]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Primitive Variables and Wave Speeds               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Handle potential division by zero in dry states
    uL = zeros(size(hL)); idxL = hL > eps_flux; uL(idxL) = huL(idxL) ./ hL(idxL); % [m/s] Left velocity
    uR = zeros(size(hR)); idxR = hR > eps_flux; uR(idxR) = huR(idxR) ./ hR(idxR); % [m/s] Right velocity

    % Wave speeds (celerities)
    cL = sqrt(g * max(hL, 0)); % [m/s] Ensure non-negative depth for sqrt
    cR = sqrt(g * max(hR, 0)); % [m/s]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Physical Fluxes                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FL = core.utils.physical_flux(wL, cfg); % [m^2/s; m^3/s^2]
    FR = core.utils.physical_flux(wR, cfg); % [m^2/s; m^3/s^2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate One-Sided Local Speeds of Propagation             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute maximum and minimum wave speeds at the interface
    % These are the eigenvalues of the Jacobian: u ± sqrt(g*h)
    a_plus = max(uR + cR, uL + cL);  % Maximum right-going wave speed
    a_minus = min(uR - cR, uL - cL); % Minimum left-going wave speed
    
    % Ensure a_plus >= 0 and a_minus <= 0 for consistency
    a_plus = max(a_plus, 0);
    a_minus = min(a_minus, 0);
    
    % Handle the case where both speeds are 0 (avoid division by zero)
    idx_zero_speed = abs(a_plus - a_minus) < eps_flux;
    if any(idx_zero_speed)
        a_plus(idx_zero_speed) = eps_flux;
        a_minus(idx_zero_speed) = -eps_flux;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Central-Upwind Flux Formulation                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The central-upwind flux is based on the Kurganov-Noelle-Petrova (KNP) scheme.
    % F = (a_plus * FL - a_minus * FR) / (a_plus - a_minus) + 
    %     (a_plus * a_minus) / (a_plus - a_minus) * (wR - wL)
    
    % Calculate the denominator (avoid division by zero handled earlier)
    delta_a = a_plus - a_minus; % Nx1
    
    % Calculate the flux terms using element-wise operations for broadcasting
    term1 = (a_plus ./ delta_a) .* FL;   % (Nx1 ./ Nx1) .* Nx2 => Nx2
    term2 = (a_minus ./ delta_a) .* FR;  % (Nx1 ./ Nx1) .* Nx2 => Nx2
    term3 = ((a_plus .* a_minus) ./ delta_a) .* (wR - wL); % ((Nx1.*Nx1)./Nx1) .* Nx2 => Nx2
    
    % Combine terms for the KNP flux
    F = term1 - term2 + term3; % Nx2
    
    % Handle the purely upwind cases (where the formula might be less stable or simplifies)
    % When a_minus >= 0, all waves move rightward, use left flux
    idx_right = a_minus >= 0;
    F(idx_right, :) = FL(idx_right, :);
    
    % When a_plus <= 0, all waves move leftward
    idx_left = a_plus <= 0;
    F(idx_left, :) = FR(idx_left, :);
    
    % Handle any NaN values that might have occurred despite safeguards
    if any(isnan(F(:)))
        warning('CentralUpwind:NaNFlux', 'NaN detected in Central-Upwind flux calculation. Replacing with 0.');
        F(isnan(F)) = 0;
    end

end