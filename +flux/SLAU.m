%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/SLAU.m
%
% Purpose:
%   Computes the SLAU (Simple Low-dissipation AUSM) numerical flux for the
%   1D Non-Linear Shallow Water Equations (NSW). The SLAU flux is a variant
%   of the AUSM family with significantly reduced numerical dissipation,
%   making it well-suited for capturing vortical structures and maintaining
%   accuracy in both low and high Froude number regimes.
%
% Syntax:
%   F = SLAU(wL, wR, cfg)
%
% Inputs:
%   wL  - [1 x 2, double] State vector [H, HU] on the left side of the interface.
%   wR  - [1 x 2, double] State vector [H, HU] on the right side of the interface.
%   cfg - [struct] Configuration structure. Required fields: cfg.phys.g, cfg.phys.dry_tolerance.
%
% Outputs:
%   F   - [1 x 2, double] SLAU numerical flux vector [F_H, F_HU] across the interface.
%
% Dependencies:
%   None (standalone flux function, but expects correct cfg.phys.g and cfg.phys.dry_tolerance).
%
% References:
%   - Shima, E., & Kitamura, K. (2011). Parameter-Free Simple Low-Dissipation
%     AUSM-Family Scheme for All Speeds. AIAA Journal, 49(8), 1693-1709.
%   - Adapted for Shallow Water Equations by following the analogy between
%     compressible flow and shallow water equations.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
% Date:   22 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = SLAU(wL, wR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract Parameters and State Variables                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = cfg.phys.g;        % Acceleration due to gravity
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

    % Calculate wave speeds (celerity) for left and right states
    cL = sqrt(g * max(hL, 0)); % [m/s] Ensure non-negative depth for sqrt
    cR = sqrt(g * max(hR, 0)); % [m/s]

    % Calculate "Mach numbers" (Froude numbers for shallow water)
    ML = uL ./ (cL + eps_flux); % Avoid division by zero
    MR = uR ./ (cR + eps_flux);

    % Pressure terms (analogous to pressure in gas dynamics)
    pL = 0.5 * g * hL.^2; % [m^3/s^2]
    pR = 0.5 * g * hR.^2; % [m^3/s^2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SLAU-specific Parameters and Functions                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate average values
    h_avg = 0.5 * (hL + hR);                        % [m] Average depth
    c_avg = sqrt(g * max(h_avg, eps_flux));         % [m/s] Average wave speed
    u_avg = 0.5 * (uL + uR);                        % [m/s] Average velocity
    
    % Calculate local Froude number for shock detection
    Fr_local = abs(u_avg) ./ (c_avg + eps_flux);
    
    % SLAU blending function - decreases numerical dissipation in smooth regions
    % and increases it near discontinuities
    chi = min(1.0, Fr_local);
    
    % Pressure diffusion coefficient that adapts based on Froude number
    f_p = chi .* (1.0 - chi.^2);
    
    % Interface velocity for upwinding
    u_interface = 0.5 * (uL + uR) - 0.5 * (pR - pL) ./ (h_avg .* c_avg + eps_flux);
    
    % Modified SLAU mass flux calculation
    % This reduces dissipation compared to standard AUSM variants
    mdot = h_avg .* u_interface .* (1 - f_p);
    
    % Diffusion coefficient for interface velocity 
    alpha = 0.1875; % SLAU parameter (3/16)
    
    % Velocity diffusion term
    vel_diff = alpha * abs(pR - pL) ./ (h_avg .* c_avg.^2 + eps_flux);
    
    % Interface mass flux with dissipation control
    % The form ensures reduced dissipation in smooth flow regions
    mass_flux = mdot + vel_diff .* h_avg .* c_avg .* chi .* sign(u_interface);
    
    % Pressure blending for momentum flux - more accurate than simple average
    p_interface = 0.5 * (pL + pR) - 0.5 * f_p .* (pR - pL);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute SLAU Fluxes                                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate SLAU fluxes components
    % Mass flux (continuity equation)
    F_h = mass_flux;
    
    % Momentum flux with blended pressure term
    % Upwinded convective part plus pressure term
    F_hu = zeros(size(F_h));
    
    % Upwind the convective part based on interface mass flux
    pos_mdot = mass_flux > 0;
    F_hu(pos_mdot) = mass_flux(pos_mdot) .* uL(pos_mdot);
    F_hu(~pos_mdot) = mass_flux(~pos_mdot) .* uR(~pos_mdot);
    
    % Add pressure term to momentum flux
    F_hu = F_hu + p_interface;
    
    % Assemble final flux vector
    F = [F_h, F_hu];

    % Handle any NaN values that might have occurred despite safeguards
    if any(isnan(F(:)))
        warning('SLAU:NaNFlux', 'NaN detected in SLAU flux calculation. Replacing with 0.');
        F(isnan(F)) = 0;
    end

end