%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/Rusanov.m
%
% Purpose:
%   Computes the Rusanov (local Lax-Friedrichs) numerical flux for the 1D
%   Non-Linear Shallow Water Equations (NSW). The Rusanov flux is a simple,
%   robust, and very diffusive approximate Riemann solver that averages the
%   physical fluxes and adds a dissipation term proportional to the maximum
%   local wave speed at the interface.
%
% Syntax:
%   Phi = Rusanov(vL, vR, cfg)
%
% Inputs:
%   vL  - [1 x 2, double] State vector [H, HU] on the left side of the interface.
%   vR  - [1 x 2, double] State vector [H, HU] on the right side of the interface.
%   cfg - [struct] Configuration structure. Required field: cfg.phys.g (gravity).
%
% Outputs:
%   Phi - [1 x 2, double] Rusanov numerical flux vector [Phi_H, Phi_HU].
%
% Dependencies:
%   Requires +core/+utils/physical_flux.m function.
%   Expects correct cfg.phys.g.
%
% References:
%   - Toro, E. F. (2009). Riemann solvers and numerical methods for fluid dynamics:
%     A practical introduction (3rd ed.). Springer. (Chapter 6.4)
%   - LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems.
%     Cambridge University Press. (Chapter 6)
%
% Note:
%   This is essentially the same as the Lax-Friedrichs flux (LF.m), as both use
%   the maximum local wave speed for dissipation. Rusanov is often used as a
%   baseline for robustness.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Phi = Rusanov(vL, vR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract Parameters and State Variables                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = cfg.phys.g;       % [m/s^2] Acceleration due to gravity
    eps_flux = 1e-10;     % Tolerance for numerical stability & dry state

    % Ensure inputs are row vectors if they are vectors
    if isvector(vL); vL = vL(:)'; end
    if isvector(vR); vR = vR(:)'; end

    % Extract states H and HU from left and right vectors
    Hl = vL(:,1); HuL = vL(:,2); % [m], [m^2/s]
    Hr = vR(:,1); HuR = vR(:,2); % [m], [m^2/s]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Primitive Variables (Velocity)                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Handle potential division by zero in dry states
    uL = zeros(size(Hl)); idxL = Hl > eps_flux; uL(idxL) = HuL(idxL) ./ Hl(idxL); % [m/s] Left velocity
    uR = zeros(size(Hr)); idxR = Hr > eps_flux; uR(idxR) = HuR(idxR) ./ Hr(idxR); % [m/s] Right velocity

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate Maximum Local Wave Speed (S_max)                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The Rusanov flux adds diffusion proportional to the maximum characteristic
    % speed S_max = max(|u| + c) at the interface.
    cL = sqrt(g * Hl); % [m/s] Left wave celerity
    cR = sqrt(g * Hr); % [m/s] Right wave celerity
    S_max = max(abs(uL) + cL, abs(uR) + cR); % [m/s] Max speed estimate
    % Ensure S_max has the correct dimensions if vL/vR are matrices (should not happen here)
    if ~isscalar(S_max)
         S_max = max(S_max, [], 2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Physical Fluxes F(vL) and F(vR)                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FL = core.utils.physical_flux(vL, cfg); % [m^2/s; m^3/s^2]
    FR = core.utils.physical_flux(vR, cfg); % [m^2/s; m^3/s^2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Rusanov Flux                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Formula: Phi_Rusanov = 0.5 * (F(vL) + F(vR)) - 0.5 * S_max * (vR - vL)
    % It averages the physical fluxes and subtracts a diffusion term.
    Phi = 0.5 * (FL + FR) - 0.5 * S_max .* (vR - vL);

end