%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/LF.m
%
% Purpose:
%   Computes the Lax-Friedrichs (LF) numerical flux for the 1D Non-Linear
%   Shallow Water Equations (NSW). The LF flux is a simple and robust scheme
%   that introduces numerical diffusion proportional to the maximum local
%   wave speed to achieve stability. It averages the physical fluxes and adds
%   a diffusive term based on the jump in the state vector.
%
% Syntax:
%   Phi = LF(vL, vR, cfg)
%
% Inputs:
%   vL  - [1 x 2, double] State vector [H, HU] on the left side of the interface.
%   vR  - [1 x 2, double] State vector [H, HU] on the right side of the interface.
%   cfg - [struct] Configuration structure. Required fields: cfg.phys.g, cfg.phys.dry_tolerance.
%
% Outputs:
%   Phi - [1 x 2, double] Lax-Friedrichs numerical flux vector [Phi_H, Phi_HU].
%
% Dependencies:
%   Requires +core/+utils/physical_flux.m function.
%   Expects correct cfg.phys.g and cfg.phys.dry_tolerance.
%
% References:
%   - LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems.
%     Cambridge University Press. (Chapter 6)
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Phi = LF(vL, vR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract Parameters and State Variables                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = cfg.phys.g;       % [m/s^2] Acceleration due to gravity
    eps_flux = cfg.phys.dry_tolerance;     % Tolerance for numerical stability & dry state

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
    % Estimate Maximum Wave Speed (alpha)                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The LF scheme adds diffusion proportional to the maximum characteristic
    % speed alpha = max(|u| + c) around the interface.
    cL = sqrt(g * Hl); % [m/s] Left wave celerity
    cR = sqrt(g * Hr); % [m/s] Right wave celerity
    alpha = max(abs(uL) + cL, abs(uR) + cR); % [m/s] Max speed estimate
    % Ensure alpha has the correct dimensions if vL/vR are matrices (should not happen here)
    if ~isscalar(alpha)
         alpha = max(alpha, [], 2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Physical Fluxes F(vL) and F(vR)                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FL = utils.physical_flux(vL, cfg); % [m^2/s; m^3/s^2]
    FR = utils.physical_flux(vR, cfg); % [m^2/s; m^3/s^2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Lax-Friedrichs Flux                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Formula: Phi_LF = 0.5 * (F(vL) + F(vR)) - 0.5 * alpha * (vR - vL)
    % It averages the physical fluxes and subtracts a diffusion term.
    Phi = 0.5 * (FL + FR) - 0.5 * alpha .* (vR - vL);

    % --- Explicitly zero out flux for fully dry interfaces ---
    is_dry_interface = (Hl <= eps_flux) & (Hr <= eps_flux);
    if any(is_dry_interface)
        Phi(is_dry_interface, :) = 0; % Set both components [Phi_H, Phi_HU] to 0
    end
    % --- End of modification ---

end