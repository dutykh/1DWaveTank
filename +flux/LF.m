function Phi = LF(vL, vR, cfg)

% LF Lax-Friedrichs numerical flux for Nonlinear Shallow Water Equations.
%   Phi = LF(vL, vR, cfg) calculates the numerical flux across an
%   interface using the Lax-Friedrichs scheme.
%
%   Reference:
%       LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems.
%       Cambridge University Press. (Chapter 6)
%
%   Inputs:
%       vL  - State vector [H; HU] at the left of the interface.
%       vR  - State vector [H; HU] at the right of the interface.
%       cfg - Configuration structure (must contain phys.g and mesh.dx).
%
%   Outputs:
%       Phi - Lax-Friedrichs numerical flux vector [Phi_H; Phi_HU].

    % Extract gravity from config
    g = cfg.phys.g;

    % Define epsilon for numerical stability
    eps_flux = 1e-10;

    % Extract states H and HU from left and right vectors
    Hl = vL(:,1); HuL = vL(:,2);
    Hr = vR(:,1); HuR = vR(:,2);

    % Calculate primitive variable U (velocity) on left and right
    % Avoid division by zero for dry states
    uL = zeros(size(Hl)); idxL = Hl > eps_flux; uL(idxL) = HuL(idxL) ./ Hl(idxL);
    uR = zeros(size(Hr)); idxR = Hr > eps_flux; uR(idxR) = HuR(idxR) ./ Hr(idxR);

    % Calculate wave speeds (celerity) on left and right
    cL = sqrt(g * Hl);
    cR = sqrt(g * Hr);

    % Estimate maximum wave speed (alpha) for the numerical viscosity term
    alpha = max(abs(uL) + cL, abs(uR) + cR);
    % Ensure alpha has the correct dimensions for element-wise operation
    alpha = max(alpha, [], 2); % Take row-wise max if vL/vR are matrices

    % Calculate physical fluxes on the left and right
    FL = core.utils.physical_flux(vL, cfg);
    FR = core.utils.physical_flux(vR, cfg);

    % Calculate the Lax-Friedrichs flux
    Phi = 0.5 * (FL + FR) - 0.5 * alpha .* (vR - vL);

end