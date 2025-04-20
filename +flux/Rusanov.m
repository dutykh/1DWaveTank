function Phi = Rusanov(vL, vR, cfg)

    % Rusanov (Local Lax-Friedrichs) numerical flux for NSW Equations.
    %   Phi = Rusanov(vL, vR, cfg) calculates the numerical flux across an
    %   interface using the Rusanov scheme.
    %
    %   Reference:
    %       Toro, E. F. (2009). Riemann solvers and numerical methods for fluid dynamics:
    %       A practical introduction (3rd ed.). Springer. (Chapter 6.4)
    %       LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems.
    %       Cambridge University Press. (Chapter 6)
    %
    %   Note: This is essentially the same as the Lax-Friedrichs flux implemented
    %   in LF.m, as both use the maximum local wave speed for dissipation.
    %
    %   Inputs:
    %       vL  - State vector [H; HU] at the left of the interface.
    %       vR  - State vector [H; HU] at the right of the interface.
    %       cfg - Configuration structure (must contain phys.g).
    %
    %   Outputs:
    %       Phi - Rusanov numerical flux vector [Phi_H; Phi_HU].

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

    % Estimate maximum local wave speed (S_max) for the numerical viscosity term
    S_max = max(abs(uL) + cL, abs(uR) + cR);
    % Ensure S_max has the correct dimensions for element-wise operation
    S_max = max(S_max, [], 2); % Take row-wise max if vL/vR are matrices

    % Calculate physical fluxes on the left and right
    FL = core.utils.physical_flux(vL, cfg);
    FR = core.utils.physical_flux(vR, cfg);

    % Calculate the Rusanov flux
    Phi = 0.5 * (FL + FR) - 0.5 * S_max .* (vR - vL);

end