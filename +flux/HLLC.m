function Phi = HLLC(vL, vR, cfg)

    % HLLC Harten-Lax-van Leer-Contact numerical flux for NSW Equations.
    %   Phi = HLLC(vL, vR, cfg) calculates the numerical flux across an
    %   interface using the HLLC approximate Riemann solver.
    %
    %   Reference:
    %       Toro, E. F. (2009). Riemann solvers and numerical methods for fluid dynamics:
    %       A practical introduction (3rd ed.). Springer. (Chapter 10)
    %
    %       Batten, P., Clarke, N., Lambert, C., & Causon, D. M. (1997).
    %       On the choice of wave speeds for the HLLC Riemann solver.
    %       SIAM Journal on Scientific Computing, 18(6), 1553-1570.
    %
    %   Inputs:
    %       vL  - State vector [H; HU] at the left of the interface.
    %       vR  - State vector [H; HU] at the right of the interface.
    %       cfg - Configuration structure (must contain phys.g).
    %
    %   Outputs:
    %       Phi - HLLC numerical flux vector [Phi_H; Phi_HU].

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

    % Estimate wave speeds SL and SR (Davis estimates, as in HLL)
    SL = min(uL - cL, uR - cR);
    SR = max(uL + cL, uR + cR);

    % Calculate physical fluxes on the left and right
    FL = core.utils.physical_flux(vL, cfg);
    FR = core.utils.physical_flux(vR, cfg);

    % --- HLLC specific calculations ---

    % Estimate contact wave speed S_star (Toro, Eq. 10.60, simplified for NSW)
    % This uses pressures p = 0.5*g*H^2
    pL = 0.5 * g * Hl.^2;
    pR = 0.5 * g * Hr.^2;
    S_star = (pR - pL + HuL.*(SL - uL) - HuR.*(SR - uR)) ./ ...
             (Hl.*(SL - uL) - Hr.*(SR - uR) + eps_flux); % Add eps to avoid division by zero

    % Calculate intermediate states vL_star, vR_star (Toro, Eq. 10.38 / 10.58)
    % Note: v = [H; HU]
    factorL = Hl .* (SL - uL) ./ (SL - S_star + eps_flux);
    vL_star = zeros(size(vL));
    vL_star(:,1) = factorL; % H_L_star
    vL_star(:,2) = factorL .* S_star; % (HU)_L_star

    factorR = Hr .* (SR - uR) ./ (SR - S_star + eps_flux);
    vR_star = zeros(size(vR));
    vR_star(:,1) = factorR; % H_R_star
    vR_star(:,2) = factorR .* S_star; % (HU)_R_star

    % Calculate intermediate fluxes FL_star, FR_star (Toro, Eq. 10.37)
    FL_star = FL + SL .* (vL_star - vL);
    FR_star = FR + SR .* (vR_star - vR);

    % --- HLLC Flux Selection ---
    Phi = zeros(size(vL));

    % Region L: 0 <= SL
    idxL_region = SL >= 0;
    Phi(idxL_region,:) = FL(idxL_region,:);

    % Region R: SR <= 0
    idxR_region = SR <= 0;
    Phi(idxR_region,:) = FR(idxR_region,:);

    % Region L_star: SL < 0 <= S_star
    idxLstar_region = (SL < 0) & (S_star >= 0) & ~idxL_region & ~idxR_region;
    Phi(idxLstar_region,:) = FL_star(idxLstar_region,:);

    % Region R_star: S_star < 0 < SR
    idxRstar_region = (S_star < 0) & (SR > 0) & ~idxL_region & ~idxR_region;
    Phi(idxRstar_region,:) = FR_star(idxRstar_region,:);

end