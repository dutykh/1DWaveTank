function Phi = HLL(vL, vR, cfg)

% HLL Harten-Lax-van Leer numerical flux for Nonlinear Shallow Water Equations.
%   Phi = HLL(vL, vR, cfg) calculates the numerical flux across an
%   interface using the HLL approximate Riemann solver.
%
%   Reference:
%       Harten, A., Lax, P. D., & van Leer, B. (1983).
%       On upstream differencing and Godunov-type schemes for hyperbolic conservation laws.
%       SIAM Review, 25(1), 35-61.
%       Toro, E. F. (2009). Riemann solvers and numerical methods for fluid dynamics:
%       A practical introduction (3rd ed.). Springer. (Chapter 10)
%
%   Inputs:
%       vL  - State vector [H; HU] at the left of the interface.
%       vR  - State vector [H; HU] at the right of the interface.
%       cfg - Configuration structure (must contain phys.g).
%
%   Outputs:
%       Phi - HLL numerical flux vector [Phi_H; Phi_HU].

    % Extract gravity from config
    g = cfg.phys.g;

    % Define epsilon for numerical stability
    eps_flux = 1e-10;

    % Extract states H and HU from left and right vectors
    Hl = vL(:,1); HuL = vL(:,2);
    Hr = vR(:,1); HuR = vR(:,2);

    % Calculate primitive variable U (velocity) on left and right
    uL = zeros(size(Hl)); idxL = Hl > eps_flux; uL(idxL) = HuL(idxL) ./ Hl(idxL);
    uR = zeros(size(Hr)); idxR = Hr > eps_flux; uR(idxR) = HuR(idxR) ./ Hr(idxR);

    % Calculate wave speeds (celerity) on left and right
    cL = sqrt(g * Hl);
    cR = sqrt(g * Hr);

    % Estimate wave speeds SL and SR (Davis estimates)
    % Fastest left-going wave speed
    SL = min(uL - cL, uR - cR);
    % Fastest right-going wave speed
    SR = max(uL + cL, uR + cR);

    % Calculate physical fluxes on the left and right
    FL = core.utils.physical_flux(vL, cfg);
    FR = core.utils.physical_flux(vR, cfg);

    % HLL Flux Calculation
    % Initialize flux vector
    Phi = zeros(size(vL));

    % Case 1: SR <= 0 (All waves move left, flux is FR)
    idx1 = SR <= 0;
    Phi(idx1,:) = FR(idx1,:);

    % Case 2: SL >= 0 (All waves move right, flux is FL)
    idx2 = SL >= 0;
    Phi(idx2,:) = FL(idx2,:);

    % Case 3: SL < 0 < SR (Waves move apart, standard HLL flux)
    % We only need to calculate for indices where neither idx1 nor idx2 is true
    idx3 = ~idx1 & ~idx2;
    if any(idx3)
        SR_idx3 = SR(idx3);
        SL_idx3 = SL(idx3);
        FL_idx3 = FL(idx3,:);
        FR_idx3 = FR(idx3,:);
        vL_idx3 = vL(idx3,:);
        vR_idx3 = vR(idx3,:);
        
        denominator = SR_idx3 - SL_idx3;
        % Avoid division by zero if SR is very close to SL
        denominator(denominator < eps_flux) = eps_flux; 
        
        Phi(idx3,:) = (SR_idx3 .* FL_idx3 - SL_idx3 .* FR_idx3 + SR_idx3 .* SL_idx3 .* (vR_idx3 - vL_idx3)) ./ denominator;
    end

end
