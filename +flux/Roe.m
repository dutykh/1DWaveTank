function Phi = Roe(vL, vR, cfg)

    % Roe numerical flux for the Nonlinear Shallow Water Equations.
    %   Phi = Roe(vL, vR, cfg) calculates the numerical flux across an
    %   interface using the Roe approximate Riemann solver.
    %
    %   Reference:
    %       Roe, P. L. (1981). Approximate Riemann solvers, parameter vectors,
    %       and difference schemes. Journal of Computational Physics, 43(2), 357-372.
    %       Toro, E. F. (2009). Riemann solvers and numerical methods for fluid dynamics:
    %       A practical introduction (3rd ed.). Springer. (Chapter 11)
    %       LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems.
    %       Cambridge University Press. (Chapter 15)
    %
    %   Note: This basic implementation does not include an entropy fix, which might
    %         be necessary to handle sonic points correctly and prevent expansion shocks.
    %
    %   Inputs:
    %       vL  - State vector [H; HU] at the left of the interface.
    %       vR  - State vector [H; HU] at the right of the interface.
    %       cfg - Configuration structure (must contain phys.g).
    %
    %   Outputs:
    %       Phi - Roe numerical flux vector [Phi_H; Phi_HU].

    % Extract gravity from config
    g = cfg.phys.g;

    % Define epsilon for numerical stability (especially for dry states)
    eps_flux = 1e-10;

    % Extract states H and HU from left and right vectors
    Hl = vL(:,1); HuL = vL(:,2);
    Hr = vR(:,1); HuR = vR(:,2);

    % Calculate primitive variable U (velocity) on left and right
    uL = zeros(size(Hl)); idxL = Hl > eps_flux; uL(idxL) = HuL(idxL) ./ Hl(idxL);
    uR = zeros(size(Hr)); idxR = Hr > eps_flux; uR(idxR) = HuR(idxR) ./ Hr(idxR);

    % --- Roe Averages ---
    sqrtHl = sqrt(Hl);
    sqrtHr = sqrt(Hr);
    H_tilde = 0.5 * (Hl + Hr); % Arithmetic mean H often used for NSW Roe
    u_tilde = (sqrtHl .* uL + sqrtHr .* uR) ./ (sqrtHl + sqrtHr + eps_flux); % Weighted average velocity
    c_tilde = sqrt(g * H_tilde); % Roe average wave speed

    % Handle cases where denominator is zero (Hl=Hr=0)
    idx_zero_H = (sqrtHl + sqrtHr) < eps_flux;
    u_tilde(idx_zero_H) = 0;
    c_tilde(idx_zero_H) = 0;

    % --- Eigenvalues of Roe Matrix ---
    lambda1_tilde = u_tilde - c_tilde;
    lambda2_tilde = u_tilde + c_tilde;

    % --- Right Eigenvectors of Roe Matrix (columns of R_tilde) ---
    % r1 = [1; u_tilde - c_tilde] = [1; lambda1_tilde]
    % r2 = [1; u_tilde + c_tilde] = [1; lambda2_tilde]

    % --- Wave Strengths (alpha = L_tilde * delta_v) ---
    delta_v = vR - vL; % Jump in conserved variables [dH; dHu]
    dH = delta_v(:,1);
    dHu = delta_v(:,2);

    % alpha = L_tilde * delta_v, where L = R^-1
    % alpha1 = ( (u_tilde + c_tilde)*dH - dHu ) / (2*c_tilde)
    % alpha2 = (-(u_tilde - c_tilde)*dH + dHu ) / (2*c_tilde)
    inv_2c = 0.5 ./ (c_tilde + eps_flux); % Avoid division by zero if c_tilde is zero
    alpha1 = inv_2c .* (lambda2_tilde .* dH - dHu);
    alpha2 = inv_2c .* (-lambda1_tilde .* dH + dHu);

    % --- Matrix Absolute Value times delta_v ---
    % |A_tilde|*delta_v = sum(|lambda_k| * alpha_k * r_k)
    %                  = |lambda1|*alpha1*r1 + |lambda2|*alpha2*r2
    abs_lambda1 = abs(lambda1_tilde);
    abs_lambda2 = abs(lambda2_tilde);

    term1_H  = abs_lambda1 .* alpha1;
    term1_Hu = abs_lambda1 .* alpha1 .* lambda1_tilde;
    term2_H  = abs_lambda2 .* alpha2;
    term2_Hu = abs_lambda2 .* alpha2 .* lambda2_tilde;

    absA_deltaV = [term1_H + term2_H, term1_Hu + term2_Hu];

    % --- Physical Fluxes ---
    FL = core.utils.physical_flux(vL, cfg);
    FR = core.utils.physical_flux(vR, cfg);

    % --- Roe Flux ---
    Phi = 0.5 * (FL + FR) - 0.5 * absA_deltaV;

end