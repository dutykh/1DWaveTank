%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/Roe.m
%
% Purpose:
%   Computes the Roe numerical flux for the 1D Non-Linear Shallow Water
%   Equations (NSW). The Roe flux is a highly regarded approximate Riemann
%   solver that uses Roe-averaged states to construct a linearized problem,
%   providing sharp resolution for shocks and contacts when implemented correctly.
%
% Syntax:
%   Phi = Roe(vL, vR, cfg)
%
% Inputs:
%   vL  - [1 x 2, double] State vector [H, HU] on the left side of the interface.
%   vR  - [1 x 2, double] State vector [H, HU] on the right side of the interface.
%   cfg - [struct] Configuration structure. Required fields: cfg.phys.g, cfg.phys.dry_tolerance.
%
% Outputs:
%   Phi - [1 x 2, double] Roe numerical flux vector [Phi_H, Phi_HU].
%
% Dependencies:
%   Requires +core/+utils/physical_flux.m function.
%   Expects correct cfg.phys.g and cfg.phys.dry_tolerance.
%
% References:
%   - Roe, P. L. (1981). Approximate Riemann solvers, parameter vectors,
%     and difference schemes. Journal of Computational Physics, 43(2), 357-372.
%   - Toro, E. F. (2009). Riemann solvers and numerical methods for fluid dynamics:
%     A practical introduction (3rd ed.). Springer. (Chapter 11)
%   - LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems.
%     Cambridge University Press. (Chapter 15)
%
% Note:
%   This basic implementation does not include an entropy fix (e.g., Harten's
%   entropy fix), which might be necessary to handle sonic points correctly
%   and prevent expansion shocks in certain scenarios.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Phi = Roe(vL, vR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract Parameters and State Variables                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = cfg.phys.g;        % [m/s^2] Acceleration due to gravity
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
    % Calculate Roe Averages                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Roe averages define a state (u_tilde, c_tilde) such that the jump
    % condition delta(F) = A_roe(v_tilde) * delta(v) is satisfied exactly.
    sqrtHl = sqrt(Hl); % Square root weights for averaging
    sqrtHr = sqrt(Hr);

    % Roe average velocity
    u_tilde = (sqrtHl .* uL + sqrtHr .* uR) ./ (sqrtHl + sqrtHr + eps_flux); % [m/s]

    % Roe average celerity - often simplified for NSW
    % Note: H_tilde = 0.5*(Hl+Hr) is a common simplification for NSW Roe average
    H_tilde = 0.5 * (Hl + Hr);
    c_tilde = sqrt(g * H_tilde); % [m/s]

    % Handle cases where denominator sqrtHl + sqrtHr is zero (Hl=Hr=0)
    idx_zero_H = (sqrtHl + sqrtHr) < eps_flux;
    u_tilde(idx_zero_H) = 0;
    c_tilde(idx_zero_H) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Eigenvalues and Wave Strengths of Roe Matrix               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Eigenvalues of the Roe matrix A_roe(v_tilde)
    lambda1_tilde = u_tilde - c_tilde; % [m/s]
    lambda2_tilde = u_tilde + c_tilde; % [m/s]

    % Calculate Wave Strengths (alpha = L_tilde * delta_v)
    % These represent the projection of the jump delta_v onto the eigenvectors.
    delta_v = vR - vL; % Jump in conserved variables [dH; dHu]
    dH = delta_v(:,1);   % [m]
    dHu = delta_v(:,2);  % [m^2/s]

    % Formulas for alpha1, alpha2 derived from L_tilde * delta_v
    inv_2c = 0.5 ./ (c_tilde + eps_flux); % Avoid division by zero if c_tilde is zero
    alpha1 = inv_2c .* (lambda2_tilde .* dH - dHu); % Strength of wave 1
    alpha2 = inv_2c .* (-lambda1_tilde .* dH + dHu); % Strength of wave 2

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Roe Dissipation Term                              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The Roe flux dissipation term is 0.5 * |A_tilde| * delta_v
    % |A_tilde| * delta_v = sum_k |lambda_k_tilde| * alpha_k * r_k_tilde
    % where r_k_tilde are the right eigenvectors: r1=[1; lam1], r2=[1; lam2]

    abs_lambda1 = abs(lambda1_tilde);
    abs_lambda2 = abs(lambda2_tilde);

    % Component-wise calculation of the dissipation term vector
    dissipation_H  = abs_lambda1 .* alpha1 .* 1 + abs_lambda2 .* alpha2 .* 1;
    dissipation_Hu = abs_lambda1 .* alpha1 .* lambda1_tilde + abs_lambda2 .* alpha2 .* lambda2_tilde;

    absA_deltaV = [dissipation_H, dissipation_Hu]; % [m; m^2/s]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Physical Fluxes F(vL) and F(vR)                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FL = utils.physical_flux(vL, cfg); % [m^2/s; m^3/s^2]
    FR = utils.physical_flux(vR, cfg); % [m^2/s; m^3/s^2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Roe Flux                                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Formula: Phi_Roe = 0.5 * (F(vL) + F(vR)) - 0.5 * |A_tilde| * delta_v
    Phi = 0.5 * (FL + FR) - 0.5 * absA_deltaV;

end