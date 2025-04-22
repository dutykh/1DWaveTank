%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/PVM.m
%
% Purpose:
%   Computes the PVM (Polynomial Viscosity Matrix) numerical flux for the 1D
%   Non-Linear Shallow Water Equations (NSW). The PVM approach uses polynomial 
%   forms of the viscosity matrix, providing a flexible framework that can 
%   recover many classical schemes (Roe, Rusanov, etc.) and create new schemes 
%   with specific dissipation properties by changing the polynomial function.
%
% Syntax:
%   F = PVM(wL, wR, cfg)
%
% Inputs:
%   wL  - [Nx2, double] State vector [H, HU] on the left side of the interface.
%   wR  - [Nx2, double] State vector [H, HU] on the right side of the interface.
%   cfg - [struct] Configuration structure containing:
%           cfg.phys.g: [double] Acceleration due to gravity [m/s^2].
%           cfg.phys.dry_tolerance: [double] Tolerance for numerical stability & dry state.
%           cfg.pvm_degree: [integer, optional] Degree of polynomial (default: 1).
%                           1: Recovers Roe-type scheme
%                           0: Recovers Rusanov-type scheme (local Lax-Friedrichs)
%                           2: Higher accuracy with less dissipation
%
% Outputs:
%   F   - [Nx2, double] PVM numerical flux vector [F_H, F_HU] across the interface.
%
% Dependencies:
%   Requires +core/+utils/physical_flux.m function.
%   Expects correct cfg.phys.g and cfg.phys.dry_tolerance.
%
% References:
%   - Castro, M.J., Gallardo, J.M., Marquina, A. (2014). A class of incomplete
%     Riemann solvers based on uniform polynomial approximations to the absolute
%     value function. Journal of Scientific Computing, 61(1), 40-57.
%   - Castro, M.J., Fernández-Nieto, E.D. (2012). A class of computationally fast
%     first order finite volume solvers: PVM methods. SIAM Journal on Scientific
%     Computing, 34(4), A2173-A2196.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
% Date:   22 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = PVM(wL, wR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract Parameters and State Variables                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = cfg.phys.g;       % [m/s^2] Acceleration due to gravity
    eps_flux = cfg.phys.dry_tolerance;     % Tolerance for numerical stability & dry state
    
    % Determine polynomial degree for PVM method
    if isfield(cfg, 'pvm_degree')
        p_degree = cfg.pvm_degree;
    else
        p_degree = 1;     % Default: degree 1 (Roe-type)
    end
    
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
    % Calculate Roe Averages and Jacobian Eigenstructure          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Square roots for Roe averaging
    hL_sqrt = sqrt(max(hL, eps_flux));
    hR_sqrt = sqrt(max(hR, eps_flux));
    denom = hL_sqrt + hR_sqrt;
    
    % Roe-averaged velocity
    u_roe = (hL_sqrt .* uL + hR_sqrt .* uR) ./ denom;
    
    % Roe-averaged depth and wave speed
    h_roe = 0.5 * (hL + hR);
    c_roe = sqrt(g * h_roe);
    
    % Eigenvalues of the Jacobian matrix at the Roe state
    lambda1 = u_roe - c_roe;  % First eigenvalue (left-going wave)
    lambda2 = u_roe + c_roe;  % Second eigenvalue (right-going wave)
    
    % Maximum wave speed (used for scaling)
    lambda_max = max(abs([lambda1, lambda2]), [], 2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute Physical Fluxes                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FL = core.utils.physical_flux(wL, cfg);
    FR = core.utils.physical_flux(wR, cfg);
    
    % Average flux
    F_avg = 0.5 * (FL + FR);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct PVM Viscosity Matrix Based on Polynomial Degree   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize matrices for PVM calculation
    num_points = size(wL, 1);
    A_abs = zeros(num_points, 2, 2); % Viscosity matrix for each interface point
    
    % Right eigenvector matrix R and its inverse R^(-1)
    % For shallow water: R = [1, 1; lambda1, lambda2]
    % and R^(-1) = [0.5*(1+lambda2/c_roe), -0.5/c_roe; -0.5*(1+lambda1/c_roe), 0.5/c_roe]
    
    % State difference
    delta_w = wR - wL;
    
    for i = 1:num_points
        % Skip computation for identical states (zero jump)
        if all(abs(delta_w(i,:)) < eps_flux)
            continue;
        end
        
        % Local eigenvalues
        lam1 = lambda1(i);
        lam2 = lambda2(i);
        c = c_roe(i);
        
        % Normalized eigenvalues for polynomial
        if lambda_max(i) > eps_flux
            mu1 = lam1 / lambda_max(i);
            mu2 = lam2 / lambda_max(i);
        else
            mu1 = 0;
            mu2 = 0;
        end
        
        % Apply polynomial approximation to |λ|
        % P_0(λ) = 1 (Rusanov/LLF)
        % P_1(λ) = |λ| (Roe)
        % P_2(λ) = λ² for |λ| ≤ 1 (less dissipative)
        
        if p_degree == 0
            % Rusanov/Local Lax-Friedrichs: constant polynomial
            abs_mu1 = 1.0;
            abs_mu2 = 1.0;
        elseif p_degree == 1
            % Roe: linear polynomial (exact |λ|)
            abs_mu1 = abs(mu1);
            abs_mu2 = abs(mu2);
        elseif p_degree == 2
            % Quadratic polynomial: P_2(λ) = λ² for |λ| ≤ 1
            abs_mu1 = min(mu1^2, abs(mu1));
            abs_mu2 = min(mu2^2, abs(mu2));
        else
            % Default to Roe for unsupported degrees
            abs_mu1 = abs(mu1);
            abs_mu2 = abs(mu2);
        end
        
        % Scale back to original eigenvalue magnitudes
        abs_lam1 = abs_mu1 * lambda_max(i);
        abs_lam2 = abs_mu2 * lambda_max(i);
        
        % Construct diagonal matrix of polynomial-approximated |λ|
        Lambda_abs = diag([abs_lam1, abs_lam2]);
        
        % Right eigenvector matrix
        R = [1, 1; lam1, lam2];
        
        % Inverse of right eigenvector matrix
        % Using explicit formula for 2x2 case to avoid numerical issues
        det_R = lam2 - lam1; % Should be equal to 2*c
        if abs(det_R) < eps_flux
            det_R = eps_flux * sign(det_R); % Avoid division by zero
        end
        
        R_inv = [lam2, -1; -lam1, 1] / det_R;
        
        % Compute viscosity matrix |A| = R|Λ|R⁻¹
        A_abs(i,:,:) = R * Lambda_abs * R_inv;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute Final PVM Flux                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PVM flux formula: F_pvm = F_avg - 0.5 * |A| * (wR - wL)
    
    % Initialize flux
    F = zeros(size(F_avg));
    
    % Apply viscosity term to each interface point
    for i = 1:num_points
        % Extract viscosity matrix for this point
        A_abs_i = squeeze(A_abs(i,:,:));
        
        % Calculate viscosity term
        visc_term = 0.5 * A_abs_i * delta_w(i,:)';
        
        % Compute final flux
        F(i,:) = F_avg(i,:) - visc_term';
    end
    
    % Handle special cases where all waves move in one direction
    % When all waves move right, use left flux
    idx_right = lambda1 >= 0;
    F(idx_right,:) = FL(idx_right,:);
    
    % When all waves move left, use right flux
    idx_left = lambda2 <= 0;
    F(idx_left,:) = FR(idx_left,:);
    
    % Handle any NaN values
    if any(isnan(F(:)))
        warning('PVM:NaNFlux', 'NaN detected in PVM flux calculation. Replacing with 0.');
        F(isnan(F)) = 0;
    end

end