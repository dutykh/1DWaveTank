%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/FVCF.m
%
% Purpose:
%   Calculates the numerical flux using the FVCF (Flux Vector and
%   Characteristic Flux) scheme for the 1D Non-Linear Shallow Water Equations
%   (NSW). This scheme is a Roe-type method based on flux differences rather
%   than state differences for the dissipation term.
%
% Syntax:
%   Phi = FVCF(vL, vR, cfg)
%
% Inputs:
%   vL  - [Nx2, double] State vector [H, HU] on the left side of interfaces.
%   vR  - [Nx2, double] State vector [H, HU] on the right side of interfaces.
%   cfg - [struct] Configuration structure. Required fields: cfg.phys.g, cfg.phys.dry_tolerance.
%
% Outputs:
%   Phi - [Nx2, double] FVCF numerical flux vector [Phi_H, Phi_HU] across interfaces.
%
% Dependencies:
%   Requires +core/+utils/physical_flux.m function.
%   Expects correct cfg.phys.g and cfg.phys.dry_tolerance.
%
% References:
%   - Based on user-provided reference implementation.
%   - Ghidaglia, J.-M., Kumbaro, A., & Le Coq, G. (2001).
%     On the numerical solution to two fluid models via a cell centered finite volume method.
%     European Journal of Mechanics - B/Fluids, 20(6), 841-867. (Likely related to naming)
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Phi = FVCF(vL, vR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract Parameters and State Variables                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = cfg.phys.g;        % [m/s^2] Acceleration due to gravity
    eps_flux = cfg.phys.dry_tolerance;      % Tolerance for numerical stability & dry state handling
    N = size(vL, 1);       % Number of interfaces

    % Ensure inputs are Nx2
    if size(vL, 2) ~= 2 || size(vR, 2) ~= 2
        error('FVCF:InputDim', 'Input state vectors vL and vR must have 2 columns.');
    end

    % Extract states H and HU as Nx1 column vectors
    Hl = vL(:,1); HuL = vL(:,2); % [m], [m^2/s]
    Hr = vR(:,1); HuR = vR(:,2); % [m], [m^2/s]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Primitive Variables (Velocity) - Result Nx1       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uL = zeros(N, 1); idxL = Hl > eps_flux; uL(idxL) = HuL(idxL) ./ Hl(idxL); % [m/s] Left velocity
    uR = zeros(N, 1); idxR = Hr > eps_flux; uR(idxR) = HuR(idxR) ./ Hr(idxR); % [m/s] Right velocity

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Roe Weights and Averages - Result Nx1             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Hls = sqrt(max(Hl, 0)); % Weights for Roe averaging (ensure non-negative H)
    Hrs = sqrt(max(Hr, 0)); % Weights for Roe averaging

    % Averaged states on faces (using reference variable names)
    mu1 = 0.5 * (Hl + Hr); % [m] Averaged depth (arithmetic mean)
    % Roe averaged velocity
    denominator_u = Hls + Hrs;
    idx_zero_denom_u = denominator_u < eps_flux;
    mu2 = zeros(N, 1);
    mu2(~idx_zero_denom_u) = (Hls(~idx_zero_denom_u) .* uL(~idx_zero_denom_u) + Hrs(~idx_zero_denom_u) .* uR(~idx_zero_denom_u)) ./ denominator_u(~idx_zero_denom_u); 
    % Handle case where both H are near zero: set velocity to 0
    % mu2(idx_zero_denom_u) = 0; % Already initialized to 0
    
    % Averaged wave speed
    cm = sqrt(g * max(mu1, 0)); % [m/s]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Eigenvalue Sign Coefficients - Result Nx1         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda1 = mu2 + cm;
    lambda2 = mu2 - cm;
    s1 = sign(lambda1);
    s2 = sign(lambda2);

    denominator_sd = cm;
    idx_zero_denom_sd = denominator_sd < eps_flux;
    sd = zeros(N, 1);
    sd(~idx_zero_denom_sd) = 0.5 * (s2(~idx_zero_denom_sd) - s1(~idx_zero_denom_sd)) ./ denominator_sd(~idx_zero_denom_sd);
    % Handle case where cm is near zero (implies H is near zero): sd is ill-defined, set to 0?
    sd(idx_zero_denom_sd) = 0; 
    
    sp = 0.5 * (s1 + s2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Physical Fluxes and Differences - Result Nx2     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FL = core.utils.physical_flux(vL, cfg); % Nx2
    FR = core.utils.physical_flux(vR, cfg); % Nx2

    Phi_avg = 0.5 * (FL + FR); % Nx2 - Average flux
    fd = 0.5 * (FL - FR);      % Nx2 - Flux difference (note sign convention from ref)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply FVCF Correction                                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize flux with the average flux
    Phi = Phi_avg;

    % Apply corrections based on reference formula (vectorized)
    % Phi(:,1) = Phi(:,1) + (sd.*mu2 + sp).*fd(:,1) - sd.*fd(:,2);
    % Phi(:,2) = Phi(:,2) + sd.*(mu2.^2 - cm.^2).*fd(:,1) + (sp - sd.*mu2).*fd(:,2);

    % Apply corrections column-wise for clarity
    correction1 = (sd .* mu2 + sp) .* fd(:,1) - sd .* fd(:,2);
    correction2 = sd .* (mu2.^2 - cm.^2) .* fd(:,1) + (sp - sd .* mu2) .* fd(:,2);
    
    Phi(:,1) = Phi(:,1) + correction1;
    Phi(:,2) = Phi(:,2) + correction2;
    
    % Ensure no NaN values crept in
    if any(isnan(Phi(:)))
        warning('FVCF:NaNFlux', 'NaN detected in FVCF flux calculation.');
        Phi(isnan(Phi)) = 0; % Replace NaN with 0 as a fallback
    end

end % FVCF