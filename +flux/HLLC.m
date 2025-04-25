%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/HLLC.m
%
% Purpose:
%   Computes the HLLC (Harten-Lax-van Leer-Contact) numerical flux for the
%   1D Non-Linear Shallow Water Equations (NSW). The HLLC flux improves upon
%   HLL by explicitly incorporating the contact wave (and shear waves in
%   multi-D), leading to better resolution of contact discontinuities.
%   It assumes a three-wave structure (left, contact, right).
%
% Syntax:
%   Phi = HLLC(vL, vR, cfg)
%
% Inputs:
%   vL  - [Nx2, double] State vector [H, HU] on the left side of the interface.
%   vR  - [Nx2, double] State vector [H, HU] on the right side of the interface.
%   cfg - [struct] Configuration structure. Required fields: cfg.phys.g, cfg.phys.dry_tolerance.
%
% Outputs:
%   Phi - [Nx2, double] HLLC numerical flux vector [Phi_H, Phi_HU] across the interface.
%
% Dependencies:
%   Requires +core/+utils/physical_flux.m function.
%   Expects correct cfg.phys.g and cfg.phys.dry_tolerance.
%
% References:
%   - Toro, E. F. (2009). Riemann solvers and numerical methods for fluid dynamics:
%     A practical introduction (3rd ed.). Springer. (Chapter 10)
%   - Batten, P., Clarke, N., Lambert, C., & Causon, D. M. (1997).
%     On the choice of wave speeds for the HLLC Riemann solver.
%     SIAM Journal on Scientific Computing, 18(6), 1553-1570.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Phi = HLLC(vL, vR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract Parameters and State Variables                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = cfg.phys.g;       % [m/s^2] Acceleration due to gravity
    dry_tolerance = cfg.phys.dry_tolerance; % Physical threshold for dry states
    epsilon = cfg.numerics.epsilon;         % Numerical tolerance for stability
    N = size(vL, 1);      % Number of interfaces

    % Ensure inputs are Nx2
    if size(vL, 2) ~= 2 || size(vR, 2) ~= 2
        error('HLLC:InputDim', 'Input state vectors vL and vR must have 2 columns.');
    end

    % Work with 2xN column vectors internally for states
    vL_col = vL'; % Transpose to 2xN
    vR_col = vR'; % Transpose to 2xN

    % Extract states H and HU as 1xN row vectors
    Hl = vL_col(1,:); HuL = vL_col(2,:); % [m], [m^2/s]
    Hr = vR_col(1,:); HuR = vR_col(2,:); % [m], [m^2/s]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Primitive Variables (Velocity) - Result 1xN       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uL = zeros(1, N); idxL = Hl > dry_tolerance; uL(idxL) = HuL(idxL) ./ Hl(idxL); % [m/s] Left velocity
    uR = zeros(1, N); idxR = Hr > dry_tolerance; uR(idxR) = HuR(idxR) ./ Hr(idxR); % [m/s] Right velocity

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate Wave Speeds (Signal Velocities) SL and SR - Result 1xN %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ensure H is non-negative before sqrt
    Hl_safe = max(Hl, 0);
    Hr_safe = max(Hr, 0);
    cL = sqrt(g * Hl_safe); % [m/s] Left wave celerity 
    cR = sqrt(g * Hr_safe); % [m/s] Right wave celerity
    
    % Standard Davis estimates
    SL_davis = min(uL - cL, uR - cR); % [m/s] Min characteristic speed
    SR_davis = max(uL + cL, uR + cR); % [m/s] Max characteristic speed
    
    % Simple Entropy Fix (Ensures 0 speed is included if characteristic speeds span 0)
    SL = min(SL_davis, 0);
    SR = max(SR_davis, 0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Physical Fluxes F(vL) and F(vR) - Result 2xN      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assuming physical_flux takes Nx2 and returns Nx2, we transpose
    FL_col = utils.physical_flux(vL, cfg)'; % Transpose to 2xN
    FR_col = utils.physical_flux(vR, cfg)'; % Transpose to 2xN

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HLLC Specific Calculations: Contact Speed & Intermediate States %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % --- Estimate contact wave speed S_star --- Result 1xN
    pL = 0.5 * g * Hl.^2; % [Pa equiv] Pressure term left
    pR = 0.5 * g * Hr.^2; % [Pa equiv] Pressure term right
    denominator = (Hl.*(SL - uL) - Hr.*(SR - uR));
    denom_zero_idx = abs(denominator) < epsilon;
    denominator(denom_zero_idx) = epsilon .* sign(denominator(denom_zero_idx)); % Avoid division by zero
    denominator(denom_zero_idx & denominator == 0) = epsilon; % Handle exact zero case
    S_star = (pR - pL + HuL.*(SL - uL) - HuR.*(SR - uR)) ./ denominator; % [m/s]

    % --- Calculate intermediate ('star') states vL* and vR* --- Result 2xN
    denomL_star = SL - S_star;
    denomR_star = SR - S_star;
    denomL_zero_idx = abs(denomL_star) < epsilon;
    denomR_zero_idx = abs(denomR_star) < epsilon;
    denomL_star(denomL_zero_idx) = epsilon .* sign(denomL_star(denomL_zero_idx));
    denomL_star(denomL_zero_idx & denomL_star == 0) = epsilon;
    denomR_star(denomR_zero_idx) = epsilon .* sign(denomR_star(denomR_zero_idx));
    denomR_star(denomR_zero_idx & denomR_star == 0) = epsilon;

    factorL = Hl .* (SL - uL) ./ denomL_star; % 1xN
    factorR = Hr .* (SR - uR) ./ denomR_star; % 1xN

    vL_star = [factorL; factorL .* S_star]; % Creates 2xN
    vR_star = [factorR; factorR .* S_star]; % Creates 2xN

    % --- Calculate intermediate fluxes FL* and FR* --- Result 2xN
    % FL_star = FL + SL .* (vL_star - vL) -> Dimensions: (2xN) + (1xN) .* (2xN - 2xN)
    % Need to broadcast SL and SR for element-wise multiplication
    FL_star = FL_col + repmat(SL, 2, 1) .* (vL_star - vL_col);
    FR_star = FR_col + repmat(SR, 2, 1) .* (vR_star - vR_col);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HLLC Flux Selection based on Wave Speeds                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Select the appropriate flux based on the region defined by SL, S_star, SR.

    Phi_col = zeros(2, N); % Initialize internal flux as 2xN

    % Define logical indices for the regions (1xN)
    idxL_region = SL >= 0;
    idxR_region = SR <= 0;
    idxLstar_region = (SL < 0) & (S_star >= 0);
    idxRstar_region = (S_star < 0) & (SR > 0);

    % Assign fluxes based on regions (assigning columns)
    if any(idxL_region)
        Phi_col(:, idxL_region) = FL_col(:, idxL_region);
    end
    if any(idxLstar_region)
        Phi_col(:, idxLstar_region) = FL_star(:, idxLstar_region);
    end
    if any(idxRstar_region)
        Phi_col(:, idxRstar_region) = FR_star(:, idxRstar_region);
    end
    if any(idxR_region)
        Phi_col(:, idxR_region) = FR_col(:, idxR_region);
    end
    
    % Ensure no NaN values crept in (e.g. from division issues despite checks)
    if any(isnan(Phi_col(:)))
        warning('HLLC:NaNFlux', 'NaN detected in HLLC flux calculation.');
        % Consider adding debugging here, e.g., find(isnan(Phi_col))
        Phi_col(isnan(Phi_col)) = 0; % Replace NaN with 0 as a fallback?
    end

    % Transpose final flux back to Nx2 format for output
    Phi = Phi_col';

end