%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/HLLE.m
%
% Purpose:
%   Computes the HLLE (Harten-Lax-van Leer-Einfeldt) numerical flux for the
%   1D Non-Linear Shallow Water Equations (NSW). The HLLE flux is a variant
%   of the HLL approximate Riemann solver that uses Einfeldt's wave speed
%   estimates for improved accuracy and robustness, especially for shallow
%   water equations.
%
% Syntax:
%   F = HLLE(wL, wR, cfg)
%
% Inputs:
%   wL  - [Nx2, double] State vector [H, HU] on the left side of the interface.
%   wR  - [Nx2, double] State vector [H, HU] on the right side of the interface.
%   cfg - [struct] Configuration structure. Required fields: cfg.phys.g, cfg.phys.dry_tolerance.
%
% Outputs:
%   F   - [Nx2, double] HLLE numerical flux vector [F_H, F_HU] across the interface.
%
% Dependencies:
%   Requires +core/+utils/physical_flux.m function.
%   Expects correct cfg.phys.g and cfg.phys.dry_tolerance.
%
% References:
%   - Einfeldt, B. (1988). On Godunov-Type Methods for Gas Dynamics.
%     SIAM Journal on Numerical Analysis, 25(2), 294-318.
%   - Toro, E. F. (2009). Riemann solvers and numerical methods for fluid dynamics:
%     A practical introduction (3rd ed.). Springer. (Chapter 10)
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
% Date:   22 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = HLLE(wL, wR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract Parameters and State Variables                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = cfg.phys.g;        % [m/s^2] Acceleration due to gravity
    dry_tolerance = cfg.phys.dry_tolerance; % Physical threshold for dry states
    epsilon = cfg.numerics.epsilon;         % Numerical tolerance for stability

    % Ensure inputs are properly formatted
    if isvector(wL) && length(wL) == 2; wL = wL(:)'; end
    if isvector(wR) && length(wR) == 2; wR = wR(:)'; end

    % Extract states H and HU from left and right vectors
    Hl = wL(:,1); HuL = wL(:,2); % [m], [m^2/s]
    Hr = wR(:,1); HuR = wR(:,2); % [m], [m^2/s]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Primitive Variables (Velocity)                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Handle potential division by zero in dry states
    uL = zeros(size(Hl)); idxL = Hl > dry_tolerance; uL(idxL) = HuL(idxL) ./ Hl(idxL); % [m/s] Left velocity
    uR = zeros(size(Hr)); idxR = Hr > dry_tolerance; uR(idxR) = HuR(idxR) ./ Hr(idxR); % [m/s] Right velocity

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute Wave Speeds for Einfeldt's Estimates                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate wave speeds for both states
    cL = sqrt(g * max(Hl, 0)); % [m/s] Left wave celerity (ensure non-negative depth)
    cR = sqrt(g * max(Hr, 0)); % [m/s] Right wave celerity

    % Compute Roe averages for Einfeldt's wave speed estimates
    Hls = sqrt(max(Hl, 0)); % Square root of water depth (left)
    Hrs = sqrt(max(Hr, 0)); % Square root of water depth (right)
    
    % Handle potential division by zero
    denominator = Hls + Hrs;
    idx_zero_denom = denominator < epsilon;
    
    % Roe-averaged velocity
    u_roe = zeros(size(Hl));
    u_roe(~idx_zero_denom) = (Hls(~idx_zero_denom) .* uL(~idx_zero_denom) + ...
                            Hrs(~idx_zero_denom) .* uR(~idx_zero_denom)) ./ ...
                            denominator(~idx_zero_denom);
    
    % Roe-averaged depth and celerity
    h_roe = 0.5 * (Hl + Hr); % Arithmetic mean for water depth
    c_roe = sqrt(g * max(h_roe, 0)); % Celerity based on averaged depth

    % Einfeldt's estimates for the wave speeds - more accurate than standard HLL
    SL = min(uL - cL, uR - cR);
    SR = max(uL + cL, uR + cR);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Physical Fluxes F(wL) and F(wR)                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FL = utils.physical_flux(wL, cfg); % [m^2/s; m^3/s^2]
    FR = utils.physical_flux(wR, cfg); % [m^2/s; m^3/s^2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HLLE Flux Calculation based on Wave Speed Estimates         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize flux vector
    F = zeros(size(wL));

    % --- Case 1: All waves move left (SR <= 0) ---
    idx1 = SR <= 0;
    F(idx1,:) = FR(idx1,:);

    % --- Case 2: All waves move right (SL >= 0) ---
    idx2 = SL >= 0;
    F(idx2,:) = FL(idx2,:);

    % --- Case 3: Waves move apart (SL < 0 < SR) --- Standard HLLE formula ---
    % F_HLLE = (SR*FL - SL*FR + SL*SR*(wR - wL)) / (SR - SL)
    idx3 = ~idx1 & ~idx2;
    if any(idx3)
        SR_idx3 = SR(idx3);
        SL_idx3 = SL(idx3);
        FL_idx3 = FL(idx3,:);
        FR_idx3 = FR(idx3,:);
        wL_idx3 = wL(idx3,:);
        wR_idx3 = wR(idx3,:);
        
        denominator = SR_idx3 - SL_idx3;
        % Avoid division by zero if SR happens to equal SL (although unlikely with estimates used)
        denominator(denominator < epsilon) = epsilon; 
        
        F(idx3,:) = (SR_idx3 .* FL_idx3 - SL_idx3 .* FR_idx3 + SR_idx3 .* SL_idx3 .* (wR_idx3 - wL_idx3)) ./ denominator;
    end

    % Handle any NaN values that might have occurred despite safeguards
    if any(isnan(F(:)))
        warning('HLLE:NaNFlux', 'NaN detected in HLLE flux calculation. Replacing with 0.');
        F(isnan(F)) = 0;
    end

end