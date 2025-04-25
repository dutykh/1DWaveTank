%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/HLL.m
%
% Purpose:
%   Computes the HLL (Harten-Lax-van Leer) numerical flux for the 1D
%   Non-Linear Shallow Water Equations (NSW). The HLL flux is an approximate
%   Riemann solver that considers a two-wave structure (left- and right-going
%   waves) separated by a contact discontinuity. It is known for its
%   robustness, especially for strong shocks.
%
% Syntax:
%   Phi = HLL(vL, vR, cfg)
%
% Inputs:
%   vL  - [1 x 2, double] State vector [H, HU] on the left side of the interface.
%   vR  - [1 x 2, double] State vector [H, HU] on the right side of the interface.
%   cfg - [struct] Configuration structure. Required fields: cfg.phys.g, cfg.phys.dry_tolerance.
%
% Outputs:
%   Phi - [1 x 2, double] HLL numerical flux vector [Phi_H, Phi_HU] across the interface.
%
% Dependencies:
%   None (standalone flux function, but expects correct cfg.phys.g and cfg.phys.dry_tolerance).
%
% References:
%   - Harten, A., Lax, P. D., & van Leer, B. (1983).
%     On upstream differencing and Godunov-type schemes for hyperbolic conservation laws.
%     SIAM Review, 25(1), 35-61.
%   - Toro, E. F. (2009). Riemann solvers and numerical methods for fluid dynamics:
%     A practical introduction (3rd ed.). Springer. (Chapter 10)
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Phi = HLL(vL, vR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract Parameters and State Variables                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = cfg.phys.g;        % [m/s^2] Acceleration due to gravity
    dry_tolerance = cfg.phys.dry_tolerance; % Physical threshold for dry states
    epsilon = cfg.numerics.epsilon;         % Numerical tolerance for stability

    % Extract states H and HU from left and right vectors
    Hl = vL(:,1); HuL = vL(:,2); % [m], [m^2/s]
    Hr = vR(:,1); HuR = vR(:,2); % [m], [m^2/s]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Primitive Variables (Velocity)                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Handle potential division by zero in dry states
    uL = zeros(size(Hl)); idxL = Hl > dry_tolerance; uL(idxL) = HuL(idxL) ./ Hl(idxL); % [m/s] Left velocity
    uR = zeros(size(Hr)); idxR = Hr > dry_tolerance; uR(idxR) = HuR(idxR) ./ Hr(idxR); % [m/s] Right velocity

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate Wave Speeds (Signal Velocities) SL and SR          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HLL requires estimates for the fastest left-going (SL) and
    % right-going (SR) wave speeds emerging from the Riemann problem
    % at the interface. Various estimates exist; Davis estimates are common.
    cL = sqrt(g * Hl); % [m/s] Left wave celerity
    cR = sqrt(g * Hr); % [m/s] Right wave celerity

    % Davis estimates (or similar conservative estimates)
    SL = min(uL - cL, uR - cR); % [m/s] Min characteristic speed (fastest left)
    SR = max(0, max(uL + cL, uR + cR)); % [m/s] Max characteristic speed (fastest right)
    % Other estimates (e.g., based on Roe averages) could also be used.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Physical Fluxes F(vL) and F(vR)                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FL = utils.physical_flux(vL, cfg); % [m^2/s; m^3/s^2]
    FR = utils.physical_flux(vR, cfg); % [m^2/s; m^3/s^2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HLL Flux Calculation based on Wave Speed Estimates         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The HLL flux depends on the estimated signal speeds SL and SR relative to zero.
    % Formula: See Toro (2009), Eq. 10.28.
    % Phi_HLL = ( SR*FL - SL*FR + SL*SR*(vR - vL) ) / ( SR - SL )  if SL < 0 < SR
    % Phi_HLL = FL                                                if 0 <= SL
    % Phi_HLL = FR                                                if SR <= 0

    % Initialize flux vector
    Phi = zeros(size(vL));

    % --- Case 1: All waves move left (SR <= 0) ---
    idx1 = SR <= 0;
    Phi(idx1,:) = FR(idx1,:);

    % --- Case 2: All waves move right (SL >= 0) ---
    idx2 = SL >= 0;
    Phi(idx2,:) = FL(idx2,:);

    % --- Case 3: Waves move apart (SL < 0 < SR) --- Standard HLL formula ---
    idx3 = ~idx1 & ~idx2;
    if any(idx3)
        SR_idx3 = SR(idx3);
        SL_idx3 = SL(idx3);
        FL_idx3 = FL(idx3,:);
        FR_idx3 = FR(idx3,:);
        vL_idx3 = vL(idx3,:);
        vR_idx3 = vR(idx3,:);
        
        denominator = SR_idx3 - SL_idx3;
        % Avoid division by zero if SR happens to equal SL (although unlikely with estimates used)
        denominator(denominator < epsilon) = epsilon; 
        
        Phi(idx3,:) = (SR_idx3 .* FL_idx3 - SL_idx3 .* FR_idx3 + SR_idx3 .* SL_idx3 .* (vR_idx3 - vL_idx3)) ./ denominator;
    end

end