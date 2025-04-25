%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +flux/FORCE.m
%
% Purpose:
%   Calculates the numerical flux using the FORCE (First ORder CEntred)
%   scheme for the 1D Non-Linear Shallow Water Equations (NSW).
%   The FORCE scheme is an arithmetic average of the Lax-Friedrichs (LF)
%   and Richtmyer (RI) fluxes, aiming to balance dissipation and accuracy.
%
% Syntax:
%   F_num = FORCE(wL, wR, cfg)
%
% Inputs:
%   wL   - [N x 2, double] State vector [H, HU] on the left side of the interface.
%   wR   - [N x 2, double] State vector [H, HU] on the right side of the interface.
%   cfg  - [struct] Configuration structure. Required fields:
%            cfg.phys.g: [double] Acceleration due to gravity [m/s^2].
%            cfg.phys.dry_tolerance: [double] Tolerance for identifying dry states [m].
%            cfg.time.cfl: [double] CFL number (used to approximate dt/dx).
%
% Outputs:
%   F_num - [N x 2, double] Numerical flux vector [F_H, F_HU] across the interface.
%
% Dependencies:
%   None (standalone flux function, but expects correct cfg fields).
%
% References:
%   - Toro, E. F. (2009). Riemann Solvers and Numerical Methods for
%     Fluid Dynamics: A Practical Introduction (3rd ed.). Springer. (Chapter 6)
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F_num = FORCE(wL, wR, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Persistent Variables and Initialization                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use persistent variables for g, dry_tolerance and cfl to avoid repeated checks/access
    % within the loop calling this function. Initialize on first call.
    persistent g dry_tolerance cfl;

    if isempty(g) || isempty(dry_tolerance) || isempty(cfl)
        if isfield(cfg, 'phys') && isfield(cfg.phys, 'g')
            g = cfg.phys.g;     % [m/s^2] Gravity
        else
            g = 9.81; warning('FORCE:UsingDefault', 'Using default g = 9.81 m/s^2');
        end
        if isfield(cfg, 'phys') && isfield(cfg.phys, 'dry_tolerance')
            dry_tolerance = cfg.phys.dry_tolerance; % [m] Tolerance for identifying dry states
        else
            dry_tolerance = 1e-6; warning('FORCE:UsingDefault', 'Using default dry_tolerance = 1e-6 m');
        end
        if isfield(cfg, 'time') && isfield(cfg.time, 'cfl')
            cfl = cfg.time.cfl; % [unitless] CFL number
        else
            cfl = 0.9; warning('FORCE:UsingDefault', 'Using default CFL = 0.9');
        end
    end

    N = size(wL, 1); % Number of interfaces

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested **Vectorized** Physical Flux Function                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define F(w) = [HU; HU^2 + 0.5*g*H^2]
    function F = physical_flux_vec(w_mat) % Input N x 2 matrix
        h_loc = w_mat(:, 1); % N x 1
        q_loc = w_mat(:, 2); % N x 1
        
        F = zeros(N, 2); % Initialize N x 2 flux matrix
        F(:, 1) = q_loc; % Mass flux is always discharge [m^2/s]
        
        wet_idx = h_loc > dry_tolerance; % Logical index for wet cells
        u_loc = zeros(N, 1);
        u_loc(wet_idx) = q_loc(wet_idx) ./ h_loc(wet_idx); % [m/s] Velocity for wet cells
        
        % Momentum flux [m^3/s^2]
        F(wet_idx, 2) = q_loc(wet_idx) .* u_loc(wet_idx) + 0.5 * g * h_loc(wet_idx).^2;
        % Dry cells already have F(:, 2) = 0 from initialization
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Left/Right States and Maximum Wave Speed         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hL = wL(:, 1); qL = wL(:, 2); % N x 1
    hR = wR(:, 1); qR = wR(:, 2); % N x 1

    uL = zeros(N, 1); cL = zeros(N, 1);
    wetL = hL > dry_tolerance;
    uL(wetL) = qL(wetL) ./ hL(wetL);       % [m/s] Left velocity
    cL(wetL) = sqrt(g * hL(wetL));  % [m/s] Left wave speed

    uR = zeros(N, 1); cR = zeros(N, 1);
    wetR = hR > dry_tolerance;
    uR(wetR) = qR(wetR) ./ hR(wetR);       % [m/s] Right velocity
    cR(wetR) = sqrt(g * hR(wetR));  % [m/s] Right wave speed

    % Estimate maximum signal speed S_max = max(|u| + c) for Lax-Friedrichs [N x 1]
    S_max = max(abs(uL) + cL, abs(uR) + cR); % [m/s]
    S_max(S_max < cfg.numerics.epsilon) = cfg.numerics.epsilon; % Avoid division by zero if flow is stagnant

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Component Fluxes (LF and RI)                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate physical fluxes at left and right states (using vectorized function)
    F_L = physical_flux_vec(wL); % N x 2
    F_R = physical_flux_vec(wR); % N x 2

    % 1. Lax-Friedrichs (LF) Flux: [N x 2]
    %    F_LF = 0.5 * [F(wL) + F(wR) - S_max .* (wR - wL)]
    F_LF = 0.5 * (F_L + F_R - S_max .* (wR - wL)); % Broadcasting: (N x 1) .* (N x 2) -> (N x 2)

    % 2. Richtmyer (RI) Flux:
    %    F_RI = F(w_half), where w_half is an intermediate state.
    %    w_half = 0.5*(wL + wR) - 0.5*(dt/dx)*(F(wR) - F(wL))
    %    Approximate dt/dx using the CFL number: dt/dx approx cfl / S_max
    dtdx_approx = cfl ./ S_max; % [s/m] [N x 1]
    w_half = 0.5 * (wL + wR) - 0.5 * dtdx_approx .* (F_R - F_L); % [N x 2]. Broadcasting: (N x 1) .* (N x 2)
    F_RI = physical_flux_vec(w_half); % [N x 2]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate FORCE Flux (Average of LF and RI)                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % F_FORCE = 0.5 * (F_LF + F_RI)
    F_num = 0.5 * (F_LF + F_RI); % [N x 2]

end
