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
%   wL   - [1 x 2, double] State vector [H, HU] on the left side of the interface.
%   wR   - [1 x 2, double] State vector [H, HU] on the right side of the interface.
%   cfg  - [struct] Configuration structure. Required fields:
%            cfg.phys.g: [double] Acceleration due to gravity [m/s^2].
%            cfg.phys.dry_tolerance: [double] Tolerance for identifying dry states [m].
%            cfg.time.cfl: [double] CFL number (used to approximate dt/dx).
%
% Outputs:
%   F_num - [1 x 2, double] Numerical flux vector [F_H, F_HU] across the interface.
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested Physical Flux Function                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define F(w) = [HU; HU^2 + 0.5*g*H^2]
    function F = physical_flux(w_vec)
        h_loc = w_vec(1);
        q_loc = w_vec(2);
        if h_loc <= dry_tolerance
            F = [0; 0]; % No flux in dry state
        else
            u_loc = q_loc / h_loc; % [m/s] Velocity
            F = [q_loc; q_loc*u_loc + 0.5*g*h_loc^2]; % [m^2/s; m^3/s^2]
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Left/Right States and Maximum Wave Speed         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wL = wL(:); % Ensure column vector [HL; qL]
    wR = wR(:); % Ensure column vector [HR; qR]
    hL = wL(1); qL = wL(2);
    hR = wR(1); qR = wR(2);

    uL = 0; cL = 0;
    if hL > dry_tolerance
        uL = qL / hL;       % [m/s] Left velocity
        cL = sqrt(g * hL);  % [m/s] Left wave speed
    end

    uR = 0; cR = 0;
    if hR > dry_tolerance
        uR = qR / hR;       % [m/s] Right velocity
        cR = sqrt(g * hR);  % [m/s] Right wave speed
    end

    % Estimate maximum signal speed S_max = max(|u| + c) for Lax-Friedrichs
    S_max = max(abs(uL) + cL, abs(uR) + cR); % [m/s]
    if S_max < 1e-9
        S_max = 1e-9; % Avoid division by zero if flow is stagnant
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Component Fluxes (LF and RI)                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate physical fluxes at left and right states
    F_L = physical_flux(wL);
    F_R = physical_flux(wR);

    % 1. Lax-Friedrichs (LF) Flux:
    %    F_LF = 0.5 * [F(wL) + F(wR) - S_max*(wR - wL)]
    F_LF = 0.5 * (F_L + F_R - S_max * (wR - wL));

    % 2. Richtmyer (RI) Flux:
    %    F_RI = F(w_half), where w_half is an intermediate state predicted
    %    using a Lax-Wendroff type step.
    %    w_half = 0.5*(wL + wR) - 0.5*(dt/dx)*(F(wR) - F(wL))
    %    Approximate dt/dx using the CFL number: dt/dx approx cfl / S_max
    dtdx_approx = cfl / S_max; % [s/m]
    w_half = 0.5 * (wL + wR) - 0.5 * dtdx_approx * (F_R - F_L);
    F_RI = physical_flux(w_half);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate FORCE Flux (Average of LF and RI)                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % F_FORCE = 0.5 * (F_LF + F_RI)
    F_num = 0.5 * (F_LF + F_RI);
    F_num = F_num(:)'; % Ensure output is a row vector [F_H, F_HU]

end
