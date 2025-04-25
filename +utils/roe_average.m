function [H_roe, U_roe, C_roe] = roe_average(wL, wR, cfg)
% ROE_AVERAGE Calculates the Roe average state for the shallow water equations.
%
% Syntax:
%   [H_roe, U_roe, C_roe] = roe_average(wL, wR, cfg)
%
% Inputs:
%   wL  - [1 x 2] Left state vector [H_L, HU_L]
%   wR  - [1 x 2] Right state vector [H_R, HU_R]
%   cfg - [struct] Configuration structure containing cfg.phys.g and cfg.phys.dry_tolerance
%
% Outputs:
%   H_roe - [scalar] Roe averaged depth
%   U_roe - [scalar] Roe averaged velocity
%   C_roe - [scalar] Roe averaged wave speed (sqrt(g*H_roe))
%
% Reference:
%   Toro, E. F. (2009). Riemann Solvers and Numerical Methods for Fluid Dynamics.
%   Springer. (Chapter 10)

    g = cfg.phys.g;
    dry_tol = cfg.phys.dry_tolerance;

    HL = wL(1);
    HR = wR(1);
    HUL = wL(2);
    HUR = wR(2);

    % Handle dry states to avoid division by zero
    UL = 0; UR = 0;
    if HL > dry_tol
        UL = HUL / HL;
    end
    if HR > dry_tol
        UR = HUR / HR;
    end

    % Calculate Roe averages
    sqrt_HL = sqrt(HL);
    sqrt_HR = sqrt(HR);
    
    % Simple average for H_roe (alternatives exist, but this is common)
    H_roe = 0.5 * (HL + HR); 
    
    % Roe average velocity U_roe
    if (sqrt_HL + sqrt_HR) < dry_tol % Both essentially dry
        U_roe = 0.0;
    else
        U_roe = (sqrt_HL * UL + sqrt_HR * UR) / (sqrt_HL + sqrt_HR);
    end

    % Roe average wave speed C_roe
    if H_roe < dry_tol
        C_roe = 0.0;
    else
        C_roe = sqrt(g * H_roe);
    end
    
end
