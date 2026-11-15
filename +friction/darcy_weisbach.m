%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +friction/darcy_weisbach.m
%
% Purpose:
%   Implements the Darcy-Weisbach friction law for shallow water flows. The
%   momentum equation of the depth-averaged flow over a wide channel reads
%       (HU)_t + (HU^2 + g H^2/2)_x = - g H z_x - g H S_f,
%   where the friction slope of the Darcy-Weisbach law is
%       S_f = f U|U| / (8 g H)                                  [-],
%   f being the dimensionless Darcy friction factor and H the hydraulic
%   radius of a wide channel (hydraulic diameter 4H). This function returns
%   the friction term of the momentum equation,
%       - g H S_f = - f U|U| / 8                                [m^2/s^2].
%   In uniform flow down a slope S_0, S_f = S_0 gives the normal velocity
%   U = (8 g H S_0 / f)^(1/2).
%
% Syntax:
%   friction_term = darcy_weisbach(H, HU, g, cfg)
%
% Inputs:
%   H            - [N x 1, double] Water depth at each cell [m].
%   HU           - [N x 1, double] Discharge at each cell [m^2/s].
%   g            - [scalar, double] Gravitational acceleration [m/s^2].
%   cfg          - [struct] Configuration structure containing:
%                  cfg.phys.darcy_f: [double] Darcy friction factor (constant value).
%                  cfg.phys.f_calculation: [string] Method to calculate f ('constant',
%                                         'colebrook_white').
%                  cfg.phys.dry_tolerance: [double] Threshold for dry cells [m].
%
% Outputs:
%   friction_term - [N x 1, double] Friction term - g H S_f of the momentum
%                   equation, added to d(HU)/dt [m^2/s^2] (zero in dry cells).
%
% References:
%   - Castro-Orgaz, O., Hager, W.H. (2019). Shallow Water Hydraulics.
%     Springer, Section 3.3.2, equation (3.10): S_f = f U^2 / (8 g R_h), with
%     R_h = H in a wide channel.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: 23 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function friction_term = darcy_weisbach(H, HU, g, cfg)
    % Initialize output vector
    friction_term = zeros(size(HU));

    % Get dry tolerance
    if isfield(cfg.phys, 'dry_tolerance')
        dry_tol = cfg.phys.dry_tolerance;
    else
        dry_tol = 1e-6; % Default
    end

    % Wet cells (H > dry_tol > 0, so that no division by zero occurs)
    wet = H > dry_tol;
    Hw = H(wet);
    U = HU(wet) ./ Hw;                     % [m/s] Velocity

    % Determine how to calculate friction factor f
    if isfield(cfg.phys, 'f_calculation') && strcmpi(cfg.phys.f_calculation, 'colebrook_white')
        % Colebrook-White formula, cell by cell
        f = friction.colebrook_white(Hw, U, cfg);
    else
        % Use constant f value (default)
        if isfield(cfg.phys, 'darcy_f') && ~isempty(cfg.phys.darcy_f)
            f = cfg.phys.darcy_f * ones(size(Hw));
        else
            warning('Darcy friction factor not specified. Using default f = 0.02');
            f = 0.02 * ones(size(Hw));
        end
    end

    % Friction slope of the Darcy-Weisbach law and friction term of the momentum equation
    Sf = f .* abs(U) .* U ./ (8 * g * Hw); % [-]
    friction_term(wet) = -g * Hw .* Sf;    % [m^2/s^2]
end
