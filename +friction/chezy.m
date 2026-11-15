%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +friction/chezy.m
%
% Purpose:
%   Implements the Chézy friction law for shallow water flows. The momentum
%   equation of the depth-averaged flow over a wide channel reads
%       (HU)_t + (HU^2 + g H^2/2)_x = - g H z_x - g H S_f,
%   where the friction slope of Chézy's law is
%       S_f = U|U| / (C^2 H)                                    [-],
%   C being the Chézy coefficient [m^(1/2)/s] and H the hydraulic radius of a
%   wide channel. This function returns the friction term of the momentum
%   equation,
%       - g H S_f = - g U|U| / C^2                              [m^2/s^2].
%   In uniform flow down a slope S_0, S_f = S_0 gives the normal velocity
%   U = C (H S_0)^(1/2).
%
% Syntax:
%   friction_term = chezy(H, HU, g, cfg)
%
% Inputs:
%   H            - [N x 1, double] Water depth at each cell [m].
%   HU           - [N x 1, double] Discharge at each cell [m^2/s].
%   g            - [scalar, double] Gravitational acceleration [m/s^2].
%   cfg          - [struct] Configuration structure containing:
%                  cfg.phys.chezy_C: [double] Chézy coefficient [m^(1/2)/s].
%                  cfg.phys.dry_tolerance: [double] Threshold for dry cells [m].
%
% Outputs:
%   friction_term - [N x 1, double] Friction term - g H S_f of the momentum
%                   equation, added to d(HU)/dt [m^2/s^2] (zero in dry cells).
%
% References:
%   - Castro-Orgaz, O., Hager, W.H. (2019). Shallow Water Hydraulics.
%     Springer, Section 3.3.2: the Chézy coefficient C = (8g/f)^(1/2) and the
%     Darcy-Weisbach friction slope S_f = f U^2/(8 g R_h), equation (3.10),
%     whence S_f = U^2/(C^2 R_h), with R_h = H in a wide channel.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: 23 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function friction_term = chezy(H, HU, g, cfg)

    % Initialize output vector
    friction_term = zeros(size(HU));

    % Get Chézy coefficient from config
    if isfield(cfg.phys, 'chezy_C') && ~isempty(cfg.phys.chezy_C)
        C = cfg.phys.chezy_C;
    else
        warning('Chézy coefficient not specified. Using default C = 50 m^(1/2)/s');
        C = 50; % Default value
    end

    % Get dry tolerance
    if isfield(cfg.phys, 'dry_tolerance') && ~isempty(cfg.phys.dry_tolerance)
        dry_tol = cfg.phys.dry_tolerance;
    else
        dry_tol = 1e-6; % Default
    end

    % Wet cells (H > dry_tol > 0, so that no division by zero occurs)
    wet = H > dry_tol;
    Hw = H(wet);
    U = HU(wet) ./ Hw;                     % [m/s] Velocity

    % Friction slope of Chézy's law and friction term of the momentum equation
    Sf = abs(U) .* U ./ (C^2 * Hw);        % [-]
    friction_term(wet) = -g * Hw .* Sf;    % [m^2/s^2]

end
