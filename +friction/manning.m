%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +friction/manning.m
%
% Purpose:
%   Implements the Manning (Gauckler-Manning-Strickler) friction law for
%   shallow water flows. The momentum equation of the depth-averaged flow
%   over a wide channel reads
%       (HU)_t + (HU^2 + g H^2/2)_x = - g H z_x - g H S_f,
%   where the friction slope of Manning's law is
%       S_f = n^2 U|U| / H^(4/3)                                [-],
%   n being Manning's roughness coefficient [s/m^(1/3)] and H the hydraulic
%   radius of a wide channel. This function returns the friction term of the
%   momentum equation,
%       - g H S_f = - g n^2 U|U| / H^(1/3)                      [m^2/s^2].
%   In uniform flow down a slope S_0, S_f = S_0 gives the normal velocity
%   U = H^(2/3) S_0^(1/2) / n.
%
% Syntax:
%   friction_term = manning(H, HU, g, cfg)
%
% Inputs:
%   H            - [N x 1, double] Water depth at each cell [m].
%   HU           - [N x 1, double] Discharge at each cell [m^2/s].
%   g            - [scalar, double] Gravitational acceleration [m/s^2].
%   cfg          - [struct] Configuration structure containing:
%                  cfg.phys.manning_n: [double] Manning's roughness coefficient [s/m^(1/3)].
%                  cfg.phys.dry_tolerance: [double] Threshold for dry cells [m].
%
% Outputs:
%   friction_term - [N x 1, double] Friction term - g H S_f of the momentum
%                   equation, added to d(HU)/dt [m^2/s^2] (zero in dry cells).
%
% References:
%   - Castro-Orgaz, O., Hager, W.H. (2019). Shallow Water Hydraulics.
%     Springer, Section 3.3.2, equation (3.16): S_f = n^2 U^2 / R_h^(4/3)
%     in SI units, with R_h = H in a wide channel.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: 23 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function friction_term = manning(H, HU, g, cfg)

    % Initialize output vector
    friction_term = zeros(size(HU));

    % Get Manning's coefficient from config
    if isfield(cfg.phys, 'manning_n') && ~isempty(cfg.phys.manning_n)
        n = cfg.phys.manning_n;
    else
        warning('Manning coefficient not specified. Using default n = 0.03 s/m^(1/3)');
        n = 0.03; % Default value for a clean, straight channel
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

    % Friction slope of Manning's law and friction term of the momentum equation
    Sf = n^2 * abs(U) .* U ./ Hw.^(4/3);   % [-]
    friction_term(wet) = -g * Hw .* Sf;    % [m^2/s^2]

end
