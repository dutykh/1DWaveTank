%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +friction/manning.m
%
% Purpose:
%   Implements the Manning friction model for shallow water flows.
%   Formula: Sf = n²|U|U/H^(4/3) where n is Manning's roughness coefficient.
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
%   friction_term - [N x 1, double] Friction source term for momentum equation [m/s^2].
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
    
    % Identify wet cells (only apply friction to wet cells)
    wet_indices = H > dry_tol;

    % Use global epsilon from config for division-by-zero protection
    epsilon = cfg.numerics.epsilon;
    
    % Calculate velocity in wet cells (avoid division by zero)
    U = zeros(size(H));
    U(wet_indices) = HU(wet_indices) ./ H(wet_indices);
    
    % Apply Manning formula: Sf = n²|U|U/H^(4/3)
    % Source term needs to be negative (resistance to flow)
    % Final equation: -g*n²*|U|*U/H^(4/3)
    friction_term(wet_indices) = -g * n^2 * abs(U(wet_indices)) .* U(wet_indices) ./ ...
                                (H(wet_indices).^(4/3) + epsilon); % Use epsilon from config to prevent division by zero

end