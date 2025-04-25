%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +friction/darcy_weisbach.m
%
% Purpose:
%   Implements the Darcy-Weisbach friction model for shallow water flows.
%   Formula: Sf = f|U|U/(8gH) where f is the Darcy friction factor.
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
%   friction_term - [N x 1, double] Friction source term for momentum equation [m/s^2].
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
    
    % Identify wet cells (only apply friction to wet cells)
    wet_indices = H > dry_tol;
    
    % Calculate velocity in wet cells
    U = zeros(size(H));
    U(wet_indices) = HU(wet_indices) ./ H(wet_indices);
    
    % Determine how to calculate friction factor f
    if isfield(cfg.phys, 'f_calculation') && strcmpi(cfg.phys.f_calculation, 'colebrook_white')
        % Use Colebrook-White formula to calculate f for each cell
        f = zeros(size(H));
        f(wet_indices) = friction.colebrook_white(H(wet_indices), U(wet_indices), cfg);
    else
        % Use constant f value (default)
        if isfield(cfg.phys, 'darcy_f') && ~isempty(cfg.phys.darcy_f)
            f = cfg.phys.darcy_f * ones(size(H));
        else
            warning('Darcy friction factor not specified. Using default f = 0.02');
            f = 0.02 * ones(size(H));
        end
    end
    
    % Apply Darcy-Weisbach formula: Sf = f|U|U/(8gH)
    % Source term needs to be negative (resistance to flow)
    friction_term(wet_indices) = -g * f(wet_indices) .* abs(U(wet_indices)) .* U(wet_indices) ./ ...
                               (8 * H(wet_indices) + cfg.numerics.epsilon); % Use config epsilon to prevent division by zero
end
