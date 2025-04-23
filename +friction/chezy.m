%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +friction/chezy.m
%
% Purpose:
%   Implements the Chézy friction model for shallow water flows.
%   Formula: Sf = U²/(C²H) where C is the Chézy coefficient.
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
%   friction_term - [N x 1, double] Friction source term for momentum equation [m/s^2].
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
    
    % Identify wet cells (only apply friction to wet cells)
    wet_indices = H > dry_tol;
    
    % Calculate velocity in wet cells (avoid division by zero)
    U = zeros(size(H));
    U(wet_indices) = HU(wet_indices) ./ H(wet_indices);
    
    % Apply Chézy formula: Sf = U²/(C²H)
    % Source term needs to be negative (resistance to flow)
    % Final equation: -g*|U|*U/(C²*H)
    friction_term(wet_indices) = -g * abs(U(wet_indices)) .* U(wet_indices) ./ ...
                                (C^2 * H(wet_indices) + 1e-10); % Add small constant to prevent division by zero
                                
    % Return momentum source term (applied to HU)

end