%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +ic/dry_dam_break.m
%
% Purpose:
%   Generates an initial condition for the 1DWaveTank code representing a
%   dry dam break problem. The water is initially at rest with a discontinuity
%   in depth: wet on the left side and dry on the right side.
%
% Syntax:
%   w0 = dry_dam_break(xc, cfg)
%
% Inputs:
%   xc    - [vector, double] Cell center coordinates [m].
%   cfg   - [struct] Configuration structure containing:
%           cfg.dry_dam_break.h_L [double] Water depth left of dam (default: 0.5 m)
%           cfg.dry_dam_break.x_dam [double] Dam location (default: domain center)
%           cfg.phys.dry_tolerance [double] Threshold for dry states (default: 1e-6 m)
%
% Outputs:
%   w0    - [2N x 1, double] Flattened initial state vector [H0; HU0].
%
% Dependencies:
%   None (standalone IC function).
%
% References:
%   - Standard test case for shallow water solvers with wetting and drying.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
% Date:   April 26, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w0 = dry_dam_break(xc, cfg)

    % Extract parameters or use defaults
    if isfield(cfg, 'dry_dam_break') && isfield(cfg.dry_dam_break, 'h_L')
        h_L = cfg.dry_dam_break.h_L;
    else
        h_L = 0.5; % Default water depth on left side [m]
    end
    
    if isfield(cfg, 'dry_dam_break') && isfield(cfg.dry_dam_break, 'x_dam')
        x_dam = cfg.dry_dam_break.x_dam;
    else
        % Default dam location at domain center
        x_dam = 0.5 * (min(xc) + max(xc));
    end
    
    % Use a small epsilon for dry region
    if isfield(cfg, 'phys') && isfield(cfg.phys, 'dry_tolerance')
        epsilon = cfg.phys.dry_tolerance;
    else
        epsilon = 1e-6; % Default dry tolerance
    end
    
    % Number of spatial cells
    N = length(xc);
    
    % Initialize water depth vector
    H0 = zeros(N, 1);
    
    % Set water depth based on dam position
    for i = 1:N
        if xc(i) < x_dam
            H0(i) = h_L; % Left of dam (wet)
        else
            H0(i) = epsilon; % Right of dam (dry)
        end
    end
    
    % Initialize discharge to zero (water at rest initially)
    HU0 = zeros(N, 1);
    
    % Combine into flattened state vector
    w0 = [H0; HU0];

end