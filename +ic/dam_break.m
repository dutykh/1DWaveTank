%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +ic/dam_break.m
%
% Purpose:
%   Generates an initial condition for the 1DWaveTank code representing a
%   classical dam break problem. The water is initially at rest with a
%   discontinuity in depth (higher water level on left side).
%
% Syntax:
%   w0 = dam_break(xc, cfg)
%
% Inputs:
%   xc    - [vector, double] Cell center coordinates [m].
%   cfg   - [struct] Configuration structure containing:
%           cfg.dam_break.h_L [double] Water depth left of dam (default: 0.8 m)
%           cfg.dam_break.h_R [double] Water depth right of dam (default: 0.5 m)
%           cfg.dam_break.x_dam [double] Dam location (default: domain center)
%
% Outputs:
%   w0    - [2N x 1, double] Flattened initial state vector [H0; HU0].
%
% Dependencies:
%   None (standalone IC function).
%
% References:
%   - Standard test case for shallow water solvers.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   April 24, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w0 = dam_break(xc, cfg)
    % Extract parameters or use defaults
    if isfield(cfg, 'dam_break') && isfield(cfg.dam_break, 'h_L')
        h_L = cfg.dam_break.h_L;
    else
        h_L = 0.8; % Default water depth on left side [m]
    end
    
    if isfield(cfg, 'dam_break') && isfield(cfg.dam_break, 'h_R')
        h_R = cfg.dam_break.h_R;
    else
        h_R = 0.5; % Default water depth on right side [m]
    end
    
    if isfield(cfg, 'dam_break') && isfield(cfg.dam_break, 'x_dam')
        x_dam = cfg.dam_break.x_dam;
    else
        % Default dam location at domain center
        x_dam = 0.5 * (min(xc) + max(xc));
    end
    
    % Number of spatial cells
    N = length(xc);
    
    % Initialize water depth vector
    H0 = zeros(N, 1);
    
    % Set water depth based on dam position
    for i = 1:N
        if xc(i) < x_dam
            H0(i) = h_L; % Left of dam
        else
            H0(i) = h_R; % Right of dam
        end
    end
    
    % Initialize discharge to zero (water at rest initially)
    HU0 = zeros(N, 1);
    
    % Combine into flattened state vector
    w0 = [H0; HU0];
end
