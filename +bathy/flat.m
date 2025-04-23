%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +bathy/flat.m
%
% Purpose:
%   Implements a flat (constant-depth) bathymetry for the 1DWaveTank code.
%   Returns a vector of constant still-water depths for all spatial points.
%
% Syntax:
%   h = flat(x, cfg)
%
% Inputs:
%   x   - [vector, double] Spatial locations (cell centers or nodes) [m].
%   cfg - [struct] Configuration structure. Should contain cfg.param.H0 (reference depth).
%
% Outputs:
%   h   - [vector, double] Bathymetry (water depth) at each x [m].
%
% Dependencies:
%   Uses cfg.param.H0 (if not present, defaults to 0.5 m).
%
% References:
%   - Standard practice in shallow water modeling.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = flat(x, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Flat bathymetry: all points have the same still-water depth %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is the simplest possible bathymetry, representing a tank
    % or channel with a perfectly horizontal, flat bottom.
    %
    % The depth value H0 is taken from cfg.param.H0 if provided,
    % otherwise defaults to 0.5 m (a typical test value).

    if ~isfield(cfg.param, "H0")
        cfg.param.H0 = 0.5; % [m] Default still water depth
    end

    h = cfg.param.H0 * ones(size(x)); % [m] Constant depth at all x

end % flat