%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +bathy/flat.m
%
% Purpose: Implements a flat bottom elevation, meaning the bottom elevation
%          is constant. By default, this is 0.0 m, representing the z=0 datum.
%          This function now returns a bottom elevation (e.g., 0 by default,
%          or configurable via cfg.bathy_params.flat_elevation or
%          cfg.param.flat_bottom_elevation). The parameter cfg.param.H0
%          is NO LONGER used by this function to define a depth.
%
% Syntax:
%   h = flat(x, cfg)
%
% Inputs:
%   x   - [vector, double] Spatial locations (cell centers or nodes) [m].
%   cfg - [struct] Configuration structure.
%
% Outputs:
%   h   - [vector, double] Bottom elevation $z_b(x)$ at each x [m].
%
% Dependencies:
%   None.
%
% References:
%   - Standard practice in shallow water modeling.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = flat(x, cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Flat bathymetry: all points have the same bottom elevation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is the simplest possible bathymetry, representing a tank
    % or channel with a perfectly horizontal, flat bottom.
    %
    % This function returns a constant bottom elevation. By default, this is 0.0,
    % representing a flat bottom at the z=0 datum. This can be configured
    % via cfg.bathy_params.flat_elevation if a different constant elevation is needed.
    % The parameter cfg.param.H0 is NO LONGER used by this function to define a depth.

    % Default flat bottom elevation (e.g., at z=0 datum)
    bottom_elevation = 0.0; 

    % Allow overriding default via a new configuration parameter if needed
    if isfield(cfg, 'bathy_params') && isfield(cfg.bathy_params, 'flat_elevation')
        bottom_elevation = cfg.bathy_params.flat_elevation;
    elseif isfield(cfg, 'param') && isfield(cfg.param, 'flat_bottom_elevation') % Alternative check
         bottom_elevation = cfg.param.flat_bottom_elevation;
    end
    bottom_elevation = cfg.bathy_params.flat_elevation;
elseif isfield(cfg, 'param') && isfield(cfg.param, 'flat_bottom_elevation') % Alternative check
     bottom_elevation = cfg.param.flat_bottom_elevation;
end

h = bottom_elevation * ones(size(x)); % [m] Constant bottom elevation

end % flat