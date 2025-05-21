function h = sloping_beach(cfg, x)

    % SLOPING_BEACH Creates a bathymetry with flat bottom and sloping beach.
%
% Purpose:
%   Defines a bathymetry with a flat portion up to 2/3 of the channel length
%   followed by a constantly sloping beach for runup simulations.
%
% Syntax:
%   h = sloping_beach(cfg, x)
%
% Inputs:
%   cfg    - [struct] Configuration structure containing bathymetry parameters
%   x      - [vector] Spatial coordinates [m]
%
% Required parameters in cfg.bathy.params:
%   L      - [scalar] Channel length [m]
%   h0     - [scalar] Water depth in flat region [m]
%   slope  - [scalar] Beach slope (positive value) [-]
%
% Outputs:
%   h - [vector] Bottom elevation at positions x [m]
%
% Description:
%   The bathymetry consists of a flat portion (h = 0) up to 2/3 of the 
%   channel length, followed by a constantly sloping beach. The slope is 
%   defined as positive (upward).
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
% Date:   21 May 2025

    % Extract parameters from the config structure
    if isfield(cfg, 'bathy') && isfield(cfg.bathy, 'params')
        params = cfg.bathy.params;
    else
        error('sloping_beach:MissingParameters', 'Required bathymetry parameters not found in cfg.bathy.params');
    end
    
    % Get specific parameters
    L = params.L;        % Channel length
    h0 = params.h0;      % Water depth in flat region
    slope = params.slope; % Beach slope (positive value)
    
    % Transition point from flat to sloping bottom (at 2/3 of channel length)
    x_transition = 2*L/3;
    
    % Initialize bottom elevation (negative values for underwater)
    % According to the notation:
    % - Still water level is at y = 0
    % - Bottom elevation is negative for underwater portions
    % - h0 is the water depth in the flat region
    % - Bottom should be at y = -1 in the flat region
    h = zeros(size(x));
    
    % Flat bottom region (underwater, so negative elevation)
    flat_region = x <= x_transition;
    h(flat_region) = -1.0;  % Bottom at y = -1 (underwater)
    
    % Sloping beach region
    sloping_region = x > x_transition;
    % Start at -1.0 and slope upward
    h(sloping_region) = -1.0 + slope * (x(sloping_region) - x_transition);
    
    % No need to enforce non-negativity since we want to allow negative elevations
    % for underwater portions

end