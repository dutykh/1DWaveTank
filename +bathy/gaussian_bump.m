function b = gaussian_bump(cfg, x)

%GAUSSIAN_BUMP Defines a Gaussian bump bathymetry.
%   b = gaussian_bump(cfg, x) returns the bathymetry profile b(x)
%   representing a flat bottom with a Gaussian bump centered in the domain.
%   The bathymetry 'b' represents the bottom elevation relative to z=0.
%   Still water level is at z = cfg.h0. Water depth is eta - b.
%
%   Inputs:
%       cfg - Configuration structure containing parameters. Expected fields:
%             cfg.h0                - Still water depth (required).
%             cfg.L                 - Domain length (required).
%             cfg.bathy_bump_center - Center position of the bump (optional, default: L/2).
%             cfg.bathy_bump_height - Height of the bump (optional, default: 0.2 * h0).
%             cfg.bathy_bump_width  - Characteristic width (std dev) of the bump (optional, default: L/10).
%             cfg.min_depth         - Minimum allowed water depth (optional, used for warning/capping).
%       x   - Vector of spatial coordinates.
%
%   Output:
%       b   - Vector of bathymetry elevation values at coordinates x.
%
%   Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi, UAE)
%   Date:   2025-05-12

% --- Input Validation & Default Parameters ---
% Robustly extract h0 and L (accepts config.h0 or config.param.H0, config.L or config.domain.xmax-xmin)
if isfield(cfg, 'h0') && ~isempty(cfg.h0)
    h0 = cfg.h0;
elseif isfield(cfg, 'param') && isfield(cfg.param, 'H0') && ~isempty(cfg.param.H0)
    h0 = cfg.param.H0;
else
    error('gaussian_bump:MissingParameter', 'Still water depth cfg.h0 (or cfg.param.H0) must be provided.');
end
if isfield(cfg, 'L') && ~isempty(cfg.L)
    L = cfg.L;
elseif isfield(cfg, 'domain') && isfield(cfg.domain, 'xmax') && isfield(cfg.domain, 'xmin')
    L = cfg.domain.xmax - cfg.domain.xmin;
else
    error('gaussian_bump:MissingParameter', 'Domain length cfg.L (or cfg.domain.xmax-xmin) must be provided.');
end

% Use these robust values throughout

% Default parameters
default_center = L / 2;
default_height = 0.2 * h0;
default_width  = L / 10;
default_min_depth = 1e-6; % Default minimum depth to avoid division by zero etc.

bump_center = default_center;
if isfield(cfg, 'bathy_bump_center') && ~isempty(cfg.bathy_bump_center)
    bump_center = cfg.bathy_bump_center;
end

bump_height = default_height;
if isfield(cfg, 'bathy_bump_height') && ~isempty(cfg.bathy_bump_height)
    bump_height = cfg.bathy_bump_height;
    if bump_height < 0
        warning('gaussian_bump:InvalidHeight', 'Specified bump height is negative. Using absolute value.');
        bump_height = abs(bump_height);
    end
end

bump_width = default_width;
if isfield(cfg, 'bathy_bump_width') && ~isempty(cfg.bathy_bump_width)
    bump_width = cfg.bathy_bump_width;
     if bump_width <= 0
        warning('gaussian_bump:InvalidWidth', 'Specified bump width must be positive. Using default.');
        bump_width = default_width;
    end
end

min_depth = default_min_depth;
 if isfield(cfg, 'min_depth') && ~isempty(cfg.min_depth)
    min_depth = cfg.min_depth;
end
% --- End Validation ---

% --- Calculate Bathymetry ---
% Gaussian bump centered at 'bump_center' with height 'bump_height'
% and standard deviation 'bump_width'.
% b(x) = bump_height * exp(-(x - bump_center)^2 / (2 * bump_width^2))
b = -bump_height * exp(-((x - bump_center).^2) / (2 * bump_width^2));
% --- End Calculation ---

% --- Check for excessive bump height ---
% Water depth at rest over the bump is h0 - b. Ensure it's > min_depth.
max_b = max(b);
if max_b >= cfg.h0 - min_depth
     warning('gaussian_bump:PotentialDryArea', ...
            ['Maximum bump height (%.2f) results in water depth <= min_depth (%.2e) ',...
             'at some points for still water level h0=%.2f. Simulation might become unstable.'], ...
             max_b, min_depth, cfg.h0);
    % Optional: Cap the bump height to ensure minimum depth
    % scaling_factor = (cfg.h0 - min_depth) / max_b;
    % b = b * scaling_factor;
    % fprintf('Bump height capped to %.2f to maintain minimum depth.\n', max(b));
end
% --- End Check ---

end