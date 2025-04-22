%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +bathy/flat_bathymetry.m
%
% Purpose:
%   Defines a flat bathymetry profile (z=0) for the 1DWaveTank simulation.
%
% Syntax:
%   b = flat_bathymetry(xc, cfg)
%
% Inputs:
%   xc    - [vector, double] Cell center coordinates [m].
%   cfg   - [struct] Configuration structure (unused for flat bathymetry,
%             but included for consistency with other bathy functions).
%
% Outputs:
%   b     - [vector, double] Bathymetry elevation vector [m] (all zeros),
%           same size as xc.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
%         Implemented by Cascade (Codeium)
% Date:   22 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = flat_bathymetry(xc, cfg) %#ok<INUSD> cfg is unused
    % Return a vector of zeros with the same size as the input xc vector.
    % This represents a flat bottom at elevation z = 0.
    b = zeros(size(xc));
end
