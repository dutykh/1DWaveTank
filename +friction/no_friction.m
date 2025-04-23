%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +friction/no_friction.m
%
% Purpose:
%   Default friction model that applies no friction (returns zeros).
%   Used when no friction effects are desired in the simulation.
%
% Syntax:
%   friction_term = no_friction(H, HU, g, cfg)
%
% Inputs:
%   H            - [N x 1, double] Water depth at each cell [m].
%   HU           - [N x 1, double] Discharge at each cell [m^2/s].
%   g            - [scalar, double] Gravitational acceleration [m/s^2].
%   cfg          - [struct] Configuration structure (not used in this model).
%
% Outputs:
%   friction_term - [N x 1, double] Friction source term for momentum equation [m/s^2].
%                   Returns zeros (no friction).
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: 23 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function friction_term = no_friction(H, HU, g, cfg)

    % Return zero friction for all cells
    friction_term = zeros(size(HU));

end