%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +core/+utils/uniform.m
%
% Purpose:
%   Generates a uniform 1D finite volume mesh (cell centers, cell edges,
%   and spatial step size) based on the specified domain limits and number
%   of cells. This is a common utility for setting up the spatial discretization
%   for finite volume methods.
%
% Syntax:
%   [xc, dx, x_edge] = uniform(domain, N)
%
% Inputs:
%   domain - [struct] Structure defining the domain boundaries.
%            Required fields:
%              domain.xmin: [double] Left boundary coordinate [m].
%              domain.xmax: [double] Right boundary coordinate [m].
%   N      - [integer] Number of finite volume cells (control volumes) to create.
%
% Outputs:
%   xc     - [N x 1, double] Vector of cell center coordinates [m].
%   dx     - [scalar, double] Constant spatial step size (cell width) [m].
%   x_edge - [N+1 x 1, double] Vector of cell edge coordinates [m].
%
% Dependencies:
%   None (uses built-in MATLAB functions).
%
% References:
%   - LeVeque, R.J. (2002). Finite Volume Methods for Hyperbolic Problems.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xc, dx, x_edge] = uniform(domain, N)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Validation                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isstruct(domain) || ~isfield(domain, 'xmin') || ~isfield(domain, 'xmax')
        error('uniform:InputError', 'Input ''domain'' must be a structure with ''xmin'' and ''xmax'' fields.');
    end
    if ~isscalar(N) || ~isnumeric(N) || N <= 0 || mod(N,1) ~= 0
        error('uniform:InputError', 'Input ''N'' must be a positive integer scalar.');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Mesh Properties                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the N+1 cell edges uniformly distributed across the domain.
    x_edge = linspace(domain.xmin, domain.xmax, N+1)'; % [N+1 x 1, m] Cell edges

    % Calculate the constant spatial step size (cell width).
    dx = (domain.xmax - domain.xmin) / N; % [scalar, m] Spatial step

    % Calculate the N cell centers, located halfway between adjacent edges.
    % xc(i) = 0.5 * (x_edge(i) + x_edge(i+1))
    xc = x_edge(1:N) + 0.5 * dx; % [N x 1, m] Cell centres

end