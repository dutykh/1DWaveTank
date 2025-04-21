%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +ic/lake_at_rest.m
%
% Purpose:
%   Generates an initial condition for the 1DWaveTank code representing a
%   lake at rest (still water). The water depth is constant everywhere and
%   the initial velocity/discharge is zero.
%
% Syntax:
%   w0 = lake_at_rest(xc, param)
%
% Inputs:
%   xc    - [vector, double] Cell center coordinates [m].
%   param - [struct] Parameters for the initial condition (optional):
%              param.H0 [double] Background still water depth (default: 0.5)
%
% Outputs:
%   w0    - [2N x 1, double] Flattened initial state vector [H0; HU0], where
%            H0 = param.H0 (water depth), HU0 = 0 (zero velocity everywhere).
%
% Dependencies:
%   None (standalone IC function, but expects xc and param as described).
%
% References:
%   - Standard test case for shallow water solvers.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w0 = lake_at_rest(xc, param)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameter Defaults and Input Handling                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 2 || isempty(param), param = struct(); end
    if ~isfield(param, 'H0'), param.H0 = 0.5; end      % [m] Background depth

    N = numel(xc);
    H0 = param.H0 * ones(N, 1);    % [m] Constant water depth
    HU0 = zeros(N, 1);             % [m^2/s] Zero initial discharge (still water)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assemble Initial State Vector                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial state: H = H0, HU = 0 (lake at rest)
    w0 = [H0; HU0];                % [2N x 1]

end