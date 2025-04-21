%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +ic/gaussian_bump.m
%
% Purpose:
%   Generates an initial condition for the 1DWaveTank code with a Gaussian
%   bump in the water surface elevation (eta) and zero initial velocity.
%   The background bathymetry is assumed flat.
%
% Syntax:
%   w0 = gaussian_bump(xc, param)
%
% Inputs:
%   xc    - [vector, double] Cell center coordinates [m].
%   param - [struct] Parameters for the Gaussian bump (fields are optional):
%              param.a      [double] Amplitude of bump (default: 0.25)
%              param.lambda [double] Decay rate of bump (default: 0.1)
%              param.x0     [double] Center position of bump (default: domain center)
%              param.H0     [double] Background still water depth (default: 0.5)
%
% Outputs:
%   w0    - [2N x 1, double] Flattened initial state vector [H0; HU0], where
%            H0 = h + eta (water depth), HU0 = 0 (zero velocity everywhere).
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

function w0 = gaussian_bump(xc, param)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameter Defaults and Input Handling                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 2 || isempty(param), param = struct(); end
    if ~isfield(param, 'a'),      param.a = 0.25; end           % [m] Amplitude
    if ~isfield(param, 'lambda'), param.lambda = 0.1; end       % [1/m^2] Decay rate
    if ~isfield(param, 'H0'),    param.H0 = 0.5; end            % [m] Background depth
    if ~isfield(param, 'x0')
        if ~isempty(xc)
            param.x0 = mean([min(xc), max(xc)]);                % [m] Center at domain middle
        else
            param.x0 = 0;
        end
    end

    N = numel(xc);
    h = param.H0 * ones(size(xc)); % [m] Flat bathymetry
    a = param.a;                   % [m] Amplitude
    lambda = param.lambda;         % [1/m^2] Decay rate
    x0 = param.x0;                 % [m] Center position

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Gaussian Bump in Surface Elevation                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % eta(x) = a * exp(-lambda * (x - x0)^2)
    eta = a * exp(-lambda * (xc - x0).^2); % [m]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assemble Initial State Vector                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial state: H = h + eta, HU = 0 (still water velocity)
    H0 = h + eta;              % [m] Initial water depth
    H0 = max(H0, 0);           % Ensure non-negative water depth (no dry cells)
    HU0 = zeros(N, 1);         % [m^2/s] Zero initial discharge

    % Flatten state vector: [H; HU]
    w0 = [H0(:); HU0];            % [2N x 1]

end