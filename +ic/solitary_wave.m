%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +ic/solitary_wave.m
%
% Purpose:
%   Generates an initial condition for the 1DWaveTank code representing a
%   solitary wave on a flat bottom.
%
% Syntax:
%   w0 = solitary_wave(xc, cfg)
%
% Inputs:
%   xc    - [vector, double] Cell center coordinates [m].
%   cfg   - [struct] Configuration structure. Must contain:
%           cfg.domain.xmin - [double] Left domain boundary
%           cfg.domain.xmax - [double] Right domain boundary
%           cfg.phys.g      - [double] Gravity acceleration
%         Optional parameters (in cfg.ic_param or cfg.param):
%           cfg.ic_param.h0 or cfg.param.h0 - [double] Still water depth (default: 0.5)
%           cfg.ic_param.a  or cfg.param.a  - [double] Wave amplitude (default: 0.2)
%
% Outputs:
%   w0    - [2N x 1, double] Flattened initial state vector [H0; HU0].
%
% Formulas:
%   H(x) = h0 + eta(x)
%   eta(x) = a*sech(0.5*k*(x - x0))^2
%   u(x) = c*eta(x)/H(x)
%   c = sqrt(g*(h0 + a))
%   k = sqrt((3*a)/(h0 + a))/h0
%
% Dependencies:
%   None (standalone IC function).
%
% References:
%   - Boussinesq theory for solitary waves.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   22 April 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w0 = solitary_wave(xc, cfg)
    % disp(['Size of xc in solitary_wave: ', mat2str(size(xc))]); % DEBUG: Print size of input xc

    % --- Parameters ---
    g = cfg.phys.g; % Gravity

    % Get h0 (still water depth) - check ic_param first, then param, then default
    if isfield(cfg, 'ic_param') && isfield(cfg.ic_param, 'h0')
        h0 = cfg.ic_param.h0;
    elseif isfield(cfg, 'param') && isfield(cfg.param, 'h0')
        h0 = cfg.param.h0;
    else
        h0 = 0.5; % Default still water depth
        fprintf('Using default h0 = %.2f for solitary wave IC.\n', h0);
    end

    % Get a (wave amplitude) - check ic_param first, then param, then default
    if isfield(cfg, 'ic_param') && isfield(cfg.ic_param, 'a')
        a = cfg.ic_param.a;
    elseif isfield(cfg, 'param') && isfield(cfg.param, 'a')
        a = cfg.param.a;
    else
        a = 0.2; % Default amplitude
        fprintf('Using default a = %.2f for solitary wave IC.\n', a);
    end

    % Center position
    x0 = 0.5 * (cfg.domain.xmin + cfg.domain.xmax);

    % --- Derived parameters ---
    c = sqrt(g*(h0 + a));     % Wave celerity
    k = sqrt(3*a / (h0^3 * (1 + a/h0))); % Corrected k based on common form sqrt(3a/h0^3)
    % Simplified k = sqrt((3*a)/(h0 + a))/h0; % User provided formula

    % --- Calculate Initial Profile ---
    N = numel(xc);
    eta = a * sech(0.5 * k * (xc - x0)).^2; % Free surface elevation profile
    eta = eta(:);                            % Ensure eta is a column vector
    H0 = h0 + eta;                         % Total water depth H = h0 + eta
    H0 = H0(:);                            % Ensure H0 is a column vector

    % Calculate velocity U
    U0 = zeros(N, 1);
    wet_indices = H0 > 1e-6; % Avoid division by zero in dry cells
    U0(wet_indices) = c * eta(wet_indices) ./ H0(wet_indices);

    HU0 = H0 .* U0; % Initial discharge

    % --- Format Output ---
    w0 = [H0(:); HU0(:)]; % Flattened column vector [H; HU]

end
