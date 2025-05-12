%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +ic/lake_at_rest.m
%
% Purpose:
%   Generates an initial condition for the 1DWaveTank code representing a
%   lake at rest (still water) with specified level over bathymetry.
%
% Syntax:
%   w0 = lake_at_rest(cfg)
%
% Description:
%   Sets the initial state corresponding to a flat free surface at elevation z = cfg.h0
%   over the bathymetry defined by cfg.bathyHandle. The water is initially at rest (zero velocity).
%
% Inputs:
%   cfg   - [struct] The complete simulation configuration structure. Must contain:
%             cfg.h0           - Still water level (free surface elevation) [m]
%             cfg.mesh.xc      - Cell center coordinates [m]
%             cfg.bathyHandle  - Function handle to calculate bathymetry b(x) [m]
%             cfg.phys.dry_tolerance - Minimum allowed water depth [m]
%
% Outputs:
%   w0    - [2N x 1, double] Flattened initial state vector [H0; HU0], where
%            H0 = max(cfg.h0 - b(xc), cfg.phys.dry_tolerance), HU0 = 0.
%
% Dependencies:
%   Requires the bathymetry function specified by cfg.bathyHandle.
%
% References:
%   - Standard test case for shallow water solvers.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   2025-05-12 (Refactored)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w0 = lake_at_rest(cfg)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameter Extraction and Validation                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Validate required fields in cfg
    if ~isfield(cfg, 'h0') || isempty(cfg.h0)
        error('lake_at_rest:MissingParameter', 'Still water level cfg.h0 must be provided in the configuration.');
    end
    if ~isfield(cfg, 'mesh') || ~isfield(cfg.mesh, 'xc')
        error('lake_at_rest:MissingParameter', 'Mesh coordinates cfg.mesh.xc must be provided.');
    end
     if ~isfield(cfg, 'bathyHandle') || ~isa(cfg.bathyHandle, 'function_handle')
        error('lake_at_rest:MissingParameter', 'Bathymetry function handle cfg.bathyHandle must be provided.');
    end
    if ~isfield(cfg, 'phys') || ~isfield(cfg.phys, 'dry_tolerance') || isempty(cfg.phys.dry_tolerance)
        warning('lake_at_rest:MissingParameter', 'Dry tolerance cfg.phys.dry_tolerance not found. Using default 1e-6.');
        cfg.phys.dry_tolerance = 1e-6;
    end

    xc = cfg.mesh.xc;
    h0 = cfg.h0;
    dry_tol = cfg.phys.dry_tolerance;
    N = numel(xc);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Initial State                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Calculate bathymetry
    b = cfg.bathyHandle(cfg, xc); % Calculate bottom elevation b(x)
    b = b(:); % Ensure column vector
    
    % 2. Calculate initial water depth H0 = h0 - b
    H0 = h0 - b;    % [m] Water depth for flat surface at z=h0
    
    % 3. Apply dry tolerance
    H0 = max(H0, dry_tol); % Ensure depth is at least dry_tolerance
    
    % 4. Initial discharge is zero
    HU0 = zeros(N, 1);             % [m^2/s] Zero initial discharge (still water)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assemble Initial State Vector                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w0 = [H0; HU0];                % [2N x 1]

end