function w0 = lake_at_rest(cfg)

    %LAKE_AT_REST Initial condition for a lake at rest.
    %   w0 = LAKE_AT_REST(cfg) returns the initial state vector w0
    %   representing a lake at rest (H = h, HU = 0). The bathymetry 'h' is
    %   evaluated using cfg.bathyHandle at the cell centres cfg.mesh.xc.
    %   The state vector w0 is flattened [H1;..;HN; HU1;..;HUN].
    
    N = cfg.mesh.N;
    xc = cfg.mesh.xc;
    
    % Evaluate bathymetry at cell centers
    h = cfg.bathyHandle(xc); % Get water depth at rest from bathymetry function
    
    % Ensure bathymetry is non-negative
    h = max(h, 0);
    
    % Initial state: H = h, HU = 0
    H0 = h;
    HU0 = zeros(N, 1);
    
    % Flatten state vector: [H; HU]
    w0 = [H0; HU0];

end