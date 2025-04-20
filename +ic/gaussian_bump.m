function w0 = gaussian_bump(cfg)

    %GAUSSIAN_BUMP Initial condition with a Gaussian bump on the water surface.
    %   w0 = GAUSSIAN_BUMP(cfg) returns the initial state vector w0
    %   representing H = h + a*exp(-lambda*(x - x0)^2), HU = 0.
    %   Bathymetry 'h' is from cfg.bathyHandle.
    %   Parameters for the bump are taken from cfg.ic.param using
    %   core.utils.get_param:
    %       cfg.ic.param.a      (amplitude, default 0.25)
    %       cfg.ic.param.lambda (decay rate, default 0.1)
    %       cfg.ic.param.x0     (center position, default domain center)
    %
    %   Outputs:
    %       w0 - Flattened initial state vector [H0; HU0].

    N = cfg.mesh.N;
    xc = cfg.mesh.xc;
    xmin = cfg.domain.xmin;
    xmax = cfg.domain.xmax;

    % --- Get parameters with defaults using the utility function ---
    if isfield(cfg, 'ic') && isfield(cfg.ic, 'param')
        params = cfg.ic.param;
    else
        params = struct(); % Use defaults if params substructure not provided
    end

    % Use core.utils.get_param to retrieve parameters
    a = core.utils.get_param(params, 'a', 0.25);           % Amplitude
    lambda = core.utils.get_param(params, 'lambda', 0.1);  % Decay rate
    x0 = core.utils.get_param(params, 'x0', (xmin + xmax) / 2); % Center

    % Evaluate bathymetry at cell centers
    h = cfg.bathyHandle(xc, cfg);
    h = max(h, 0); % Ensure non-negative bathymetry

    % Calculate Gaussian perturbation
    eta = a * exp(-lambda * (xc - x0).^2);

    % Initial state: H = h + eta, HU = 0
    H0 = h + eta;
    H0 = max(H0, 0); % Ensure non-negative water depth
    HU0 = zeros(N, 1);

    % Flatten state vector: [H; HU]
    w0 = [H0; HU0];

end