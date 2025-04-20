function w0 = gaussian_bump(xc, param)

    %GAUSSIAN_BUMP Initial condition with a Gaussian bump on the water surface.
    %   w0 = GAUSSIAN_BUMP(xc, param) returns the initial state vector w0
    %   representing H = h + a*exp(-lambda*(x - x0)^2), HU = 0.
    %   Bathymetry h is assumed flat with depth param.H0 (default 0.5).
    %   Parameters for the bump:
    %       param.a      (amplitude, default 0.25)
    %       param.lambda (decay rate, default 0.1)
    %       param.x0     (center position, default domain center)
    %       param.H0     (background depth, default 0.5)
    %   Outputs:
    %       w0 - Flattened initial state vector [H0; HU0].

    if nargin < 2 || isempty(param), param = struct(); end
    if ~isfield(param, 'a'),      param.a = 0.25; end
    if ~isfield(param, 'lambda'), param.lambda = 0.1; end
    if ~isfield(param, 'H0'),    param.H0 = 0.5; end
    if ~isfield(param, 'x0')
        if ~isempty(xc)
            param.x0 = mean([min(xc), max(xc)]);
        else
            param.x0 = 0;
        end
    end

    N = numel(xc);
    h = param.H0 * ones(size(xc));
    a = param.a;
    lambda = param.lambda;
    x0 = param.x0;

    % Calculate Gaussian perturbation
    eta = a * exp(-lambda * (xc - x0).^2);

    % Initial state: H = h + eta, HU = 0
    H0 = h + eta;
    H0 = max(H0, 0); % Ensure non-negative water depth
    HU0 = zeros(N, 1);

    % Flatten state vector: [H; HU]
    w0 = [H0; HU0];

end