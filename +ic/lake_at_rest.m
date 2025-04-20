function w0 = lake_at_rest(xc, param)
%LAKE_AT_REST Initial condition for a lake at rest.
%   w0 = LAKE_AT_REST(xc, param) returns the initial state vector w0
%   representing a lake at rest (H = H0, HU = 0). H0 is taken from param.H0 (default 0.5).
%   The state vector w0 is flattened [H1;..;HN; HU1;..;HUN].

    if nargin < 2 || isempty(param), param = struct(); end
    if ~isfield(param, 'H0'), param.H0 = 0.5; end

    N = numel(xc);
    H0 = param.H0 * ones(N, 1);
    HU0 = zeros(N, 1);
    w0 = [H0; HU0];

end