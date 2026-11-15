function limited_slope = sweby(delta_minus, delta_plus, cfg_or_beta)
% SWEBY Parameterized TVD limiter
%
% Sweby's limiter incorporates a parameter beta which allows it
% to be adjusted between minmod (beta=1) and superbee (beta=2):
%     phi(r) = max(0, min(beta r, 1), min(r, beta)),   1 ⩽ beta ⩽ 2,
% equation (3.17) of Sweby (1984), with phi(r)/r = phi(1/r). With
% r = delta_plus/delta_minus, the limited slope phi(r) delta_minus is
% written below without the quotient r,
%     slope = s max(min(beta |delta_plus|, |delta_minus|), min(|delta_plus|, beta |delta_minus|))
% when delta_minus delta_plus > 0 (s the common sign), and 0 otherwise.
%
% Syntax:
%   limited_slope = sweby(delta_minus, delta_plus)
%   limited_slope = sweby(delta_minus, delta_plus, cfg)
%   limited_slope = sweby(delta_minus, delta_plus, beta)
%
% Inputs:
%   delta_minus - Backward difference (q_i - q_{i-1})
%   delta_plus  - Forward difference (q_{i+1} - q_i)
%   cfg_or_beta - (optional) the configuration, as the MUSCL reconstruction
%                 passes it to every limiter, whose field
%                 cfg.reconstruct.sweby_beta gives beta; or beta itself.
%                 beta = 1.5 when neither gives it.
%                 beta=1: equivalent to minmod
%                 beta=2: equivalent to superbee
%
% Output:
%   limited_slope - Limited slope
%
% References:
%   - P. K. Sweby, High resolution schemes using flux limiters for
%     hyperbolic conservation laws, SIAM J. Numer. Anal. 21(5), 995-1011
%     (1984), equation (3.17).
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: April 24, 2025

% Parameter beta: from the configuration, or given directly, or 1.5
beta = 1.5; % Good compromise between minmod and superbee
if nargin >= 3
    if isnumeric(cfg_or_beta) && isscalar(cfg_or_beta)
        beta = cfg_or_beta;
    elseif isstruct(cfg_or_beta) && isfield(cfg_or_beta, 'reconstruct') ...
            && isfield(cfg_or_beta.reconstruct, 'sweby_beta') && ~isempty(cfg_or_beta.reconstruct.sweby_beta)
        beta = cfg_or_beta.reconstruct.sweby_beta;
    end
end

% Check if slopes have different signs
if delta_minus * delta_plus <= 0
    % Opposite signs or one is zero - return zero
    limited_slope = 0;
else
    % Same sign: phi(r) delta_minus with r = delta_plus/delta_minus, in a form
    % that needs no division
    a = abs(delta_minus);
    b = abs(delta_plus);
    limited_slope = sign(delta_minus) * max(min(beta*b, a), min(b, beta*a));
end
end
