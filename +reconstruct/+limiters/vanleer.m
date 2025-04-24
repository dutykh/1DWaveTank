function limited_slope = vanleer(delta_minus, delta_plus)
% VANLEER van Leer slope limiter for TVD schemes
%
% The van Leer limiter uses a harmonic mean of the slope estimates,
% which produces smoother solutions than minmod but is less aggressive
% than superbee.
%
% Syntax:
%   limited_slope = vanleer(delta_minus, delta_plus)
%
% Inputs:
%   delta_minus - Backward difference (q_i - q_{i-1})
%   delta_plus  - Forward difference (q_{i+1} - q_i)
%
% Output:
%   limited_slope - Limited slope

% Check if slopes have different signs
if delta_minus * delta_plus <= 0
    % Opposite signs or one is zero - return zero
    limited_slope = 0;
else
    % Same sign - use harmonic mean
    limited_slope = 2 * delta_minus * delta_plus / (delta_minus + delta_plus);
end
end
