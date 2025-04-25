function limited_slope = superbee(delta_minus, delta_plus)
% SUPERBEE Superbee slope limiter for TVD schemes
%
% The superbee limiter is the most aggressive TVD limiter, 
% choosing the maximum slope that won't violate monotonicity.
% Author: Dr. Denys Dutykh
%
% Syntax:
%   limited_slope = superbee(delta_minus, delta_plus)
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
    % Same sign - compute the two candidate slopes
    slope1 = sign(delta_minus) * min(abs(delta_minus), 2*abs(delta_plus));
    slope2 = sign(delta_plus) * min(2*abs(delta_minus), abs(delta_plus));
    
    % Choose the maximum magnitude
    if abs(slope1) > abs(slope2)
        limited_slope = slope1;
    else
        limited_slope = slope2;
    end
end
end
