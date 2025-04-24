function limited_slope = mc(delta_minus, delta_plus)
% MC Monotonized Central slope limiter for TVD schemes
%
% The MC limiter provides a good balance between accuracy and stability,
% less dissipative than minmod but more stable than superbee.
%
% Syntax:
%   limited_slope = mc(delta_minus, delta_plus)
%
% Inputs:
%   delta_minus - Backward difference (q_i - q_{i-1})
%   delta_plus  - Forward difference (q_{i+1} - q_i)
%
% Output:
%   limited_slope - Limited slope
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: April 24, 2025

% Check if slopes have different signs
if delta_minus * delta_plus <= 0
    % Opposite signs or one is zero - return zero
    limited_slope = 0;
else
    % Same sign - compute MC limited slope
    limited_slope = sign(delta_plus) * min([
        2*abs(delta_minus), 
        2*abs(delta_plus), 
        0.5*abs(delta_minus+delta_plus)
    ]);
end
end
