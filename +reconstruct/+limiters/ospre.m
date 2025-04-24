function limited_slope = ospre(delta_minus, delta_plus)
% OSPRE Smooth blended limiter for TVD schemes
%
% OSPRE limiter provides a smooth blend between minmod and superbee,
% giving good results for many practical problems.
%
% Syntax:
%   limited_slope = ospre(delta_minus, delta_plus)
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

% Small value to avoid division by zero
epsilon = 1e-12;

% Check if slopes have different signs
if delta_minus * delta_plus <= 0
    % Opposite signs or one is zero - return zero
    limited_slope = 0;
else
    % Same sign - compute r ratio
    r = delta_plus / (delta_minus + epsilon);
    
    % Compute phi(r) function for OSPRE limiter
    phi = 1.5 * (r^2 + r) / (r^2 + r + 1);
    
    % Apply limiter
    limited_slope = phi * delta_minus;
end
end
