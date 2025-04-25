function limited_slope = umist(delta_minus, delta_plus, cfg)
% UMIST Uniformly Monotonic Interpolation Scheme
%
% UMIST limiter is a smooth limiter with reduced dissipation
% for accurate solution of advection problems.
%
% Syntax:
%   limited_slope = umist(delta_minus, delta_plus)
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

% Use global epsilon from config for division-by-zero protection
epsilon = cfg.numerics.epsilon;

% Check if slopes have different signs
if delta_minus * delta_plus <= 0
    % Opposite signs or one is zero - return zero
    limited_slope = 0;
else
    % Same sign - compute r ratio
    r = delta_plus / (delta_minus + epsilon);
    
    % Compute phi(r) function for UMIST limiter
    phi = max(0, min([2*r, 0.75+0.25*r, 0.25+0.75*r, 2]));
    
    % Apply limiter
    limited_slope = phi * delta_minus;
end
end
