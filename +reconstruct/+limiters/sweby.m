function limited_slope = sweby(delta_minus, delta_plus, beta)
% SWEBY Parameterized TVD limiter
%
% Sweby's limiter incorporates a parameter beta which allows it
% to be adjusted between minmod (beta=1) and superbee (beta=2).
%
% Syntax:
%   limited_slope = sweby(delta_minus, delta_plus, beta)
%
% Inputs:
%   delta_minus - Backward difference (q_i - q_{i-1})
%   delta_plus  - Forward difference (q_{i+1} - q_i)
%   beta        - Parameter controlling limiter behavior (default 1.5)
%                 beta=1: equivalent to minmod
%                 beta=2: equivalent to superbee
%
% Output:
%   limited_slope - Limited slope
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date: April 24, 2025

% Default beta value if not provided
if nargin < 3
    beta = 1.5; % Good compromise between minmod and superbee
end

% Small value to avoid division by zero
epsilon = 1e-12;

% Check if slopes have different signs
if delta_minus * delta_plus <= 0
    % Opposite signs or one is zero - return zero
    limited_slope = 0;
else
    % Same sign - compute r ratio
    r = delta_plus / (delta_minus + epsilon);
    
    % Compute phi(r) function for Sweby's limiter
    phi = max(0, min([beta*r, beta, r]));
    
    % Apply limiter
    limited_slope = phi * delta_minus;
end
end
