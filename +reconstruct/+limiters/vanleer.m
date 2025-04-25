function limited_slope = vanleer(delta_minus, delta_plus)
% VANLEER Van Leer slope limiter for TVD schemes
%
% The Van Leer limiter is a popular choice, offering a balance between
% accuracy and oscillation control. It's smoother than minmod.
% Author: Dr. Denys Dutykh
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

% Check if slopes have opposite signs
if delta_minus * delta_plus <= 0
    limited_slope = 0;
else
    % Apply Van Leer limiter formula
    % phi(r) = (r + |r|) / (1 + |r|) = 2r / (1+r) for r>0, 0 otherwise
    % Slope = phi(r) * delta_plus
    
    % Denominator is delta_minus + delta_plus
    denominator = delta_minus + delta_plus;
    
    % Avoid division by zero if both deltas are zero (already handled by sign check)
    if denominator == 0 
        limited_slope = 0;
    else
        % Calculate the harmonic mean type expression
        limited_slope = 2 * delta_minus * delta_plus / denominator;
    end
end

end
