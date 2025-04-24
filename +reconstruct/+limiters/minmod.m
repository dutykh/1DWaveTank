function limited_slope = minmod(delta_minus, delta_plus)

% MINMOD Minmod slope limiter for TVD schemes
%
% The minmod limiter returns zero if the slopes have different signs
% and returns the minimum magnitude slope if they have the same sign.
% It is the most dissipative TVD limiter.
%
% Syntax:
%   limited_slope = minmod(delta_minus, delta_plus)
%
% Inputs:
%   delta_minus - Backward difference (q_i - q_{i-1})
%   delta_plus  - Forward difference (q_{i+1} - q_i)
%
% Output:
%   limited_slope - Limited slope

% Check signs of the slope estimates
if delta_minus * delta_plus <= 0
    % Opposite signs or one is zero - return zero
    limited_slope = 0;
else
    % Same sign - return the smaller magnitude slope
    limited_slope = sign(delta_minus) * min(abs(delta_minus), abs(delta_plus));
end

end